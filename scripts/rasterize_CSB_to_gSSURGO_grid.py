import numpy as np
import pandas as pd
import geopandas as gpd
import rasterio
from rasterio.features import rasterize
import os
from pathlib import Path


# Input and output folders for CSB
input_folder_CSB = snakemake.params.CSB_input_dir
output_folder_CSB = snakemake.params.CSB_output_dir

# Pull list of states for running the code
states = snakemake.params.states

def CSB_duplicating_fields_rasterization(state, input_path_CSB, output_path_raster, output_path_duplicating_fields):
    
    # Load CSB_dises joined datasets
    CSB_geometry = gpd.read_parquet(input_path_CSB)
    
    # Reproject geometry to the same CRS (NAD83/CONUS Albers)
    CSB_geometry = CSB_geometry.to_crs(epsg=5070)
    
    
    
    ### Clean CSB: CSB geometries contain overlaps which harms rasterization, we need to remove overlaps ###
    gdf = CSB_geometry.copy()
    gdf["geometry"] = gdf["geometry"].buffer(0)
    gdf["area"] = gdf.geometry.area

    # Overlay CSBID on itself. This will produce all overlapping polygons (including boundaries)
    CSB_overlaps = gpd.overlay(gdf, gdf, how="intersection")

    # Remove self-overlaps
    CSB_overlaps = CSB_overlaps[CSB_overlaps["CSBID_1"] != CSB_overlaps["CSBID_2"]]

    # Compute overlap ratio
    # Determine the smaller area in each pair
    CSB_overlaps["min_area"] = CSB_overlaps[["area_1", "area_2"]].min(axis=1)

    # Compute overlap fraction relative to smaller polygon
    CSB_overlaps["overlap_fraction"] = CSB_overlaps.geometry.area / CSB_overlaps["min_area"]

    # Keep only those with overlap >= 50%
    CSB_overlaps = CSB_overlaps[CSB_overlaps["overlap_fraction"] >= 0.5]

    # Keep only unique pairs (A,B) same as (B,A)
    # Make a tuple of (smaller_area_parcel, larger_area_parcel)
    CSB_overlaps["pair"] = CSB_overlaps.apply(
        lambda row: (row["CSBID_1"], row["CSBID_2"]) if row["area_1"] <= row["area_2"] else (row["CSBID_2"], row["CSBID_1"]), axis=1)
    # Drop duplicates
    CSB_overlaps = CSB_overlaps.drop_duplicates(subset="pair").reset_index(drop=True)
    # Assign smaller field IDs in each pair to CSBID_1 column and larger field IDs to CSBID_2
    CSB_overlaps.loc[:, 'CSBID_1'] = CSB_overlaps["pair"].apply(lambda x: x[0] if pd.notna(x) else None)
    CSB_overlaps.loc[:, 'CSBID_2'] = CSB_overlaps["pair"].apply(lambda x: x[1] if pd.notna(x) else None)
    # Rename CSBID_1 and CSBID_2 accordingly
    CSB_overlaps.rename(columns={"CSBID_1": "CSBID_smaller", "CSBID_2": "CSBID_larger"}, inplace = True)
    # Keep only necessary columns
    CSB_overlaps = CSB_overlaps[['CSBID_smaller', 'CSBID_larger', 'overlap_fraction']]
    # Save a separate dataset with overlapping fields
    CSB_overlaps.to_parquet(output_path_duplicating_fields, compression="zstd")
    
    # Clean CSB fields to remove overlapping fields
    # Extract all smaller parcel IDs from the pairs
    overlap_smaller_CSBIDs = CSB_overlaps["CSBID_smaller"].unique()
    # Drop rows whose CSBID is in smaller_ids and reset index
    CSB_geometry = (CSB_geometry[~CSB_geometry["CSBID"].isin(overlap_smaller_CSBIDs)]).reset_index(drop=True)
    
    
    
    ### Process mukey raster file to create pairs of CSBID CSBID: mukeys (pixel values) ###
    # CSBID CSBID → unique field integers
    id_map = {s: i for i, s in enumerate(CSB_geometry["CSBID"].unique())}
    # Unique field integers → CSBID CSBID
    reverse_id_map = {v: k for k, v in id_map.items()} 
    CSB_geometry["pid"] = CSB_geometry["CSBID"].map(id_map)

    # Rasterize CSB
    with rasterio.open(f"data/edited/Soil/gSSURGO Mukey Grid/{state}_MURASTER_30m.tif") as src:
        # Match CRS between vector CSBID and gSSURGO raster
        CSB_geometry = CSB_geometry.to_crs(src.crs)

        # Read gSSURGO raster file
        band = src.read(1)
        transform = src.transform
        height = src.height
        width = src.width

        # Rasterize fields: assign unique field integers to pixels whose centers lie inside a given polygon
        shapes = ((geom, pid) for geom, pid in zip(CSB_geometry.geometry, CSB_geometry.pid))
        
        parcel_raster = rasterize(
            shapes=shapes,
            out_shape=(height, width),
            transform=transform,
            fill=-1,            # pixels not belonging to any parcel
            all_touched=False,  # only pixels whose center is inside polygon
            dtype="int32",
            skip_invalid = False
        )

        # Save raster to disk
        out_meta = {
            "driver": "GTiff",
            "height": height,
            "width": width,
            "count": 1,
            "dtype": parcel_raster.dtype,
            "crs": src.crs,
            "transform": transform
        }

        with rasterio.open(output_path_raster, "w", **out_meta) as dst:
            dst.write(parcel_raster, 1)


for state in states:
    
    input_path_CSB = os.path.join(input_folder_CSB, f"{state}_CSB1724_CSBID_geometry.parquet")
    output_path_raster = os.path.join(output_folder_CSB, f"{state}_CSB1724_raster_to_gSSURGO_grid.tif")
    output_path_duplicating_fields = os.path.join(output_folder_CSB, f"{state}_CSB1724_duplicating_fields.parquet")
    
    CSB_duplicating_fields_rasterization(state, input_path_CSB, output_path_raster, output_path_duplicating_fields)