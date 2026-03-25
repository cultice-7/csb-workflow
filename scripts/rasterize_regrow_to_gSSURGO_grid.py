import numpy as np
import pandas as pd
import geopandas as gpd
import rasterio
from rasterio.features import rasterize
import os
from pathlib import Path
import pickle



# Input and output folders for Regrow
input_folder_Regrow = snakemake.params.regrow_input_dir
output_folder_Regrow = snakemake.params.regrow_output_dir

# Pull list of states for running the code
states = snakemake.params.states

def regrow_duplicating_fields_rasterization(state, input_path_Regrow, output_path_raster, output_path_duplicating_fields):
    
    # Load regrow_dises joined datasets
    regrow_geometry = gpd.read_parquet(input_path_Regrow)
    
    # Reproject geometry to the same CRS (NAD83/CONUS Albers)
    regrow_geometry = regrow_geometry.to_crs(epsg=5070)
    
    
    
    ### Clean Regrow: Regrow geometries contain overlaps which harms rasterization, we need to remove overlaps ###
    gdf = regrow_geometry.copy()
    gdf["geometry"] = gdf["geometry"].buffer(0)
    gdf["area"] = gdf.geometry.area

    # Overlay Regrow on itself. This will produce all overlapping polygons (including boundaries)
    regrow_overlaps = gpd.overlay(gdf, gdf, how="intersection")

    # Remove self-overlaps
    regrow_overlaps = regrow_overlaps[regrow_overlaps["field_id_1"] != regrow_overlaps["field_id_2"]]

    # Compute overlap ratio
    # Determine the smaller area in each pair
    regrow_overlaps["min_area"] = regrow_overlaps[["area_1", "area_2"]].min(axis=1)

    # Compute overlap fraction relative to smaller polygon
    regrow_overlaps["overlap_fraction"] = regrow_overlaps.geometry.area / regrow_overlaps["min_area"]

    # Keep only those with overlap >= 50%
    regrow_overlaps = regrow_overlaps[regrow_overlaps["overlap_fraction"] >= 0.5]

    # Keep only unique pairs (A,B) same as (B,A)
    # Make a tuple of (smaller_area_parcel, larger_area_parcel)
    regrow_overlaps["pair"] = regrow_overlaps.apply(
        lambda row: (row["field_id_1"], row["field_id_2"]) if row["area_1"] <= row["area_2"] else (row["field_id_2"], row["field_id_1"]), axis=1)
    # Drop duplicates
    regrow_overlaps = regrow_overlaps.drop_duplicates(subset="pair").reset_index(drop=True)
    # Assign smaller field IDs in each pair to field_id_1 column and larger field IDs to field_id_2
    regrow_overlaps.loc[:, 'field_id_1'] = regrow_overlaps["pair"].apply(lambda x: x[0] if pd.notna(x) else None)
    regrow_overlaps.loc[:, 'field_id_2'] = regrow_overlaps["pair"].apply(lambda x: x[1] if pd.notna(x) else None)
    # Rename field_id_1 and field_id_2 accordingly
    regrow_overlaps.rename(columns={"field_id_1": "field_id_smaller", "field_id_2": "field_id_larger"}, inplace = True)
    # Keep only necessary columns
    regrow_overlaps = regrow_overlaps[['field_id_smaller', 'field_id_larger', 'overlap_fraction']]
    # Save a separate dataset with overlapping fields
    regrow_overlaps.to_parquet(output_path_duplicating_fields, compression="zstd")
    
    # Clean regrow fields to remove overlapping fields
    # Extract all smaller parcel IDs from the pairs
    overlap_smaller_field_ids = regrow_overlaps["field_id_smaller"].unique()
    # Drop rows whose field_id is in smaller_ids and reset index
    regrow_geometry = (regrow_geometry[~regrow_geometry["field_id"].isin(overlap_smaller_field_ids)]).reset_index(drop=True)
    
    
    
    ### Process mukey raster file to create pairs of Regrow field_id: mukeys (pixel values) ###
    # Regrow field_id → unique field integers
    id_map = {s: i for i, s in enumerate(regrow_geometry["field_id"].unique())}
    # Unique field integers → Regrow field_id
    reverse_id_map = {v: k for k, v in id_map.items()} 
    regrow_geometry["pid"] = regrow_geometry["field_id"].map(id_map)

    # Rasterize Regrow based on mukey raster file
    with rasterio.open(f"data/edited/Soil/gSSURGO Mukey Grid/{state}_MURASTER_30m.tif") as src:
        # Match CRS between vector Regrow and gSSURGO raster
        regrow_geometry = regrow_geometry.to_crs(src.crs)

        # Read gSSURGO raster file
        mukey_raster = src.read(1)
        transform = src.transform
        height = src.height
        width = src.width

        # Rasterize fields: assign unique field integers to pixels whose centers lie inside a given polygon
        shapes = ((geom, pid) for geom, pid in zip(regrow_geometry.geometry, regrow_geometry.pid))
        
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
            
        print(f"Rasterization for {state} is complete and output files are saved.")

for state in states:
    
    input_path_Regrow = os.path.join(input_folder_Regrow, f"{state}_regrow_fieldID_geometry.parquet")
    output_path_raster = os.path.join(output_folder_Regrow, f"{state}_regrow_raster_to_gSSURGO_grid.tif")
    output_path_duplicating_fields = os.path.join(output_folder_Regrow, f"{state}_regrow_duplicating_fields.parquet")
    
    regrow_duplicating_fields_rasterization(state, input_path_Regrow, output_path_raster, output_path_duplicating_fields)