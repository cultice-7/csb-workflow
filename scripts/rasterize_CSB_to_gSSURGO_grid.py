import numpy as np
import pandas as pd
import geopandas as gpd
import rasterio
from rasterio.features import rasterize
import os
from pathlib import Path


# Import parameters from Snakemake
CSB_input_folder = snakemake.params.CSB_input_dir
CSB_checks_folder = snakemake.params.CSB_checks_dir
soil_input_folder = snakemake.params.soil_input_dir
CSB_raster_output_folder = snakemake.params.CSB_raster_output_dir
states = snakemake.params.states
CSB_years = snakemake.params.CSB_years
target_CRS = snakemake.params.target_CRS


def CSB_rasterization(state, CSB_year, CSB_input_folder, CSB_checks_folder, soil_input_folder, CSB_raster_output_folder):
    
    CSB_input_path = os.path.join(CSB_input_folder, f"{state}_CSB{CSB_year}_CSBID_geometry.parquet")
    soil_input_path = os.path.join(soil_input_folder, f"gSSURGO Mukey Grid/{state}_MURASTER_30m.tif")
    overlapping_fields_input_path = os.path.join(CSB_checks_folder, f"{state}_CSB{CSB_year}_overlapping_fields.parquet")
    CSB_output_path = os.path.join(CSB_raster_output_folder, f"{state}_CSB{CSB_year}_raster_to_gSSURGO_grid.tif")
    CSBID_pid_output_path = os.path.join(CSB_raster_output_folder, f"{state}_CSB{CSB_year}_CSBID_pid_correspondence.parquet")
    
    
    # Load CSB_dises joined datasets
    CSB_geometry = gpd.read_parquet(CSB_input_path)


    ### Clean CSB: CSB geometries contain overlaps which harms rasterization, we need to remove those overlaps ###
    CSB_overplaps = pd.read_parquet(overlapping_fields_input_path)
    # Extract all smaller parcel IDs from the pairs
    overlap_smaller_CSBIDs = CSB_overplaps["CSBID_smaller"].unique()
    # Drop rows whose CSBID is in smaller_ids and reset index
    CSB_geometry = (CSB_geometry[~CSB_geometry["CSBID"].isin(overlap_smaller_CSBIDs)]).reset_index(drop=True)
    
    # In rasterio, if multiple polygons overlap a single pixel, the value from the last polygon in the input sequence is assigned to that pixel in the final raster
    # If polygons overlap slightly, we want overlapping areas to take the value of the larger polygon
    # Sort polygons by area in ascending order so that larger polygons are processed last
    CSB_geometry = CSB_geometry.iloc[CSB_geometry.geometry.area.argsort()]
    
    
    
    ### Process mukey raster file to create mapping mukey (pixel value): CSB unique identifier ###
    # Implement mapping CSBID → unique field integers
    id_map = {s: i for i, s in enumerate(CSB_geometry["CSBID"].unique())}
    CSB_geometry["pid"] = CSB_geometry["CSBID"].map(id_map)
    CSB_geometry[['CSBID', 'pid']].to_parquet(CSBID_pid_output_path, compression="zstd")

    # Rasterize CSB
    with rasterio.open(soil_input_path) as src:
        # Match CRS between vector CSBID and gSSURGO raster
        CSB_geometry = CSB_geometry.to_crs(src.crs)

        # Properties of gSSURGO raster file
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

        with rasterio.open(CSB_output_path, "w", **out_meta) as dst:
            dst.write(parcel_raster, 1)


# Main code
# Make sure "Rasterization to gSSURGO grid" folder exists
Path(CSB_raster_output_folder).mkdir(parents=True, exist_ok=True)

for CSB_year in CSB_years:
    for state in states:
        CSB_rasterization(state, CSB_year, CSB_input_folder, CSB_checks_folder, soil_input_folder, CSB_raster_output_folder)