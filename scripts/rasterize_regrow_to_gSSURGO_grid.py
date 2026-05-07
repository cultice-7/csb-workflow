import numpy as np
import pandas as pd
import geopandas as gpd
import rasterio
from rasterio.features import rasterize
import os
from pathlib import Path


# Import parameters from Snakemake
regrow_input_folder = snakemake.params.regrow_input_dir
regrow_checks_folder = snakemake.params.regrow_checks_dir
soil_input_folder = snakemake.params.soil_input_dir
regrow_raster_output_folder = snakemake.params.regrow_raster_output_dir
states = snakemake.params.states
target_CRS = snakemake.params.target_CRS


def regrow_rasterization(state, regrow_input_folder, regrow_checks_folder, soil_input_folder, regrow_raster_output_folder):
    
    # Path to input and output fields
    regrow_input_path = os.path.join(regrow_input_folder, f"{state}_regrow_fieldID_geometry.parquet")
    soil_input_path = os.path.join(soil_input_folder, f"gSSURGO Mukey Grid/{state}_MURASTER_30m.tif")
    overlapping_fields_input_path = os.path.join(regrow_checks_folder, f"{state}_regrow_overlapping_fields.parquet")
    regrow_output_path = os.path.join(regrow_raster_output_folder, f"{state}_regrow_raster_to_gSSURGO_grid.tif")
    fieldID_pid_output_path= os.path.join(regrow_raster_output_folder, f"{state}_regrow_fieldID_pid_correspondence.parquet")
    
    
    # Load regrow shape (geometry) data
    regrow_geometry = gpd.read_parquet(regrow_input_path)
    
    ### Clean Regrow: Regrow geometries can overlap which harms rasterization, we need to remove overlaps ###
    regrow_overlaps = pd.read_parquet(overlapping_fields_input_path)
    # Extract all smaller parcel IDs from the pairs
    overlap_smaller_field_ids = regrow_overlaps["field_id_smaller"].unique()
    # Drop rows whose field_id is in smaller_ids and reset index
    regrow_geometry = (regrow_geometry[~regrow_geometry["field_id"].isin(overlap_smaller_field_ids)]).reset_index(drop=True)
  
    # In rasterio, if multiple polygons overlap a single pixel, the value from the last polygon in the input sequence is assigned to that pixel in the final raster
    # If polygons overlap slightly, we want overlapping areas to take the value of the larger polygon
    # Sort polygons by area in ascending order so that larger polygons are processed last
    regrow_geometry = regrow_geometry.iloc[regrow_geometry.geometry.area.argsort()]
    
    
    
    ### Process mukey raster file to create mapping mukey (pixel value): Regrow unique identifier ###
    # Implement mapping Regrow field_id → unique field integers (pid)
    id_map = {s: i for i, s in enumerate(regrow_geometry["field_id"].unique())}
    regrow_geometry["pid"] = regrow_geometry["field_id"].map(id_map)
    # Save Regrow-pid mapping
    regrow_geometry[['field_id', 'pid']].to_parquet(fieldID_pid_output_path, compression="zstd")

    # Rasterize Regrow based on mukey raster file
    with rasterio.open(soil_input_path) as src:
        # Match CRS between vector Regrow and gSSURGO raster
        regrow_geometry = regrow_geometry.to_crs(src.crs)

        # Properties of gSSURGO raster file
        transform = src.transform
        height = src.height
        width = src.width

        # Rasterize fields: assign unique field integers (pid) to pixels whose centers lie inside a given polygon
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

        with rasterio.open(regrow_output_path, "w", **out_meta) as dst:
            dst.write(parcel_raster, 1)
            
        print(f"Rasterization for {state} is complete and output files are saved.")


# Main code
# Make sure "Rasterization to gSSURGO grid" folder exists
Path(regrow_raster_output_folder).mkdir(parents=True, exist_ok=True)

for state in states:
    regrow_rasterization(state, regrow_input_folder, regrow_checks_folder, soil_input_folder, regrow_raster_output_folder)