import geopandas as gpd
import pandas as pd
import numpy as np
import os
from pathlib import Path
import rasterio
from rasterio.features import rasterize
import gc
  
# Import parameters from the Snakemake
regrow_input_folder = snakemake.params.regrow_input_dir
regrow_raster_input_folder = snakemake.params.regrow_raster_input_dir
regrow_checks_folder = snakemake.params.regrow_checks_dir
regrow_validation_output_folder = snakemake.params.regrow_validation_dir
CDL_input_folder = snakemake.params.cdl_input_dir
states = snakemake.params.states
years_range = snakemake.params.years


def mode_crop(x):
    x = x.dropna()

    if x.empty:
        return np.nan

    return np.bincount(x.astype(int)).argmax()


def merge_validate_regrow_cdl(state, year, regrow_table, regrow_cdl, regrow_raster_input_folder, regrow_checks_folder, CDL_input_folder):
    
    # Paths to input and output files
    regrow_raster_input_path = os.path.join(regrow_raster_input_folder, f"{state}_regrow_raster_to_CDL_grid.tif")
    cdl_input_path = os.path.join(CDL_input_folder, f"{state}_{year}_30m_cdls_clipped.tif")
    fieldID_pid_input_path = os.path.join(regrow_raster_input_folder, f"{state}_regrow_fieldID_pid_correspondence.parquet")
    overlapping_fields_input_path = os.path.join(regrow_checks_folder, f"{state}_regrow_overlapping_fields.parquet")
    
    
    
    ### Process CDL raster file to create pairs of Regrow field_id: CDL crop values ###
    # Upload the dataset with a mapping from Regrow field_id to unique field integers
    id_map = pd.read_parquet(fieldID_pid_input_path)
    id_map = id_map.set_index("pid")

    # Open CDL raster
    with rasterio.open(cdl_input_path) as src_CDL:
        # Read gSSURGO mukey raster
        CDL_raster = src_CDL.read(1)

    # Open rasterized Regrow
    with rasterio.open(regrow_raster_input_path) as src_regrow:
        # Read rasterized Regrow
        regrow_raster = src_regrow.read(1)

    # Extract pixel values for each field
    # Indices of all pixels that belong to some field
    rows, cols = np.where(regrow_raster >= 0)
    field_pids = regrow_raster[rows, cols]
    cdl_crops = CDL_raster[rows, cols]
    
    # Remove raster files from the memory
    del CDL_raster, regrow_raster
    gc.collect()

    # Assign CDL crops for each Regrow field. Compute the mode of assigned CDL crop values for each Regrow field
    fieldID_CDL_crop = (pd.DataFrame({
    "field_id": id_map.loc[field_pids, "field_id"].values,
    f"cdl_crop_{year}": cdl_crops})
    .groupby("field_id")[f"cdl_crop_{year}"]   
    .agg(mode_crop))
    
    regrow_cdl = regrow_cdl.merge(fieldID_CDL_crop, on='field_id', how="left")
    
    # Map each smaller_field_id to its corresponding larger_field_id and then copy the values from the larger fields back into the target dataset
    regrow_overlaps = pd.read_parquet(overlapping_fields_input_path)
    regrow_cdl = regrow_cdl.set_index("field_id")
    regrow_cdl.loc[regrow_overlaps["field_id_smaller"]] = regrow_cdl.loc[regrow_overlaps["field_id_larger"]].values
    regrow_cdl = regrow_cdl.reset_index()
    
    
    
    ### Compare Regrow and CDL main crop values for each Regrow field ###
    # Select the main crop columns in Regrow for a given year
    crop_cols = regrow_table.filter(regex=f"crop_{str(year)[2:]}_*").columns.tolist()
    
    regrow_table = regrow_table.set_index("field_id")
    regrow_cdl = regrow_cdl.merge(regrow_table[crop_cols], on='field_id', how='left')
    
    regrow_cdl[f"cdl_valid_{year}"] = regrow_cdl[crop_cols].eq(regrow_cdl[f"cdl_crop_{year}"], axis=0).any(axis=1).astype(int)
    
    return regrow_cdl


# Main code
# Create an output directory
os.makedirs(Path(regrow_validation_output_folder), exist_ok = True)

for state in states:
    # Paths to input and output files
    regrow_table_input_path = os.path.join(regrow_input_folder, f"{state}_regrow_table.parquet")
    regrow_cdl_output_path = os.path.join(regrow_validation_output_folder, f"{state}_regrow_cdl_validation_table.parquet")
    
    # Read Regrow input files
    regrow_table = pd.read_parquet(regrow_table_input_path)
    
    # Create regrow_cdl dataset
    regrow_cdl = regrow_table[['field_id']].copy()

    for year in years_range:
        # Path to CDL raster file
        raster_path = os.path.join(CDL_input_folder, f"{state}_{year}_30m_cdls_clipped.tif")

        print(f"Processing {state} and year {year}...")
        
        regrow_cdl = merge_validate_regrow_cdl(state, year, regrow_table, regrow_cdl, regrow_raster_input_folder, regrow_checks_folder, CDL_input_folder)
    
    # Keep only necessary columns
    regrow_cdl = regrow_cdl.loc[:, regrow_cdl.columns.str.contains("field_id|cdl", case=False)]
    
    ### Save output files ###
    # Convert crop-code columns to nullable Int16 and validity flag columns to plain int16
    cdl_crop_cols = regrow_cdl.filter(regex="^cdl_crop_").columns
    regrow_cdl[cdl_crop_cols] = regrow_cdl[cdl_crop_cols].astype("Int16")

    cdl_valid_cols = regrow_cdl.filter(regex="^cdl_valid_").columns
    regrow_cdl[cdl_valid_cols] = regrow_cdl[cdl_valid_cols].astype("int16")

    # Convert all float64 to float32
    float64_cols = regrow_cdl.select_dtypes(include="float64").columns
    regrow_cdl[float64_cols] = regrow_cdl[float64_cols].astype("float32")

    # Save the output of Regrow joined with the CDL dataset
    regrow_cdl.to_parquet(regrow_cdl_output_path, compression="zstd")
    print(f"Joining Regrow and CDL for {state} is complete, output file is saved")
    
