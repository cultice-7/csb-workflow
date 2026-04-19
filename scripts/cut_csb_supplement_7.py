import geopandas as gpd
import pandas as pd
import numpy as np
import os
from pathlib import Path
from sklearn.neighbors import KDTree
import shutil

# Import parameters from Snakemake
CSB_input_folder = snakemake.params.CSB_input_dir
CSB_output_folder = snakemake.params.CSB_output_dir
states = snakemake.params.states
CSB_years = snakemake.params.CSB_years
drop_pattern = ["elevator_nearest"]


def cut_crop_price_columns(state, year, drop_pattern, CSB_input_folder, CSB_output_folder):
    
    # Path to input and output paths
    csb_supplement_7_input_path = os.path.join(CSB_input_folder, f"{state}_CSB{year}_supplement_7_table.parquet")
    csb_supplement_7_output_path = os.path.join(CSB_output_folder, f"{state}_CSB{year}_supplement_7_table_reduced.parquet")
    
    # Load Regrow datasets
    csb_supplement_7 = pd.read_parquet(csb_supplement_7_input_path)
    print(state, year, csb_supplement_7.shape)
    
    cols_to_drop = csb_supplement_7.filter(regex="|".join(drop_pattern)).columns
    csb_supplement_7 = csb_supplement_7.drop(columns=cols_to_drop)
    
    # Convert float64 columns to float32 to save memory
    float64_cols = csb_supplement_7.select_dtypes(include=["float64"]).columns
    csb_supplement_7[float64_cols] = csb_supplement_7[float64_cols].astype("float32")
    
    print(state, year, csb_supplement_7.shape)
    csb_supplement_7.to_parquet(csb_supplement_7_output_path, index=False, compression="zstd")
    
    
# Main code
for CSB_year in CSB_years:
    for state in states:

        cut_crop_price_columns(state, CSB_year, drop_pattern, CSB_input_folder, CSB_output_folder)