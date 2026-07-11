import geopandas as gpd
import pandas as pd
import numpy as np
import os
from pathlib import Path
from sklearn.neighbors import KDTree
import shutil

# Import parameters from Snakemake
regrow_input_folder = snakemake.params.regrow_input_dir
regrow_output_folder = snakemake.params.regrow_output_dir
states = snakemake.params.states
drop_pattern = ["elevator_nearest"]


def cut_crop_price_columns(state, drop_pattern, regrow_input_folder, regrow_output_folder):
    
    # Path to input and output paths
    regrow_supplement_7_input_path = os.path.join(regrow_input_folder, f"{state}_regrow_crop_prices_table.parquet")
    regrow_supplement_7_output_path = os.path.join(regrow_output_folder, f"{state}_regrow_crop_prices_table_reduced.parquet")
    
    # Load Regrow datasets
    regrow_supplement_7 = pd.read_parquet(regrow_supplement_7_input_path)
    print(state, regrow_supplement_7.shape)
    
    cols_to_drop = regrow_supplement_7.filter(regex="|".join(drop_pattern)).columns
    regrow_supplement_7 = regrow_supplement_7.drop(columns=cols_to_drop)
    
    # Convert float64 columns to float32 to save memory
    float64_cols = regrow_supplement_7.select_dtypes(include=["float64"]).columns
    regrow_supplement_7[float64_cols] = regrow_supplement_7[float64_cols].astype("float32")
    
    print(state, regrow_supplement_7.shape)
    regrow_supplement_7.to_parquet(regrow_supplement_7_output_path, index=False, compression="zstd")
    
    
# Main code
for state in states:   

    cut_crop_price_columns(state, drop_pattern, regrow_input_folder, regrow_output_folder)