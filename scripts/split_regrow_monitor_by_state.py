import pandas as pd
import os
import numpy as np

# Import parameters from Snakemake
input_folder = snakemake.params.input_dir
output_folder = snakemake.params.output_dir
states = snakemake.params.states
states_monitor = snakemake.params.states_monitor


def split_monitor_data_by_state(monitor_merged, state, input_folder, output_folder):
    
    regrow_shape_input_path = os.path.join(input_folder, f"{state}_field_boundaries_2014-2025.parquet")
    regrow_monitor_output_path = os.path.join(output_folder, f"{state}_Monitor_data_2014-2025.parquet")
    
    # Read regrow geometry file and keep only "boundary_id" column
    regrow_shape = pd.read_parquet(regrow_shape_input_path)
    cols_to_keep = ['boundary_id']
    regrow_shape = regrow_shape[cols_to_keep] 
    
    # Split monitor data by state
    regrow_monitor_filtered = monitor_merged[monitor_merged["boundary_id"].isin(regrow_shape["boundary_id"])]
    regrow_monitor_filtered.reset_index(drop=True, inplace=True)
    
    # Save the split monitor data
    os.makedirs(output_folder, exist_ok = True)
    regrow_monitor_filtered.to_parquet(regrow_monitor_output_path, index = False, compression="zstd")


# Main script
monitor_merged = pd.DataFrame()

for state in states_monitor:
    monitor_input_path = os.path.join(input_folder, f"Monitor_data_{state}_2014-2025.parquet")
    monitor = pd.read_parquet(monitor_input_path)
    if monitor_merged.empty:
        monitor_merged = monitor.copy()
    else:
        monitor_merged = pd.concat([monitor_merged, monitor], ignore_index=True)
        
for state in states:
    split_monitor_data_by_state(monitor_merged, state, input_folder, output_folder)