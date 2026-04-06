import pandas as pd
import os
import numpy as np

# Import parameters from Snakemake
input_folder = snakemake.params.input_dir
output_folder = snakemake.params.output_dir
states = snakemake.params.states
states_monitor = snakemake.params.states_monitor


def split_monitor_data_by_state(monitor_merged, state, input_folder, output_folder):
    
    regrow_shape_input_path = os.path.join(input_folder, f"{state}_field_boundaries.parquet")
    regrow_monitor_output_path = os.path.join(output_folder, f"{state}_Monitor_data_cleaned.parquet")
    
    # Read regrow geometry file and keep only "boundary_id" column
    regrow_shape = pd.read_parquet(regrow_shape_input_path)
    cols_to_keep = ['boundary_id']
    regrow_shape = regrow_shape[cols_to_keep] 
    
    # Merge regrow geometry with monitor (tabular) data
    regrow_shape_monitor = regrow_shape.merge(monitor_merged, on="boundary_id", how="left")
    
    # Drop all rows with completely missing data
    regrow_shape_monitor = regrow_shape_monitor.dropna(subset=['cycle_end', 'crop'], how='all')
    
    # Save the split monitor data
    os.makedirs(output_folder, exist_ok = True)
    regrow_shape_monitor.to_parquet(regrow_monitor_output_path, compression="zstd")


# Main script

for state in states_monitor:
    monitor_input_path = os.path.join(input_folder, f"Monitor_data_{state}.parquet")
    monitor = pd.read_parquet(monitor_input_path)
    try:
        monitor_merged
        monitor_merged = pd.concat([monitor_merged, monitor], ignore_index=True)
    except NameError:
        monitor_merged = monitor.copy()

for state in states:
    split_monitor_data_by_state(monitor_merged, state, input_folder, output_folder)