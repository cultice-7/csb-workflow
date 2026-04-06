import geopandas as gpd
import pandas as pd
import os
import gc

# Import parameters from Snakemake
regrow_2014_2024_input_folder = snakemake.params.input_2014_2024_dir
regrow_2025_input_folder = snakemake.params.input_2025_dir
regrow_joined_output_folder = snakemake.params.output_dir
states = snakemake.params.states
states_monitor = snakemake.params.states_monitor


# Merge regrow 2014-2024 and 2025 geometries
def megre_geometry(state, regrow_2014_2024_input_folder, regrow_2025_input_folder, regrow_joined_output_folder):
    
    # Input and output file names
    input_geometry_2014_2024_path = os.path.join(regrow_2014_2024_input_folder, f"{state}_field_boundaries.geojson")
    input_geometry_2025_path = os.path.join(regrow_2025_input_folder, f"{state}_2025_boundaries.geojson")
    output_geometry_merged = os.path.join(regrow_joined_output_folder, f"{state}_field_boundaries.parquet")
    
    # Read input file
    geometry_2014_2024 = gpd.read_file(input_geometry_2014_2024_path)
    geometry_2025 = gpd.read_file(input_geometry_2025_path)
    
    # Merge two geometry datasets for 2014-2024 and 2025
    geometry_merged = geometry_2014_2024.merge(geometry_2025, on="boundary_id", how="outer", suffixes=("_1", "_2"))
    # Prefer geometry from geometry_2014_2024, fallback to geometry_2025
    geometry_merged["geometry"] = geometry_merged["geometry_1"].combine_first(geometry_merged["geometry_2"])
    geometry_merged = geometry_merged[["boundary_id", "geometry"]]
    
    # Save as a geojson file
    geometry_merged.to_parquet(output_geometry_merged, compression="zstd")


# Running vertical concatenation for raw regrow monitor datasets
def concat_monitor_data(state_monitor, regrow_2014_2024_input_folder, regrow_2025_input_folder, regrow_joined_output_folder):
    
    # Input and output file names
    input_table_2014_2024_path = os.path.join(regrow_2014_2024_input_folder, f"Monitor_data_{state_monitor}.csv")
    output_table_concatenated = os.path.join(regrow_joined_output_folder, f"Monitor_data_{state_monitor}.parquet")
    # Read input file
    table_2014_2024 = pd.read_csv(input_table_2014_2024_path)
    
    table_concatenated = table_2014_2024.copy()
    
    for state in state_monitor.split("_"):
        
        # Input file names
        input_table_2025_path = os.path.join(regrow_2025_input_folder, f"{state}_2025_Data.csv")
        # Read input file
        table_2025 = pd.read_csv(input_table_2025_path)
        
        # Vertically concatenate two monitor datasets for 2014-2024 and 2025
        table_concatenated = pd.concat([table_concatenated, table_2025], ignore_index=True)
    
    # Save as a csv file
    table_concatenated.to_parquet(output_table_concatenated, index=False, compression="zstd")


# Main script
for state in states:
    megre_geometry(state, regrow_2014_2024_input_folder, regrow_2025_input_folder, regrow_joined_output_folder)
    
for state_monitor in states_monitor:
    concat_monitor_data(state_monitor, regrow_2014_2024_input_folder, regrow_2025_input_folder, regrow_joined_output_folder)