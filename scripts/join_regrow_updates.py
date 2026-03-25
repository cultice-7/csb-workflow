import geopandas as gpd
import pandas as pd
import os
import gc

# Paths to input and output directories
input_regrow_2014_2024 = snakemake.params.input_2014_2024_dir
input_regrow_2025 = snakemake.params.input_2025_dir
output_joined_regrow = snakemake.params.output_dir

# List of states
states = snakemake.params.states

# Running vertical concatenation for raw regrow geometry datasets
for state in states:
    
    # Input and output file names
    input_geometry_2014_2024_path = os.path.join(input_regrow_2014_2024, f"{state}_field_boundaries.geojson")
    input_geometry_2025_path = os.path.join(input_regrow_2025, f"{state}_2025_boundaries.geojson")
    output_geometry_merged = os.path.join(output_joined_regrow, f"{state}_field_boundaries.geojson")
    
    # Read input file
    geometry_2014_2024 = gpd.read_file(input_geometry_2014_2024_path)
    geometry_2025 = gpd.read_file(input_geometry_2025_path)
    
    # Merge two geometry datasets for 2014-2024 and 2025
    geometry_merged = geometry_2014_2024.merge(geometry_2025, on="boundary_id", how="outer", suffixes=("_1", "_2"))
    # Prefer geometry from geometry_2014_2024, fallback to geometry_2025
    geometry_merged["geometry"] = geometry_merged["geometry_1"].combine_first(geometry_merged["geometry_2"])
    geometry_merged = geometry_merged[["boundary_id", "geometry"]]
    
    # Save as a geojson file
    geometry_merged.to_file(output_geometry_merged, driver="GeoJSON")

    # Clean memory
    del geometry_2014_2024, geometry_2025, geometry_merged
    gc.collect()


# Monitor data is represented in data files combining multiple states
states_monitor = ["OH", "MN_WI_IA_IN_IL", "MI"]
# Running vertical concatenation for raw regrow monitor datasets
for state_monitor in states_monitor:
    
    # Input and output file names
    input_table_2014_2024_path = os.path.join(input_regrow_2014_2024, f"Monitor_data_{state_monitor}.csv")
    output_table_concatenated = os.path.join(output_joined_regrow, f"Monitor_data_{state_monitor}.csv")
    # Read input file
    table_2014_2024 = pd.read_csv(input_table_2014_2024_path)
    
    table_concatenated = table_2014_2024.copy()
    
    for state in state_monitor.split("_"):
        
        # Input file names
        input_table_2025_path = os.path.join(input_regrow_2025, f"{state}_2025_Data.csv")
        # Read input file
        table_2025 = pd.read_csv(input_table_2025_path)
        
        # Vertically concatenate two monitor datasets for 2014-2024 and 2025
        table_concatenated = pd.concat([table_concatenated, table_2025], ignore_index=True)
    
    # Save as a csv file
    table_concatenated.to_csv(output_table_concatenated, index=False)        
    
    # Clean memory
    del table_2014_2024, table_2025, table_concatenated
    gc.collect()