import geopandas as gpd
import pandas as pd
import os

# Import parameters from the Snakemake
states = snakemake.params.states
target_CRS = snakemake.params.target_CRS
regrow_geometry_input_folder = snakemake.params.input_geometry_dir
regrow_table_input_folder = snakemake.params.input_table_dir
regrow_output_folder = snakemake.params.output_dir


def join_regrow_shape_table(state, target_CRS, regrow_geometry_input_folder, regrow_table_input_folder, regrow_output_folder):
    
    # Paths to input files
    input_geometry_path = os.path.join(regrow_geometry_input_folder, f"{state}_regrow_fieldID_geometry.parquet")
    input_table_path = os.path.join(regrow_table_input_folder, f"{state}_regrow_monitor_wide_coded.parquet")
    
    # Paths to output files
    regrow_shape_table_output_path = os.path.join(regrow_output_folder, f"{state}_regrow_shape_table.parquet")
    regrow_table_output_path = os.path.join(regrow_output_folder, f"{state}_regrow_table.parquet")

    # Upload Regrow tabular data (attributes)
    regrow_table = pd.read_parquet(input_table_path)
    
    # Upload Regrow shape (geometry) file
    regrow_shape = gpd.read_parquet(input_geometry_path)

    # Calculate area in US survey acres (1 acre = 4046.8564224 sq m)
    regrow_shape["area_acre"] = regrow_shape.geometry.area / 4046.8564224

    # Join shape file and table using boundary ID
    regrow_shape_table_joined = regrow_shape.merge(regrow_table, on = 'field_id', how = 'left')

    # Replace NaN with a placeholder string 
    #gdf_joined = gdf_joined.fillna('NA')

    print(regrow_shape_table_joined.shape)

    # Save regrow shape_table file
    regrow_shape_table_joined = regrow_shape_table_joined.set_crs(target_CRS)
    regrow_shape_table_joined.to_parquet(regrow_shape_table_output_path, compression="zstd")
    
    # Create a separate file with tabular data (excluding geometry)
    attribute_table = regrow_shape_table_joined.drop(columns='geometry')
    attribute_table.to_parquet(regrow_table_output_path, compression="zstd")
    
    print(f'Regrow Monitor data files for {state} are created and saved')


#---# Main script
for state in states:
    join_regrow_shape_table(state, target_CRS, regrow_geometry_input_folder, regrow_table_input_folder, regrow_output_folder)