import geopandas as gpd
import pandas as pd
import os

# Import parameters from the Snakemake
states = snakemake.params.states
target_CRS = snakemake.params.target_CRS
regrow_geometry_input_folder = snakemake.params.input_geometry_dir
regrow_output_folder = snakemake.params.output_dir


def split_clean_save_regrow_geometry(state, target_CRS, regrow_geometry_input_folder, regrow_output_folder):
    
    # Paths to input files
    input_geometry_path = os.path.join(regrow_geometry_input_folder, f"{state}_field_boundaries_2014-2025.parquet")
    
    # Paths to output files
    regrow_geometry_output_path_parquet = os.path.join(regrow_output_folder, f"{state}_regrow_fieldID_geometry.parquet")
    regrow_geometry_output_path_gpkg = os.path.join(regrow_output_folder, f"{state}_regrow_fieldID_geometry.gpkg")
    
    # Upload Regrow raw geometry file
    regrow_shape = gpd.read_parquet(input_geometry_path)

    # Reproject field geometry to NAD83/Conus Albers
    regrow_shape_reproj = regrow_shape.to_crs(target_CRS)
    
    # Buffer geometry to address invalid geometries (by reconstructing them)
    regrow_shape_reproj["geometry"] = regrow_shape_reproj["geometry"].buffer(0)

    # Rename the geometry column "boundary_id" to "field_id"
    regrow_shape_reproj.rename(columns={'boundary_id': 'field_id'}, inplace = True)
    
    # Keep only 'field_id' and 'geometry' to create a separate file where 'field_id' serves as the key
    gdf_fieldID_geomtery = regrow_shape_reproj[['field_id', 'geometry']].copy()
    # Save the geometry‑only output in both Parquet and GPKG formats
    gdf_fieldID_geomtery.to_parquet(regrow_geometry_output_path_parquet, compression="zstd")
    gdf_fieldID_geomtery.to_file(regrow_geometry_output_path_gpkg, driver = "GPKG")
    
    print(f'Regrow geometry file for {state} is created and saved')


#---# Main script
for state in states:
    split_clean_save_regrow_geometry(state, target_CRS, regrow_geometry_input_folder, regrow_output_folder)