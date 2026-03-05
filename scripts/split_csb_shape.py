import geopandas as gpd
import pandas as pd
import os
import gc

# Paths to input and output directories
input_folder = "data/edited/CSB/"
output_folder = "data/edited/CSB/"

states = snakemake.params.states
years = snakemake.params.years

for year in years:
    for state in states:
        
        # Input and outpur file names
        input_path = os.path.join(input_folder, f"{state}_CSB{year}_clipped.gpkg")
        output_path_joined = os.path.join(output_folder, f"{state}_CSB{year}_shape_table.parquet")
        output_path_fieldID_geometry_parquet = os.path.join(output_folder, f"{state}_CSB{year}_CSBID_geometry.parquet")
        output_path_fieldID_geometry_gpkg = os.path.join(output_folder, f"{state}_CSB{year}_CSBID_geometry.gpkg")
        output_path_table = os.path.join(output_folder, f"{state}_CSB{year}_table.parquet")
        
        # Read input file
        gdf = gpd.read_file(input_path)
        gdf['CSBID'] = gdf['CSBID'].astype("string")

        # Reproject to NAD83/Conus Albers
        gdf_reproj = gdf.to_crs(epsg=5070)

        print(gdf_reproj.shape)

        # Save regrow shape_table file
        gdf_reproj = gdf_reproj.set_crs(epsg=5070)
        gdf_reproj.to_parquet(output_path_joined, compression="zstd")
        # Keep only field_id and geomtery to create a separate file with field_id as a "key" and geometry only
        gdf_fieldID_geomtery = gdf_reproj[['CSBID', 'geometry']].copy()
        gdf_fieldID_geomtery.to_parquet(output_path_fieldID_geometry_parquet, compression="zstd")
        gdf_fieldID_geomtery.to_file(output_path_fieldID_geometry_gpkg, driver = "GPKG")
        # Create a file without geometry
        attribute_table = gdf_reproj.drop(columns='geometry')
        attribute_table.to_parquet(output_path_table, compression="zstd")
        
        print(f'CSB{year} data files for {state} are created and saved')
        
        # Delete dataframes from memory
        del gdf_reproj, gdf_fieldID_geomtery, attribute_table
        gc.collect()