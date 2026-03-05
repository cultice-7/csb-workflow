import geopandas as gpd
import pandas as pd
import os
import gc

# Paths to input and output directories
input_folder = "data/Regrow/"
output_folder = "data/edited/Regrow/"

# Upload Regrow table data (attributes)
regrow_table = pd.read_parquet("data/edited/Regrow/OH_MN_WI_IA_IN_IL_MI_regrow_wide_coded.parquet")

states = snakemake.params.states

for state in states:
    
    # Input and outpur file names
    input_path = os.path.join(input_folder, f"{state}_field_boundaries.geojson")
    output_path_joined = os.path.join(output_folder, f"{state}_regrow_shape_table.parquet")
    output_path_fieldID_geometry_parquet = os.path.join(output_folder, f"{state}_regrow_fieldID_geometry.parquet")
    output_path_fieldID_geometry_gpkg = os.path.join(output_folder, f"{state}_regrow_fieldID_geometry.gpkg")
    output_path_table = os.path.join(output_folder, f"{state}_regrow_table.parquet")
    
    # Read input file
    gdf = gpd.read_file(input_path)

    # Reproject to NAD83/Conus Albers
    gdf_reproj = gdf.to_crs(epsg=5070)

    # Rename the column "boundary_id"
    gdf_reproj.rename(columns={'boundary_id': 'field_id'}, inplace = True)

    # Calculate area in US survey acres (1 acre = 4046.8564224 sq m)
    gdf_reproj["area_acre"] = gdf_reproj.geometry.area / 4046.8564224

    # Join shape file and table using boundary ID
    gdf_joined = gdf_reproj.merge(regrow_table, on = 'field_id', how = 'left')

    # Replace NaN with a placeholder string 
    #gdf_joined = gdf_joined.fillna('NA')

    print(gdf_joined.shape)

    # Save regrow shape_table file
    gdf_joined = gdf_joined.set_crs(epsg=5070)
    gdf_joined.to_parquet(output_path_joined, compression="zstd")
    # Keep only field_id and geomtery to create a separate file with field_id as a "key" and geometry only
    gdf_fieldID_geomtery = gdf_joined[['field_id', 'geometry']].copy()
    gdf_fieldID_geomtery.to_parquet(output_path_fieldID_geometry_parquet, compression="zstd")
    gdf_fieldID_geomtery.to_file(output_path_fieldID_geometry_gpkg, driver = "GPKG")
    # Create a file without geometry
    attribute_table = gdf_joined.drop(columns='geometry')
    attribute_table.to_parquet(output_path_table, compression="zstd")
    
    print(f'Regrow data files for {state} are created and saved')
    
    # Delete dataframes from memory
    del gdf_joined, gdf_fieldID_geomtery, attribute_table
    gc.collect()