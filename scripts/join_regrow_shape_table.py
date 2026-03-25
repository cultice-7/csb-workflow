import geopandas as gpd
import pandas as pd
import os


def join_regrow_shape_table(state, input_geometry_path, input_table_path, output_folder):
    
    # Outpur file paths
    output_path_joined = os.path.join(output_folder, f"{state}_regrow_shape_table.parquet")
    output_path_fieldID_geometry_parquet = os.path.join(output_folder, f"{state}_regrow_fieldID_geometry.parquet")
    output_path_fieldID_geometry_gpkg = os.path.join(output_folder, f"{state}_regrow_fieldID_geometry.gpkg")
    output_path_table = os.path.join(output_folder, f"{state}_regrow_table.parquet")

    # Upload Regrow tabular data (attributes)
    regrow_table = pd.read_parquet(input_table_path)
    
    # Read input geometry file
    gdf = gpd.read_file(input_geometry_path)

    # Reproject field geometry to NAD83/Conus Albers
    gdf_reproj = gdf.to_crs(epsg=5070)

    # Rename the geometry column "boundary_id" to "field_id"
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
    
    # Keep only 'field_id' and 'geometry' to create a separate file where 'field_id' serves as the key
    gdf_fieldID_geomtery = gdf_joined[['field_id', 'geometry']].copy()
    # Save the geometry‑only output in both Parquet and GPKG formats
    gdf_fieldID_geomtery = gdf_fieldID_geomtery.set_crs(epsg=5070)
    gdf_fieldID_geomtery.to_parquet(output_path_fieldID_geometry_parquet, compression="zstd")
    gdf_fieldID_geomtery.to_file(output_path_fieldID_geometry_gpkg, driver = "GPKG")
    
    # Create a separate file with tabular data (excluding geometry)
    attribute_table = gdf_joined.drop(columns='geometry')
    attribute_table.to_parquet(output_path_table, compression="zstd")
    
    print(f'Regrow data files for {state} are created and saved')


#---# Execute the join operation between the Regrow shapefile and the associated table files
# List of states
states = snakemake.params.states

# Paths to input and output directories
input_geometry_folder = snakemake.params.input_geometry_dir
output_folder = snakemake.params.output_dir
input_table_folder = snakemake.params.input_table_dir

for state in states:

    # Input directories
    input_geometry_path = os.path.join(input_geometry_folder, f"{state}_field_boundaries.geojson")
    input_table_path = os.path.join(input_table_folder, "OH_MN_WI_IA_IN_IL_MI_regrow_wide_coded.parquet")
    
    join_regrow_shape_table(state, input_geometry_path, input_table_path, output_folder)