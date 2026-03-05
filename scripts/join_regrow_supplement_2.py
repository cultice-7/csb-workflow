import geopandas as gpd
from rasterstats import zonal_stats
import os

#---# Load required datasets
# Paths to the reprojected slope and elevation raster files
elevation_proj_path = "data/Geo/elevation/elevation_reproj.tif"
slope_proj_path = "data/Geo/elevation/slope_reproj.tif"

# Input and output folders for Regrow
input_folder_Regrow = "data/edited/Regrow/"
output_folder_Regrow = "data/edited/Regrow/"

# List of states
states = snakemake.params.states


for state in states:
    
    input_path_Regrow = os.path.join(input_folder_Regrow, f"{state}_regrow_shape_table.geojson")
    
    # Load Regrow data
    regrow_shape = gpd.read_file(input_path_Regrow)

    # Setting active geometry column
    regrow_shape = regrow_shape.set_geometry('geometry')
    
    # Reproject geometry to an equal-area CRS (NAD83/CONUS Albers)
    regrow_shape = regrow_shape.to_crs(epsg=5070)
    
    # Keep only the columns necessary for the spatial join
    cols_to_keep = ['field_id', 'geometry']
    regrow_shape = regrow_shape[cols_to_keep]
    
    
    #---# Add zonal statistics for elevation and slope
    # Set paths to output files
    output_path_geojson = os.path.join(output_folder_Regrow, f"{state}_regrow_supplement_2_spatial.geojson")
    output_path_csv = os.path.join(output_folder_Regrow, f"{state}_regrow_supplement_2_table.csv")
    
    # Zonal statistics for elevation
    try:
        print("Calculating and adding mean elevation...")
        elevation_stats = zonal_stats(regrow_shape, elevation_proj_path, stats="mean")
        elevation_means = [stat['mean'] for stat in elevation_stats]
        regrow_shape['elevation_mean'] = elevation_means
        print("Mean elevation added to attribute table.")
    except Exception as e:
        print(f"Error processing elevation: {e}")
        raise

    # Zonal statistics for slope
    try:
        print("Calculating and adding mean slope...")
        slope_stats = zonal_stats(regrow_shape, slope_proj_path, stats="mean")
        slope_means = [stat['mean'] for stat in slope_stats]
        regrow_shape['slope_mean'] = slope_means
        print("Mean slope added to attribute table.")
    except Exception as e:
        print(f"Error processing slope: {e}")
        raise
    
    # Check whether there are any missing values
    cols_to_check = ['field_id', 'elevation_mean', 'slope_mean']
    if regrow_shape[cols_to_check].isna().any().any():
        print(regrow_shape[regrow_shape[cols_to_check].isna().any(axis=1)])
    
    # Save geojson and csv files
    #regrow_shape.to_file(output_path_geojson, driver="GeoJSON")
    attribute_table = regrow_shape.drop(columns='geometry')
    print(attribute_table.shape) #check df shape
    attribute_table.to_csv(output_path_csv, index=False)
    print(f'Supplementary data 2 for {state} is created and saved')
    

