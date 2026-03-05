import geopandas as gpd
from rasterstats import zonal_stats
import os

#---# Load required datasets
# Paths to the reprojected slope and elevation raster files
elevation_proj_path = "data/Geo/elevation/elevation_reproj.tif"
slope_proj_path = "data/Geo/elevation/slope_reproj.tif"

# Input and output folders for Regrow
input_folder_CSB = "data/edited/CSB/"
output_folder_CSB = "data/edited/CSB/"

# List of states
states = snakemake.params.states


for state in states:
    
    input_path_CSB = os.path.join(input_folder_CSB, f"{state}_CSB1724_clipped.gpkg")

    # Load CSB data
    CSB_shape = gpd.read_file(input_path_CSB)

    # Setting active geometry column
    CSB_shape = CSB_shape.set_geometry('geometry')
    
    # Reproject geometry to an equal-area CRS (NAD83/CONUS Albers)
    CSB_shape = CSB_shape.to_crs(epsg=5070)
    
    # Keep only the columns necessary for the spatial join
    cols_to_keep = ['CSBID', 'geometry']
    CSB_shape = CSB_shape[cols_to_keep]
    
    
    #---# Add zonal statistics for elevation and slope
    # Set paths to output files
    output_path_geojson = os.path.join(output_folder_CSB, f"{state}_CSB1724_supplement_2_spatial.geojson")
    output_path_csv = os.path.join(output_folder_CSB, f"{state}_CSB1724_supplement_2_table.csv")
    
    # Zonal statistics for elevation
    try:
        print("Calculating mean elevation...")
        elevation_stats = zonal_stats(CSB_shape, elevation_proj_path, stats="mean")
        elevation_means = [stat['mean'] for stat in elevation_stats]
        CSB_shape['elevation_mean'] = elevation_means
        print("Mean elevation added to attribute table.")
    except Exception as e:
        print(f"Error processing elevation: {e}")
        raise

    # Zonal statistics for slope
    try:
        print("Calculating mean slope...")
        slope_stats = zonal_stats(CSB_shape, slope_proj_path, stats="mean")
        slope_means = [stat['mean'] for stat in slope_stats]
        CSB_shape['slope_mean'] = slope_means
        print("Mean slope added to attribute table.")
    except Exception as e:
        print(f"Error processing slope: {e}")
    
    # Check whether there are any missing values
    cols_to_check = ['CSBID', 'elevation_mean', 'slope_mean']
    if CSB_shape[cols_to_check].isna().any().any():
        print(CSB_shape[CSB_shape[cols_to_check].isna().any(axis=1)])
    
    # Save geojson and csv files
    #CSB_shape.to_file(output_path_geojson, driver="GeoJSON")
    attribute_table = CSB_shape.drop(columns='geometry')
    print(attribute_table.shape) #check df shape
    attribute_table.to_csv(output_path_csv, index=False)
    print(f'Supplementary data 2 for {state} is created and saved')
        