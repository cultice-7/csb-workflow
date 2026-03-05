import geopandas as gpd
import pandas as pd
import rasterio
import os
from pathlib import Path
from rasterstats import zonal_stats


# Input and output folders for CSB and road
input_folder_CSB = "data/edited/CSB/"
output_folder_CSB = "data/edited/CSB/"

states = snakemake.params.states
weather_variables = snakemake.params.weather_variables

for state in states:
    
    input_path_CSB = os.path.join(input_folder_CSB, f"{state}_CSB1724_clipped.gpkg")
    output_path_geojson = os.path.join(output_folder_CSB, f"{state}_CSB1724_supplement_6_spatial.geojson")
    output_path_csv = os.path.join(output_folder_CSB, f"{state}_CSB1724_supplement_6_table.csv")

    # Load CSB_shape datasets
    CSB_shape = gpd.read_file(input_path_CSB)
    
    # Reproject geometry to the same CRS (NAD83/CONUS Albers)
    CSB_shape = CSB_shape.to_crs(epsg=5070)
    
    # Keep only the columns necessary for the spatial join
    cols_to_keep = ['CSBID', 'geometry']
    CSB_shape = CSB_shape[cols_to_keep]
    
    # Extract parcel centroids
    centroids = CSB_shape.geometry.centroid
    centroid_coords = [(pt.x, pt.y) for pt in centroids]
    
    for variable in weather_variables:
        input_dir = Path(f"data/edited/Weather/{variable}")
        weather_files = sorted(input_dir.glob(f"{state}_prism_{variable}_us_30s_*.tif"))
        new_cols = {}
        
        for file in weather_files:
            date = file.stem.split("_")[-2]  # YYYYMM
            print(f"Adding weather variable {variable} for {state} and {date}")
            
            # Zonal statistics for weather variables
            #weather_stats = zonal_stats(CSB_shape, file, stats="mean")
            #CSB_shape[f'{variable}_mean_{date}'] = [stat['mean'] for stat in weather_stats]
            
            with rasterio.open(file) as src:
                weather_stats = list(src.sample(centroid_coords))
            
            new_cols[f'{variable}_mean_{date}'] = [stat[0] for stat in weather_stats]
            
        CSB_shape = CSB_shape.join(pd.DataFrame(new_cols, index=CSB_shape.index))

    #---# Save geojson and csv files
    #CSB_shape.to_file(output_path_geojson, driver="GeoJSON")
    attribute_table = CSB_shape.drop(columns='geometry')
    attribute_table.to_csv(output_path_csv, index=False)