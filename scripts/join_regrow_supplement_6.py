import geopandas as gpd
import pandas as pd
import rasterio
import os
from pathlib import Path
from rasterstats import zonal_stats


# Input and output folders for Regrow
input_folder_Regrow = "data/edited/Regrow/"
output_folder_Regrow = "data/edited/Regrow/"

states = snakemake.params.states
weather_variables = snakemake.params.weather_variables

for state in states:
    
    input_path_Regrow = os.path.join(input_folder_Regrow, f"{state}_regrow_shape_table.geojson")
    output_path_geojson = os.path.join(output_folder_Regrow, f"{state}_regrow_supplement_6_spatial.geojson")
    output_path_csv = os.path.join(output_folder_Regrow, f"{state}_regrow_supplement_6_table.csv")
    
    # Load regrow_dises joined datasets
    regrow_shape = gpd.read_file(input_path_Regrow)
    
    # Reproject geometry to the same CRS (NAD83/CONUS Albers)
    regrow_shape = regrow_shape.to_crs(epsg=5070)
    
    # Keep only the columns necessary for the spatial join
    cols_to_keep = ['field_id', 'geometry']
    regrow_shape = regrow_shape[cols_to_keep]
    
    # Extract parcel centroids
    centroids = regrow_shape.geometry.centroid
    centroid_coords = [(pt.x, pt.y) for pt in centroids]
    
    for variable in weather_variables:
        input_dir = Path(f"data/edited/Weather/{variable}")
        weather_files = sorted(input_dir.glob(f"{state}_prism_{variable}_us_30s_*.tif"))
        new_cols = {}
        
        for file in weather_files:
            date = file.stem.split("_")[-2]  # YYYYMM
            print(f"Adding weather variable {variable} for {state} and {date}")
            
            # Zonal statistics for weather variables
            #weather_stats = zonal_stats(regrow_shape, file, stats="mean")
            #regrow_shape[f'{variable}_mean_{date}'] = [stat['mean'] for stat in weather_stats]
            
            with rasterio.open(file) as src:
                weather_stats = list(src.sample(centroid_coords))
            
            new_cols[f'{variable}_mean_{date}'] = [stat[0] for stat in weather_stats]
            
        regrow_shape = regrow_shape.join(pd.DataFrame(new_cols, index=regrow_shape.index))

    #---# Save geojson and csv files
    #regrow_shape.to_file(output_path_geojson, driver="GeoJSON")
    attribute_table = regrow_shape.drop(columns='geometry')
    attribute_table.to_csv(output_path_csv, index=False)