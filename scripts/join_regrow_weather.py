import geopandas as gpd
import pandas as pd
import numpy as np
import rasterio
import os
from pathlib import Path
from rasterstats import zonal_stats

# Import parameters from Snakemake
regrow_input_folder = snakemake.params.regrow_input_dir
regrow_output_folder = snakemake.params.regrow_output_dir
weather_input_folder = snakemake.params.weather_input_dir
states = snakemake.params.states
weather_variables = snakemake.params.weather_variables


for state in states:
    
    regrow_input_path = os.path.join(regrow_input_folder, f"{state}_regrow_fieldID_geometry.parquet")
    output_path_spatial = os.path.join(regrow_output_folder, f"{state}_regrow_weather_spatial.parquet")
    output_path_table = os.path.join(regrow_output_folder, f"{state}_regrow_weather_table.parquet")
    
    # Load regrow_dises joined datasets
    regrow_shape = gpd.read_parquet(regrow_input_path)
    
    # Keep only the columns necessary for the spatial join
    cols_to_keep = ['field_id', 'geometry']
    regrow_shape = regrow_shape[cols_to_keep]
    
    # Extract parcel centroids
    centroids = regrow_shape.geometry.centroid
    centroid_coords = [(pt.x, pt.y) for pt in centroids]
    
    for variable in weather_variables:
        input_dir = Path(weather_input_folder) / f"{variable}"
        weather_files = sorted(input_dir.glob(f"{state}_prism_{variable}_us_30s_*.tif"))
        
        # Collect all new columns for this variable once
        new_cols = {}
        
        for file in weather_files:
            date = file.stem.split("_")[-2]  # YYYYMM
            print(f"Adding weather variable {variable} for {state} and {date}")
            
            # Zonal statistics for weather variables
            #weather_stats = zonal_stats(regrow_shape, file, stats="mean")
            #regrow_shape[f'{variable}_mean_{date}'] = [stat['mean'] for stat in weather_stats]
            
            with rasterio.open(file) as src:
                weather_stats = np.fromiter(
                    (v[0] for v in src.sample(centroid_coords)),
                    dtype=np.float32,
                    count=len(centroid_coords)
                )
            
            new_cols[f'{variable}_centroid_{date}'] = weather_stats
            
        regrow_shape = regrow_shape.join(pd.DataFrame(new_cols, index=regrow_shape.index))

    # Convert float64 columns to float32 to save memory
    float64_cols = regrow_shape.select_dtypes(include=["float64"]).columns
    regrow_shape[float64_cols] = regrow_shape[float64_cols].astype("float32")
    
    #---# Save files w/ and w/o geometry
    #regrow_shape.to_parquet(output_path_spatial, compression="zstd")
    attribute_table = regrow_shape.drop(columns='geometry')
    attribute_table.to_parquet(output_path_table, index = False, compression="zstd")