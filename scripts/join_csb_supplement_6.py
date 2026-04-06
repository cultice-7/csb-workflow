import geopandas as gpd
import pandas as pd
import rasterio
import os
from pathlib import Path
from rasterstats import zonal_stats

# Import parameters from Snakemake
CSB_input_folder = snakemake.params.CSB_input_dir
CSB_output_folder = snakemake.params.CSB_output_dir
weather_input_folder = snakemake.params.weather_input_dir
states = snakemake.params.states
CSB_years = snakemake.params.CSB_years
weather_variables = snakemake.params.weather_variables


for year in CSB_years:
    for state in states:
        
        CSB_input_path = os.path.join(CSB_input_folder, f"{state}_CSB{year}_CSBID_geometry.parquet")
        output_path_spatial = os.path.join(CSB_output_folder, f"{state}_CSB{year}_supplement_6_spatial.parquet")
        output_path_table = os.path.join(CSB_output_folder, f"{state}_CSB{year}_supplement_6_table.parquet")

        # Load CSB_shape datasets
        CSB_shape = gpd.read_parquet(CSB_input_path)
        
        # Keep only the columns necessary for the spatial join
        cols_to_keep = ['CSBID', 'geometry']
        CSB_shape = CSB_shape[cols_to_keep]
        
        # Extract parcel centroids
        centroids = CSB_shape.geometry.centroid
        centroid_coords = [(pt.x, pt.y) for pt in centroids]
        
        for variable in weather_variables:
            input_dir = Path(weather_input_folder) / f"{variable}"
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
        
        # Convert float64 columns to float32 to save memory
        float64_cols = CSB_shape.select_dtypes(include=["float64"]).columns
        CSB_shape[float64_cols] = CSB_shape[float64_cols].astype("float32")

        #---# Save files w/ and w/o geometry
        #CSB_shape.to_parquet(output_path_geojson, compression="zstd")
        attribute_table = CSB_shape.drop(columns='geometry')
        attribute_table.to_parquet(output_path_table, index=False, compression="zstd")