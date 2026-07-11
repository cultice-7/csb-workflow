import geopandas as gpd
import rasterio
from rasterio.vrt import WarpedVRT
from rasterio.warp import Resampling
from rasterio.mask import mask
from pathlib import Path
import os
import numpy as np

# Function that clip and reproject each raster file
def reproject_clip_raster(
    src_tif,
    states_gdf,
    out_tif,
    dst_crs
):

    with rasterio.open(src_tif) as src:
        # Virtual reprojection (no full raster in memory)
        with WarpedVRT(
            src,
            crs=dst_crs,
            src_nodata=-9999,
            nodata=np.nan,
            resampling=Resampling.average
        ) as vrt:

            # Mask reads ONLY needed pixels
            out_image, out_transform = mask(
                vrt,
                states_gdf.geometry,
                crop=True
            )

            out_meta = vrt.meta.copy()
            out_meta.update({
                "driver": "GTiff",
                "height": out_image.shape[1],
                "width": out_image.shape[2],
                "transform": out_transform
            })

            with rasterio.open(out_tif, "w", **out_meta) as dst:
                dst.write(out_image)

# Main code
# Import parameters from Snakemake
state_bound_folder = snakemake.params.state_bound_dir
weather_input_folder = snakemake.params.weather_input_dir
weather_output_folder = snakemake.params.weather_output_dir
states = snakemake.params.states
weather_variables = snakemake.params.weather_variables
target_CRS = snakemake.params.target_CRS

# Read file with state boundaries
state_bound = gpd.read_file(Path(state_bound_folder) / f"cb_2023_us_state_500k.shp")
# Reproject to the the target CRS
state_bound = state_bound.to_crs(target_CRS)
    
# Loop over all weather variables to reproject and clip to state boundaries
for variable in weather_variables:
    input_dir = Path(weather_input_folder) / f"{variable}"
    weather_files = sorted(input_dir.glob(f"prism_{variable}_us_30s_*.tif"))
    
    for file in weather_files:
        
        for state in states:
            
            # Extract state boundaries
            select_states_bound = state_bound[state_bound['STUSPS'] == state]
            
            # Exctract date from PRISM file name
            date = file.stem.split("_")[-1]  # YYYYMM
            
            # Path for output files
            weather_output_path = os.path.join(weather_output_folder, f"{variable}/{state}_prism_{variable}_us_30s_{date}_clipped.tif")
            Path(weather_output_path).parent.mkdir(parents=True, exist_ok=True)
            
            # Apply function to reproject and clip PRISM weather file to state boundaries
            reproject_clip_raster(
                src_tif=file,
                states_gdf = select_states_bound,
                out_tif=weather_output_path,
                dst_crs=target_CRS
            )
            
            print(f"Weather raster for {state} and {file} is successfully reprojected and clipped")