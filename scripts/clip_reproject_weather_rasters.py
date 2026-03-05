import geopandas as gpd
import rasterio
from rasterio.mask import mask
from rasterio.warp import reproject, calculate_default_transform, Resampling
from pathlib import Path

# Function that clip and reproject each raster file
def reproject_clip_raster(
    src_tif,
    states_gdf,
    out_tif_temp,
    out_tif,
    dst_crs
):
    with rasterio.open(src_tif) as src:

        # Calculate transform for reprojection to CRS 5070
        transform, width, height = calculate_default_transform(
            src.crs,
            dst_crs,
            src.width,
            src.height,
            *src.bounds
        )

        raster_meta = src.meta.copy()
        raster_meta.update({
            "crs": dst_crs,
            "transform": transform,
            "width": width,
            "height": height
        })
            
        with rasterio.open(out_tif_temp, 'w', **raster_meta) as dst:
            for i in range(1, src.count + 1):
                reproject(
                    source=rasterio.band(src, i),
                    destination=rasterio.band(dst, i),
                    src_transform=src.transform,
                    src_crs=src.crs,
                    dst_transform=transform,
                    dst_crs=dst_crs,
                    resampling=Resampling.average)

    with rasterio.open(out_tif_temp) as src_5070:
        out_image, out_transform = mask(src_5070, states_gdf.geometry, crop=True)
        
        out_meta = src_5070.meta.copy()
        out_meta.update({
            "driver": "GTiff",
            "height": out_image.shape[1],
            "width": out_image.shape[2],
            "transform": out_transform
        })
        
        with rasterio.open(out_tif, "w", **out_meta) as final_dst:
            final_dst.write(out_image)
    
    out_tif_temp.unlink()
    

state_bound = gpd.read_file("data/Census/state_bound/cb_2023_us_state_500k.shp")
states = snakemake.params.states
weather_variables = snakemake.params.weather_variables
    
for variable in weather_variables:
    input_dir = Path(f"data/Weather/{variable}")
    weather_files = sorted(input_dir.glob(f"prism_{variable}_us_30s_*.tif"))
    
    for file in weather_files:
        
        for state in states:
            
            # Extract state boundaries
            select_states_bound = state_bound[state_bound['STUSPS'] == state]

            # Reproject to the common CRS (5070)
            select_states_bound = select_states_bound.to_crs(epsg=5070)
            
            # Exctract date from PRISM file name
            date = file.stem.split("_")[-1]  # YYYYMM
            
            # Path for output files
            weather_raster_temp = Path(f"data/edited/Weather/{variable}/{state}_prism_{variable}_us_30s_{date}_temp.tif")
            weather_raster_output = Path(f"data/edited/Weather/{variable}/{state}_prism_{variable}_us_30s_{date}_clipped.tif")
            Path(weather_raster_output).parent.mkdir(parents=True, exist_ok=True)
            
            # Apply function to reproject and clip PRISM raster file
            reproject_clip_raster(
                src_tif=file,
                states_gdf = select_states_bound,
                out_tif_temp=weather_raster_temp,
                out_tif=weather_raster_output,
                dst_crs="EPSG:5070"
            )
            print(f"Raster for {state} and {file} is successfully reprojected and clipped")