import geopandas as gpd
import rasterio
from rasterio.mask import mask
from pathlib import Path

# Function that clip and reproject each raster file
def clip_raster(
    src_tif,
    state_gdf,
    output_path):
    
    with rasterio.open(src_tif) as src:

        # Convert geometry to list format expected by rasterio
        geom = [state_gdf.iloc[0]]

        # Crop raster
        out_image, out_transform = mask(
            src,
            geom,
            crop=True
        )

        out_meta = src.meta.copy()
        
        # Update metadata
        out_meta.update({
            "height": out_image.shape[1],
            "width": out_image.shape[2],
            "transform": out_transform
        })
        
        # Save cropped raster
        with rasterio.open(output_path, "w", **out_meta) as dest:
            dest.write(out_image)


state_bound = gpd.read_file("data/Census/state_bound/cb_2023_us_state_500k.shp")
states = snakemake.params.states
input_file_path = Path("data/gSSURGO/FY2026_gSSURGO_mukey_grid\MURASTER_30m.tif")
        
for state in states:
    
    # Extract state boundaries
    select_state_bound = state_bound[state_bound['STUSPS'] == state]

    # Reproject to the common CRS (5070)
    select_state_bound = select_state_bound.to_crs(epsg=5070)
    
    # 3. Create 10 km outward buffer
    state_bound_outer = select_state_bound.buffer(10000)
    
    # Path for output files
    mukey_raster_output = Path(f"data/edited/Soil/gSSURGO Mukey Grid/{state}_MURASTER_30m.tif")
    Path(mukey_raster_output).parent.mkdir(parents=True, exist_ok=True)
    
    # Apply function to clip gSSURGO mukey raster file
    clip_raster(
        src_tif=input_file_path,
        state_gdf = state_bound_outer,
        output_path=mukey_raster_output)
    
    print(f"Mukey raster for {state} is successfully clipped")