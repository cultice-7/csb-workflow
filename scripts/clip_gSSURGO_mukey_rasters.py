import geopandas as gpd
import rasterio
from rasterio.mask import mask
from pathlib import Path
import os

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

# Main code
# Import parameters from Snakemake
state_bound_folder = snakemake.params.state_bound_dir
soil_input_folder = snakemake.params.soil_input_dir
soil_output_folder = snakemake.params.soil_output_dir
states = snakemake.params.states
target_CRS = snakemake.params.target_CRS
state_buffer_m = snakemake.params.state_buffer_margin

# Read file with state boundaries
state_bound = gpd.read_file(Path(state_bound_folder) / f"cb_2023_us_state_500k.shp")

# Path to inpurt mukey raster file
mukey_input_file_path =  os.path.join(soil_input_folder, "FY2026_gSSURGO_mukey_grid/MURASTER_30m.tif")
        
for state in states:
    
    # Extract state boundaries
    select_state_bound = state_bound[state_bound['STUSPS'] == state]

    # Reproject to the the target CRS
    select_state_bound = select_state_bound.to_crs(target_CRS)
    
    # 3. Create 10 km outward buffer
    state_bound_outer = select_state_bound.buffer(state_buffer_m)
    
    # Path for output files
    mukey_raster_output = os.path.join(soil_output_folder, f"gSSURGO Mukey Grid/{state}_MURASTER_30m.tif")
    Path(mukey_raster_output).parent.mkdir(parents=True, exist_ok=True)
    
    # Apply function to clip gSSURGO mukey raster file
    clip_raster(
        src_tif=mukey_input_file_path,
        state_gdf = state_bound_outer,
        output_path=mukey_raster_output)
    
    print(f"Mukey raster for {state} is successfully clipped")