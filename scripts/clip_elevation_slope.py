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
elevation_input_folder = snakemake.params.elevation_input_dir
slope_input_folder = snakemake.params.slope_input_dir
elevation_output_folder = snakemake.params.elevation_output_dir
slope_output_folder = snakemake.params.slope_output_dir
states = snakemake.params.states
target_CRS = snakemake.params.target_CRS

# Read file with state boundaries
state_bound = gpd.read_file(Path(state_bound_folder) / f"cb_2023_us_state_500k.shp")

# Path to inpurt elevation and slope raster file
elevation_input_file_path =  os.path.join(elevation_input_folder, "elevation_reprojected.tif")
slope_input_file_path =  os.path.join(slope_input_folder, "slope_reprojected.tif")
        
for state in states:
    
    # Extract state boundaries
    select_state_bound = state_bound[state_bound['STUSPS'] == state]

    # Reproject to the the target CRS
    select_state_bound = select_state_bound.to_crs(target_CRS)
    
    # 3. Create 30 km outward buffer
    state_bound_outer = select_state_bound.buffer(30000)
    
    # Path for elevation output files
    elevation_output = os.path.join(elevation_output_folder, f"{state}_elevation_clipped.tif")
    Path(elevation_output).parent.mkdir(parents=True, exist_ok=True)
    
    # Apply function to clip elevation
    clip_raster(
        src_tif=elevation_input_file_path,
        state_gdf = state_bound_outer,
        output_path=elevation_output)
    
    print(f"Elevation for {state} is successfully clipped")
    
    # Path for slope output files
    slope_output = os.path.join(slope_output_folder, f"{state}_slope_clipped.tif")
    Path(elevation_output).parent.mkdir(parents=True, exist_ok=True)
    
    # Apply function to clip slope
    clip_raster(
        src_tif=slope_input_file_path,
        state_gdf = state_bound_outer,
        output_path=slope_output)
    
    print(f"Slope for {state} is successfully clipped")