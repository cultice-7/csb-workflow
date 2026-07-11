import geopandas as gpd
import os
import rasterio
from rasterio.mask import mask
from rasterio.windows import from_bounds
from pathlib import Path

# Import parameters from the Snakemake
state_bound_folder = snakemake.params.state_bound_dir
states = snakemake.params.states
years_range = snakemake.params.years
target_CRS = snakemake.params.target_CRS
cdl_input_folder = snakemake.params.cdl_input_dir
cdl_output_folder = snakemake.params.cdl_output_dir
state_buffer_m = snakemake.params.state_buffer_margin


def clip_cdl_using_state_boundaries(state_bound_outer, state, year, cdl_input_folder, cdl_output_folder):
    
    cdl_path = os.path.join(cdl_input_folder, f"{year}_30m_cdls", f"{year}_30m_cdls.tif")
    output_path = os.path.join(cdl_output_folder, f"{state}_{year}_30m_cdls_clipped.tif")

    print(f"Processing {state} and year {year}...")

    try:
        with rasterio.open(cdl_path) as src:
            out_image, out_transform = mask(src, state_bound_outer.geometry, crop=True)
            out_meta = src.meta.copy()
            out_meta.update({
                "driver": "GTiff",
                "height": out_image.shape[1],
                "width": out_image.shape[2],
                "transform": out_transform
            })

            with rasterio.open(output_path, "w", **out_meta) as dest:
                dest.write(out_image)

        print(f"Saved clipped raster for  {state} and year {year}")

    except FileNotFoundError:
        print(f"File not found: {cdl_path}")
    except Exception as e:
        print(f"Error processing {year}: {e}")
        raise


# Main script
# Read state boundary file
state_bound = gpd.read_file(Path(state_bound_folder) / f"cb_2023_us_state_500k.shp")
# Reproject to the target CRS
state_bound = state_bound.to_crs(target_CRS)

for state in states:
    # Extract state boundaries
    select_state_bound = state_bound[state_bound['STUSPS'] == state]
    
    # Create 10 km outward buffer
    state_bound_outer = select_state_bound.buffer(state_buffer_m)

    # Create an output directory
    os.makedirs(cdl_output_folder, exist_ok=True)

    # Loop through each year to clip CDL raster
    for year in years_range:
        clip_cdl_using_state_boundaries(state_bound_outer, state, year, cdl_input_folder, cdl_output_folder)