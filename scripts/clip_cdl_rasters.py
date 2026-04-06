import geopandas as gpd
import os
import rasterio
from rasterio.mask import mask
from rasterio.windows import from_bounds
from pathlib import Path

# Import parameters from the Snakemake
state_bound_folder = snakemake.params.state_bound_dir
states_codes = snakemake.params.states_code
years_range = snakemake.params.years
target_CRS = snakemake.params.target_CRS
cdl_input_folder = snakemake.params.cdl_input_dir
cdl_output_folder = snakemake.params.cdl_output_dir


def clip_cdl_using_state_boundaries(select_states, year, cdl_input_folder, cdl_output_folder):
    
    cdl_path = os.path.join(cdl_input_folder, f"{year}_30m_cdls", f"{year}_30m_cdls.tif")
    output_path = os.path.join(cdl_output_folder, f"{year}_30m_cdls_clipped.tif")

    print(f"Processing {year}...")

    try:
        with rasterio.open(cdl_path) as src:
            out_image, out_transform = mask(src, select_states.geometry, crop=True)
            out_meta = src.meta.copy()
            out_meta.update({
                "driver": "GTiff",
                "height": out_image.shape[1],
                "width": out_image.shape[2],
                "transform": out_transform
            })

            with rasterio.open(output_path, "w", **out_meta) as dest:
                dest.write(out_image)

        print(f"Saved clipped raster for {year}")

    except FileNotFoundError:
        print(f"File not found: {cdl_path}")
    except Exception as e:
        print(f"Error processing {year}: {e}")


# Main script
state_bound = gpd.read_file(Path(state_bound_folder) / f"cb_2023_us_state_500k.shp")
select_states = state_bound[state_bound['STATEFP'].isin(states_codes)]

# Reproject CRS to match CDL rasters
select_states = select_states.to_crs(target_CRS)

# Create an output directory
os.makedirs(cdl_output_folder, exist_ok=True)

# Loop through each year to clip CDL raster
for year in years_range:
    clip_cdl_using_state_boundaries(select_states, year, cdl_input_folder, cdl_output_folder)