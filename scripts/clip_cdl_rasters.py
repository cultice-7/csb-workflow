import geopandas as gp
import os
import rasterio
from rasterio.mask import mask

# Filter state boundary file to 3 states (OH, MI, IN)
state_bound = gp.read_file("data/Census/state_bound/cb_2018_us_state_500k.shp")
select_states = state_bound[state_bound['STATEFP'].isin(['39', '18', '26'])]

# Reproject CRS to NAD83 to match CDL rasters
select_states = select_states.to_crs(epsg=5070)

# Define input and output folders
os.makedirs("data/edited/CDL", exist_ok=True)

cdl_folder = "data/CDL/"
clipped_cdl_folder = "data/edited/CDL/"

# Loop through each year to clip CDL raster
for year in range(2014, 2025):
    cdl_path = os.path.join(cdl_folder, f"{year}_30m_cdls", f"{year}_30m_cdls.tif")
    output_path = os.path.join(clipped_cdl_folder, f"{year}_30m_cdls_clipped.tif")

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