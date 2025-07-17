#untested

import os
import requests

# Define the y and x values
y_range = range(39, 50)
x_range = range(81, 92)

# Base URL pattern
base_url = "https://prd-tnm.s3.amazonaws.com/StagedProducts/Elevation/1/TIFF/current/n{yy}w{xxx}/USGS_1_n{y}w{xxx}.tif"

# Folder to save elevation files
download_folder = "data/Geo/elevation"
os.makedirs(download_folder, exist_ok=True)

# Loop through all combinations
for y in y_range:
    for x in x_range:
        y_str = str(y).zfill(2)
        x_str = str(x).zfill(3)
        url = base_url.format(yy=y_str, xxx=x_str)
        filename = f"USGS_1_n{y_str}w{x_str}.tif"
        filepath = os.path.join(download_folder, filename)

        response = requests.get(url, stream=True)
        if response.status_code == 200:
            with open(filepath, "wb") as f:
                for chunk in response.iter_content(chunk_size=1048576):
                    f.write(chunk)

print("All files downloaded.")


import glob
import rasterio
from rasterio.merge import merge

# Find all .tif files in the folder
tif_files = glob.glob(os.path.join(download_folder, "*.tif"))

# Open all the raster files
src_files_to_mosaic = [rasterio.open(fp) for fp in tif_files]

# Merge the rasters
mosaic, out_trans = merge(src_files_to_mosaic)

# Copy metadata from one of the source files
out_meta = src_files_to_mosaic[0].meta.copy()
out_meta.update({
    "driver": "GTiff",
    "height": mosaic.shape[1],
    "width": mosaic.shape[2],
    "transform": out_trans,
    "count": mosaic.shape[0]
})

# Save the merged raster
output_path = os.path.join(download_folder, "elevation.tif")
with rasterio.open(output_path, "w", **out_meta) as dest:
    dest.write(mosaic)

print("Merged raster saved as elevation.tif")
