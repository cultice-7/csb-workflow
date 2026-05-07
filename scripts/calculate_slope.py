from osgeo import gdal
from pathlib import Path

# Import parameters from Snakemake
elevation_input_folder = snakemake.params.elevation_input_dir
slope_output_folder = snakemake.params.slope_output_dir


# Input and output paths
elevation_input = Path(elevation_input_folder) / "elevation_reprojected.tif"
slope_output = Path(slope_output_folder) / "slope_reprojected.tif"

# Ensure parent directory exists
if not slope_output.parent.exists():
    slope_output.parent.mkdir(parents=True, exist_ok=True)

# Check the CRS of input reprojected elevation file
elevation_reprojected_tif = gdal.Open(elevation_input)
crs = elevation_reprojected_tif.GetProjection()
print("The CRS of the reprojected elevation file is", crs)

# Calculate slope in degrees
gdal.DEMProcessing(slope_output, elevation_input, "slope", format="GTiff", slopeFormat="degree")

print(f"Slope raster saved to {slope_output}")
