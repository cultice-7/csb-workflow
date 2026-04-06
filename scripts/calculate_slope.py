from osgeo import gdal
from pathlib import Path

# Input and output paths
elevation_input = Path("data/Geo/elevation/elevation_reprojected.tif")
slope_output = Path("data/Geo/slope/slope_reprojected.tif")

# ensure parent directory exists
if not slope_output.parent.exists():
    slope_output.parent.mkdir(parents=True, exist_ok=True)

# Check the CRS of input reprojected elevation file
elevation_reprojected_tif = gdal.Open(elevation_input)
crs = elevation_reprojected_tif.GetProjection()
print("The CRS of the reprojected elevation file is", crs)

# Calculate slope in degrees
gdal.DEMProcessing(slope_output, elevation_input, "slope", format="GTiff", slopeFormat="degree")

print(f"Slope raster saved to {slope_output} using Albers Equal Area projection.")
