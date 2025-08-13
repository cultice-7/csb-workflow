from osgeo import gdal
import os

# Input and output paths
input_dem = "data/Geo/elevation/elevation.tif"
projected_dem = "data/Geo/elevation/elevation_projected.tif"
slope_output = "data/Geo/elevation/slope.tif"

# Reproject to an equal-area CRS (NAD83/CONUS Albers EPSG:5070)
gdal.Warp(projected_dem, input_dem, dstSRS="EPSG:5070", format="GTiff")

# Calculate slope in degrees
gdal.DEMProcessing(slope_output, projected_dem, "slope", format="GTiff", slopeFormat="degree")

if os.path.exists(projected_dem):
    os.remove(projected_dem)
    print(f"Deleted intermediate file: {projected_dem}")

print(f"Slope raster saved to {slope_output} using Albers Equal Area projection.")
