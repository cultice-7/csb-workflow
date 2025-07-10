import geopandas as gp
import pandas as pd

gdf = gp.read_file("../data/Regrow/OSU_field_boudaries.geojson")

# Reproject to NAD83/Conus Albers
gdf_reproj = gdf.to_crs(epsg=5070)

# Calculate area in US survey acres (1 acre = 4046.8726 sq m)
gdf_reproj["area_acre"] = gdf_reproj.geometry.area/4046.8726

# Load Regrow data table
csv = pd.read_csv("../data/edited/Regrow/OH_main_crop_wide_coded.csv")

# Join shape file and table using boundary ID
gdf_joined = gdf_reproj.merge(csv, on='boundary_id', how='left')

# Save shape file
gdf_joined = gdf_joined.set_crs(epsg=5070)
gdf_joined.to_file("../data/edited/Regrow/regrow_clean.geojson")