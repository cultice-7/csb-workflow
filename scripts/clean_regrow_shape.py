import geopandas as gp
import pandas as pd

gdf = gp.read_file("data/Regrow/OH_field_boundaries.geojson")

# Reproject to NAD83/Conus Albers
gdf_reproj = gdf.to_crs(epsg=5070)

# Rename the column "boundary_id"
gdf_reproj.rename(columns={'boundary_id': 'field_id'}, inplace = True)

# Calculate area in US survey acres (1 acre = 4046.8564224 sq m)
gdf_reproj["area_acre"] = gdf_reproj.geometry.area/4046.8564224

# Load Regrow data table
csv_main_crop = pd.read_csv("data/edited/Regrow/OH_main_crop_wide_coded.csv")
csv_tillage = pd.read_csv("data/edited/Regrow/OH_tillage_wide_coded.csv")
csv_cover_crop = pd.read_csv("data/edited/Regrow/OH_cover_crop_wide_coded.csv")

# Join shape file and table using boundary ID
gdf_joined = gdf_reproj.merge(csv_main_crop, on = 'field_id', how = 'left')
gdf_joined = gdf_joined.merge(csv_tillage, on = 'field_id', how = 'left')
gdf_joined = gdf_joined.merge(csv_cover_crop, on = 'field_id', how = 'left')

# Replace NaN with a placeholder string 
gdf_joined = gdf_joined.fillna('NA')

print(gdf_joined.shape)

# Save shape file
gdf_joined = gdf_joined.set_crs(epsg=5070)
gdf_joined.to_file("data/edited/Regrow/OH_regrow_clean.geojson", driver = "GeoJSON")