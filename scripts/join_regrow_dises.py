# IN DRAFT

import geopandas as gp
import pandas as pd

regrow_shape = gp.read_file("data/edited/Regrow/regrow_clean.geojson")
dises_shape = gp.read_file("data/edited/DISES/dises_consolidated.gpkg")

# Reproject both to an equal-area CRS (NAD83/CONUS Albers)
regrow_shape = regrow_shape.to_crs(epsg=5070)
dises_shape = dises_shape.to_crs(epsg=5070)

# Perform spatial join equivalent to Arcgis's largest overlap with (-)50 meters buffer.
regrow_shape['geometry'] = regrow_shape.geometry.buffer(-50)
regrow_shape = regrow_shape[regrow_shape.is_valid & ~regrow_shape.is_empty]
intersections = gp.overlay(regrow_shape, dises_shape, how='intersection')
intersections['overlap_area'] = intersections.geometry.area
largest_overlap = intersections.sort_values('overlap_area', ascending=False).drop_duplicates('regrow_id')
columns_to_merge = largest_overlap.columns.difference(['geometry'])
result = regrow_shape.merge(largest_overlap[columns_to_merge], on='boundary_id', how='left')

# Restore original geometry
result = gp.GeoDataFrame(result, geometry=regrow_shape, crs=regrow_shape.crs)

# Add match condition


# Save spatial joined Regrow shape
result.to_file("data/edited/Regrow/regrow_dises_spatialjoin.geojson", driver="GeoJSON")

# Preview result
print(result.head())
