import geopandas as gp
import pandas as pd
import numpy as np

# Load data
regrow_shape = gp.read_file("data/edited/Regrow/regrow_clean.geojson")
dises_shape = gp.read_file("data/edited/DISES/dises_consolidated.gpkg")

# Reproject both to an equal-area CRS (NAD83/CONUS Albers)
regrow_shape = regrow_shape.to_crs(epsg=5070)
dises_shape = dises_shape.to_crs(epsg=5070)

# Preserve original geometry before buffering
regrow_shape['original_geometry'] = regrow_shape.geometry

# Create a buffered copy for spatial matching
buffered = regrow_shape.copy()
buffered['geometry'] = buffered.geometry.buffer(-50)
buffered = buffered[buffered.is_valid & ~buffered.is_empty]

# Perform intersection
intersections = gp.overlay(buffered, dises_shape, how='intersection')

# Calculate overlap area
intersections['overlap_area'] = intersections.geometry.area

# Keep only the largest overlap per Regrow polygon
largest_overlap = intersections.sort_values('overlap_area', ascending=False).drop_duplicates('boundary_id')

# Merge attributes back to original Regrow data
columns_to_merge = largest_overlap.columns.difference(['geometry'])
result = regrow_shape.merge(largest_overlap[columns_to_merge], on='boundary_id', how='left', suffixes=('', '_temp'))

# Drop temporary columns
cols_to_drop = [col for col in result.columns if col.endswith('_temp')]
result.drop(columns=cols_to_drop, inplace=True)

# Restore original geometry
result = gp.GeoDataFrame(result, geometry=result['original_geometry'], crs=regrow_shape.crs)
result.drop(columns='original_geometry', inplace=True)

# Add Regrow-DISES assignment column
result['dises_assigned'] = result['overlap_area'].notna().astype(str)
result['dises_assigned'] = result['dises_assigned'].replace('-1', '1')

# Add field match conditions (1, 0, or NaN)
result['crop_match'] = np.where(
    result['field_crop'].isna(),
    np.nan,
    (
        (result['field_crop'] == result['crop2023_1']) |
        (result['field_crop'] == result['crop2023_2']) |
        (result['field_crop'] == result['crop2023_3'])
    ).astype(int)
)

result['area_match'] = np.where(
    result['field_size'].isna(),
    np.nan,
    (
        (result['area_acre'] >= 0.8 * result['field_size']) &
        (result['area_acre'] <= 1.2 * result['field_size'])
    ).astype(int)
)

# Define match_quality based on crop_match and area_match
def determine_match_quality(row):
    if pd.isna(row['crop_match']) and pd.isna(row['area_match']):
        return np.nan
    elif row['crop_match'] == 1 and row['area_match'] == 1:
        return 'A'
    elif row['crop_match'] == 1:
        return 'B_crop'
    elif row['area_match'] == 1:
        return 'B_area'
    else:
        return 'F'

result['match_quality'] = result.apply(determine_match_quality, axis=1)

# Rename DISES columns for clarity
result.rename(columns={
    'field_crop': 'field_crop_dises',
    'field_name': 'field_name_dises',
    'field_size': 'field_size_dises',
    'from_data_table': 'from_data_table_dises'
}, inplace=True)

# Preview result
print(result.head())

# Save spatial joined Regrow shape
result.to_file("data/edited/Regrow/regrow_dises_spatialjoin.geojson", driver="GeoJSON")

# Save attribute table as CSV
attribute_table = result.drop(columns='geometry')
attribute_table.to_csv("data/edited/Regrow/regrow_dises_spatialjoin_table.csv", index=False)
