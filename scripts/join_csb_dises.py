#untested

import geopandas as gp
import pandas as pd
import numpy as np

# Load data
csb1623_clipped = gp.read_file("data/edited/CSB/CSB1623_clipped.gdb")
csb1724_clipped = gp.read_file("data/edited/CSB/CSB1623_clipped.gdb")
dises_shape = gp.read_file("data/edited/DISES/dises_consolidated.gpkg")

# Reproject all to an equal-area CRS (NAD83/CONUS Albers)
csb1623_clipped = csb1623_clipped.set_crs(epsg=5070)
csb1724_clipped = csb1724_clipped.set_crs(epsg=5070)
dises_shape = dises_shape.to_crs(epsg=5070)

# Preserve original geometry before buffering
csb1623_clipped['original_geometry'] = csb1623_clipped.geometry
csb1724_clipped['original_geometry'] = csb1724_clipped.geometry

# Create a buffered copy for spatial matching
buffered_1623 = csb1623_clipped.copy()
buffered_1623['geometry'] = buffered_1623.geometry.buffer(-30)
buffered_1623 = buffered_1623[buffered_1623.is_valid & ~buffered_1623.is_empty]

buffered_1724 = csb1724_clipped.copy()
buffered_1724['geometry'] = buffered_1724.geometry.buffer(-30)
buffered_1724 = buffered_1724[buffered_1724.is_valid & ~buffered_1724.is_empty]

# Perform intersection
intersections_1623 = gp.overlay(buffered_1623, dises_shape, how='intersection')
intersections_1724 = gp.overlay(buffered_1724, dises_shape, how='intersection')

# Calculate overlap area
intersections_1623['overlap_area'] = intersections_1623.geometry.area
intersections_1724['overlap_area'] = intersections_1724.geometry.area

# Keep only the largest overlap per CSB polygon
largest_overlap_1623 = intersections_1623.sort_values('overlap_area', ascending=False).drop_duplicates('CSBID')
largest_overlap_1724 = intersections_1724.sort_values('overlap_area', ascending=False).drop_duplicates('CSBID')

# Merge attributes back to original CSB data
columns_to_merge_1623 = largest_overlap_1623.columns.difference(['geometry'])
result_1623 = csb1623_clipped.merge(largest_overlap_1623[columns_to_merge_1623], on='CSBID', how='left', suffixes=('', '_temp'))

columns_to_merge_1724 = largest_overlap_1724.columns.difference(['geometry'])
result_1724 = csb1724_clipped.merge(largest_overlap_1724[columns_to_merge_1724], on='CSBID', how='left', suffixes=('', '_temp'))

# Drop temporary columns
cols_to_drop_1623 = [col for col in result_1623.columns if col.endswith('_temp')]
result_1623.drop(columns=cols_to_drop_1623, inplace=True)

cols_to_drop_1724 = [col for col in result_1724.columns if col.endswith('_temp')]
result_1724.drop(columns=cols_to_drop_1724, inplace=True)

# Restore original geometry
result_1623 = gp.GeoDataFrame(result_1623, geometry=result_1623['original_geometry'], crs=csb1623_clipped.crs)
result_1623.drop(columns='original_geometry', inplace=True)

result_1724 = gp.GeoDataFrame(result_1724, geometry=result_1724['original_geometry'], crs=csb1724_clipped.crs)
result_1724.drop(columns='original_geometry', inplace=True)

# Add CSB-DISES assignment column
result_1623['dises_assigned'] = result_1623['overlap_area'].notna()
result_1724['dises_assigned'] = result_1724['overlap_area'].notna()

# Add field match conditions (1,0, or NaN)
result_1623['crop_match'] = np.where(
    result_1623['field_crop'].isna(),
    np.nan,
    (
        (result_1623['field_crop'] == result_1623['crop2023_1']) |
        (result_1623['field_crop'] == result_1623['crop2023_2']) |
        (result_1623['field_crop'] == result_1623['crop2023_3'])
    ).astype(int)
)

result_1623['area_match'] = np.where(
    result_1623['field_size'].isna(),
    np.nan,
    (
        (result_1623['area_acre'] >= 0.8 * result_1623['field_size']) &
        (result_1623['area_acre'] <= 1.2 * result_1623['field_size'])
    ).astype(int)
)



result_1724['crop_match'] = np.where(
    result_1724['field_crop'].isna(),
    np.nan,
    (
        (result_1724['field_crop'] == result_1724['crop2023_1']) |
        (result_1724['field_crop'] == result_1724['crop2023_2']) |
        (result_1724['field_crop'] == result_1724['crop2023_3'])
    ).astype(int)
)

result_1724['area_match'] = np.where(
    result_1724['field_size'].isna(),
    np.nan,
    (
        (result_1724['area_acre'] >= 0.8 * result_1724['field_size']) &
        (result_1724['area_acre'] <= 1.2 * result_1724['field_size'])
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

result_1623['match_quality'] = result_1623.apply(determine_match_quality, axis=1)


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

result_1724['match_quality'] = result_1724.apply(determine_match_quality, axis=1)


# Rename DISES columns for clarity
result_1623.rename(columns={
    'field_crop': 'field_crop_dises',
    'field_name': 'field_name_dises',
    'field_size': 'field_size_dises',
    'from_data_table': 'from_data_table_dises'
}, inplace=True)

result_1724.rename(columns={
    'field_crop': 'field_crop_dises',
    'field_name': 'field_name_dises',
    'field_size': 'field_size_dises',
    'from_data_table': 'from_data_table_dises'
}, inplace=True)

# Preview result
print(result_1623.head())
print(result_1724.head())


# Save spatial joined CSB shape
result_1623.to_file("data/edited/CSB/CSB1623_dises_spatialjoin.geojson", driver="GeoJSON")
result_1724.to_file("data/edited/CSB/CSB1724_dises_spatialjoin.geojson", driver="GeoJSON")

# Save attribute table as CSV
attribute_table_1623 = result_1623.drop(columns='geometry')
attribute_table_1724 = result_1724.drop(columns='geometry')

attribute_table_1623.to_csv("data/edited/CSB/CSB1623_dises_spatialjoin_table.csv", index=False)
attribute_table_1724.to_csv("data/edited/CSB/CSB1724_dises_spatialjoin_table.csv", index=False)