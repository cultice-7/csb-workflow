import geopandas as gpd
import pandas as pd
import numpy as np
import os
from shapely import intersection

# Input and output folders for Regrow
input_folder_Regrow = "data/edited/Regrow/"
output_folder_Regrow = "data/edited/Regrow/"

# Load DISES data
dises_shape = gpd.read_file("data/edited/DISES/dises_consolidated.gpkg")

# Rename all DISES columns for clarity
dises_shape = dises_shape.add_suffix('_dises')

states = snakemake.params.states

for state in states:
    
    input_path_regrow = os.path.join(input_folder_Regrow, f"{state}_regrow_shape_table.parquet")
    output_path_geospatial = os.path.join(output_folder_Regrow, f"{state}_regrow_dises_spatial.parquet")
    output_path_table = os.path.join(output_folder_Regrow, f"{state}_regrow_dises_table.parquet")
    
    # Load Regrow data
    regrow_shape = gpd.read_parquet(input_path_regrow)

    # Setting active geometry column
    regrow_shape = regrow_shape.set_geometry('geometry')
    dises_shape = dises_shape.set_geometry('geometry_dises')

    # Reproject both to an equal-area CRS (NAD83/CONUS Albers)
    regrow_shape = regrow_shape.to_crs(epsg=5070)
    dises_shape = dises_shape.to_crs(epsg=5070)

    # Preserve original geometry before buffering
    regrow_shape['original_geometry'] = regrow_shape.geometry
    dises_shape['original_geometry_dises'] = dises_shape.geometry

    # Create a buffered copy for spatial matching
    buffered = regrow_shape.copy()
    buffered['geometry'] = buffered.geometry.buffer(-10)
    buffered = buffered[buffered.is_valid & ~buffered.is_empty]

    # Perform intersection
    intersections = gpd.overlay(buffered, dises_shape, how='intersection')

    # Calculate overlap area (in acres)
    intersections['overlap_area_dises'] = intersections.geometry.area / 4046.8564224
    
    # Keep only the largest overlap per Regrow polygon
    largest_overlap = intersections.sort_values('overlap_area_dises', ascending=False).drop_duplicates('field_id')

    # Merge attributes back to original Regrow data
    columns_to_merge = largest_overlap.columns.difference(['geometry'])
    regrow_dises = regrow_shape.merge(largest_overlap[columns_to_merge], on='field_id', how='left', suffixes=('', '_temp'))

    # Drop temporary columns
    cols_to_drop = [col for col in regrow_dises.columns if col.endswith('_temp')]
    regrow_dises.drop(columns=cols_to_drop, inplace=True)

    # Add Regrow-DISES assignment column
    regrow_dises['field_assigned_dises'] = regrow_dises['overlap_area_dises'].notna().astype(str)
    regrow_dises['field_assigned_dises'] = regrow_dises['field_assigned_dises'].replace({'True': 'Y', 'False': 'N'})
    
    # Calculate overlap area between Regrow and assigned DISES fields (in acres) and its share as % of Regrow field area
    mask_overplap = regrow_dises['field_assigned_dises'] == 'Y'
    regrow_dises.loc[mask_overplap, 'overlap_area_dises'] = (
        intersection(regrow_dises.loc[mask_overplap, 'original_geometry'], regrow_dises.loc[mask_overplap, 'original_geometry_dises']).area) / 4046.8564224
    regrow_dises.loc[mask_overplap, 'overlap_area_share_dises'] = (
        intersection(regrow_dises.loc[mask_overplap, 'original_geometry'], regrow_dises.loc[mask_overplap, 'original_geometry_dises']).area
        ) / regrow_dises.loc[mask_overplap, 'original_geometry'].area
    
    # Restore original geometry
    regrow_dises.drop(columns='original_geometry_dises', inplace=True)
    regrow_dises = gpd.GeoDataFrame(regrow_dises, geometry=regrow_dises['original_geometry'], crs=regrow_shape.crs)
    regrow_dises.drop(columns='original_geometry', inplace=True)

    # Add field match conditions (1, 0, or NaN)
    regrow_dises['crop_match_dises'] = np.where(
        regrow_dises['field_crop_23_dises'].isna(),
        np.nan,
        (
            (regrow_dises['field_crop_23_dises'] == regrow_dises['crop_23_1']) |
            (regrow_dises['field_crop_23_dises'] == regrow_dises['crop_23_2'])
        ).astype(int)
    )

    regrow_dises['area_match_dises'] = np.where(
        regrow_dises['field_size_dises'].isna(),
        np.nan,
        (
            (regrow_dises['area_acre'] >= 0.75 * regrow_dises['field_size_dises']) &
            (regrow_dises['area_acre'] <= 1.25 * regrow_dises['field_size_dises'])
        ).astype(int)
    )

    # Define match_quality based on crop_match_dises and area_match_dises
    def determine_match_quality(row):
        if pd.isna(row['crop_match_dises']) and pd.isna(row['area_match_dises']):
            return np.nan
        elif row['crop_match_dises'] == 1 and row['area_match_dises'] == 1:
            return 'A'
        elif row['crop_match_dises'] == 1:
            return 'B_crop'
        elif row['area_match_dises'] == 1:
            return 'B_area'
        else:
            return 'F'

    regrow_dises['match_quality_dises'] = regrow_dises.apply(determine_match_quality, axis=1)
    
    # Keep only the necessary columns
    cols_to_keep = [c for c in regrow_dises.columns if c.endswith("_dises") or c in ["field_id", "geometry"]]
    regrow_dises = regrow_dises[cols_to_keep]
    
    # Preview result
    print(regrow_dises.head())

    # Save spatial joined Regrow_dises shape file
    #regrow_dises.to_parquet(output_path_geospatial, compression="zstd")

    # Save attribute table as CSV
    attribute_table = regrow_dises.drop(columns='geometry')
    attribute_table.to_parquet(output_path_table, compression="zstd")
    print(f"Regrow and DISES for {state} are merged and saved")
