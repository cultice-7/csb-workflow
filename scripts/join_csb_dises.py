import geopandas as gpd
import pandas as pd
import numpy as np
import os
from shapely import intersection

# Input and output folders for Regrow
input_folder_CSB = "data/edited/CSB/"
output_folder_CSB = "data/edited/CSB/"

states = snakemake.params.states
years = snakemake.params.years

# Load DISES data
dises_shape = gpd.read_file("data/edited/DISES/dises_consolidated.gpkg")

# Rename DISES columns for clarity
dises_shape = dises_shape.add_suffix('_dises')

for year in years:
    for state in states:
            
        input_path_csb = os.path.join(input_folder_CSB, f"{state}_CSB{year}_shape_table.parquet")
        output_path_geospatial = os.path.join(output_folder_CSB, f"{state}_CSB{year}_dises_spatial.parquet")
        output_path_table = os.path.join(output_folder_CSB, f"{state}_CSB{year}_dises_table.parquet")
    
        # Load CSB data
        csb_clipped = gpd.read_parquet(input_path_csb)
        
        # Setting active geometry column
        csb_clipped = csb_clipped.set_geometry('geometry')
        dises_shape = dises_shape.set_geometry('geometry_dises')

        # Reproject all to an equal-area CRS (NAD83/CONUS Albers)
        csb_clipped = csb_clipped.set_crs(epsg=5070)
        dises_shape = dises_shape.to_crs(epsg=5070)

        # Preserve original geometry before buffering
        csb_clipped['original_geometry'] = csb_clipped.geometry
        dises_shape['original_geometry_dises'] = dises_shape.geometry

        # Create a buffered copy for spatial matching
        buffered = csb_clipped.copy()
        buffered['geometry'] = buffered.geometry.buffer(-10)
        buffered = buffered[buffered.is_valid & ~buffered.is_empty]

        # Perform intersection
        intersections = gpd.overlay(buffered, dises_shape, how='intersection')

        # Calculate overlap area
        intersections['overlap_area_dises'] = intersections.geometry.area / 4046.8564224

        # Keep only the largest overlap per CSB polygon
        largest_overlap = intersections.sort_values('overlap_area_dises', ascending=False).drop_duplicates('CSBID')

        # Merge attributes back to original CSB data
        columns_to_merge = largest_overlap.columns.difference(['geometry'])
        csb_dises = csb_clipped.merge(largest_overlap[columns_to_merge], on='CSBID', how='left', suffixes=('', '_temp'))

        # Drop temporary columns
        cols_to_drop = [col for col in csb_dises.columns if col.endswith('_temp')]
        csb_dises.drop(columns=cols_to_drop, inplace=True)

        # Add CSB-DISES assignment column
        csb_dises['field_assigned_dises'] = csb_dises['overlap_area_dises'].notna().astype(str)
        csb_dises['field_assigned_dises'] = csb_dises['field_assigned_dises'].replace({'True': 'Y', 'False': 'N'})
        
        # Calculate overlap area between Regrow and assigned DISES fields (in acres) and its share as % of CSB field area
        mask_overplap = csb_dises['field_assigned_dises'] == 'Y'
        csb_dises.loc[mask_overplap, 'overlap_area_dises'] = (
            intersection(csb_dises.loc[mask_overplap, 'original_geometry'], csb_dises.loc[mask_overplap, 'original_geometry_dises']).area) / 4046.8564224
        csb_dises.loc[mask_overplap, 'overlap_area_share_dises'] = (
            intersection(csb_dises.loc[mask_overplap, 'original_geometry'], csb_dises.loc[mask_overplap, 'original_geometry_dises']).area
            ) / csb_dises.loc[mask_overplap, 'original_geometry'].area
        
        # Restore original geometry
        csb_dises.drop(columns='original_geometry_dises', inplace=True)
        csb_dises = gpd.GeoDataFrame(csb_dises, geometry=csb_dises['original_geometry'], crs=csb_clipped.crs)
        csb_dises.drop(columns='original_geometry', inplace=True)

        # Add field match conditions (1,0, or NaN)
        csb_dises['crop_match_dises'] = np.where(
            csb_dises['field_crop_23_dises'].isna(),
            np.nan,
            (
                (csb_dises['field_crop_23_dises'] == csb_dises['CDL2023'])
            ).astype(int)
        )

        csb_dises['area_match_dises'] = np.where(
            csb_dises['field_size_dises'].isna(),
            np.nan,
            (
                (csb_dises['CSBACRES'] >= 0.75 * csb_dises['field_size_dises']) &
                (csb_dises['CSBACRES'] <= 1.25 * csb_dises['field_size_dises'])
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

        csb_dises['match_quality_dises'] = csb_dises.apply(determine_match_quality, axis=1)
        
        # Keep only the necessary columns
        cols_to_keep = [c for c in csb_dises.columns if c.endswith("_dises") or c in ["CSBID", "geometry"]]
        csb_dises = csb_dises[cols_to_keep]
        
        # Preview csb_dises
        print(csb_dises.head())

        # Save spatial joined CSB_dises shape file
        #csb_dises.to_parquet(output_path_geospatial, compression="zstd")

        # Save attribute table as CSV
        attribute_table = csb_dises.drop(columns='geometry')
        attribute_table.to_parquet(output_path_table, compression="zstd")
        print(f"CSB{year} and DISES for {state} are merged and saved")
