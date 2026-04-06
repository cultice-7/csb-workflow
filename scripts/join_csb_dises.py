import geopandas as gpd
import pandas as pd
import numpy as np
import os
from shapely import intersection
from pathlib import Path

# Import parameters from Snakemake
CSB_input_folder = snakemake.params.CSB_input_dir
DISES_input_folder = snakemake.params.DISES_input_dir
CSB_DISES_output_folder = snakemake.params.CSB_DISES_output_dir
states = snakemake.params.states
years = snakemake.params.CSB_years
buffer_margin = snakemake.params.buffer_margin
area_match_lower_bound = snakemake.params["area_match_coefs"][0]
area_match_upper_bound = snakemake.params["area_match_coefs"][1]
target_CRS = snakemake.params.target_CRS


# Load DISES data
dises_shape_table = gpd.read_parquet(Path(DISES_input_folder) / "DISES_shape_table.parquet")

# Rename DISES columns for clarity
dises_shape_table = dises_shape_table.add_suffix('_dises')

for year in years:
    for state in states:
            
        csb_geometry_input_path = os.path.join(CSB_input_folder, f"{state}_CSB{year}_CSBID_geometry.parquet")
        csb_table_input_path = os.path.join(CSB_input_folder, f"{state}_CSB{year}_table.parquet")
        output_path_geospatial = os.path.join(CSB_DISES_output_folder, f"{state}_CSB{year}_dises_spatial.parquet")
        output_path_table = os.path.join(CSB_DISES_output_folder, f"{state}_CSB{year}_dises_table.parquet")
    
        # Load CSB data
        csb_shape = gpd.read_parquet(csb_geometry_input_path)
        csb_table = pd.read_parquet(csb_table_input_path)
        
        # Create CSB shape_table file by merging csb_shape and csb_table
        csb_shape_table = csb_shape.merge(csb_table, on="CSBID", how="left")
        
        # Setting active geometry column
        csb_shape_table = csb_shape_table.set_geometry('geometry')
        dises_shape_table = dises_shape_table.set_geometry('geometry_dises')

        # Reproject all to an equal-area CRS (NAD83/CONUS Albers)
        csb_shape_table = csb_shape_table.set_crs(target_CRS)
        dises_shape_table = dises_shape_table.to_crs(target_CRS)

        # Preserve original geometry before buffering
        csb_shape_table['original_geometry'] = csb_shape_table.geometry
        dises_shape_table['original_geometry_dises'] = dises_shape_table.geometry

        # Create a buffered copy for spatial matching
        buffered = csb_shape_table.copy()
        buffered['geometry'] = buffered.geometry.buffer(buffer_margin)
        buffered = buffered[buffered.is_valid & ~buffered.is_empty]

        # Perform intersection
        intersections = gpd.overlay(buffered, dises_shape_table, how='intersection')

        # Calculate overlap area
        intersections['overlap_area_dises'] = intersections.geometry.area / 4046.8564224

        # Keep only the largest overlap per CSB polygon
        largest_overlap = intersections.sort_values('overlap_area_dises', ascending=False).drop_duplicates('CSBID')

        # Merge attributes back to original CSB data
        columns_to_merge = largest_overlap.columns.difference(['geometry'])
        csb_dises = csb_shape_table.merge(largest_overlap[columns_to_merge], on='CSBID', how='left', suffixes=('', '_temp'))

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
        csb_dises = gpd.GeoDataFrame(csb_dises, geometry=csb_dises['original_geometry'], crs=csb_shape_table.crs)
        csb_dises.drop(columns='original_geometry', inplace=True)
        
        # Crop mapping to match crop categories between CSB and DISES
        # Define mapping to match "field_crop_23_dises" with CSB crop categories
        crop_mapping = {
        1: 1,
        2: 5,
        3: np.nan}
        
        # Create a temporary column in csb_dises
        csb_dises["field_crop_23_dises_mapped"] = csb_dises["field_crop_23_dises"].replace(crop_mapping)

        # Add field match conditions (1,0, or NaN)
        csb_dises['crop_match_dises'] = np.where(
            (csb_dises['field_crop_23_dises'].isna()) | (csb_dises['CDL2023'].isna()),
            np.nan,
            (
                (csb_dises['field_crop_23_dises'] == csb_dises['CDL2023'])
            ).astype(int)
        )

        csb_dises['area_match_dises'] = np.where(
            (csb_dises['field_size_dises'].isna()) | (csb_dises['CSBACRES'].isna()),
            np.nan,
            (
                (csb_dises['CSBACRES'] >= area_match_lower_bound * csb_dises['field_size_dises']) &
                (csb_dises['CSBACRES'] <= area_match_upper_bound * csb_dises['field_size_dises'])
            ).astype(int)
        )
        
        # Drop temporary column in regrow_dises
        csb_dises.drop(columns="field_crop_23_dises_mapped", inplace=True)

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
