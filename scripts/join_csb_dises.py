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
overlap_threshold = snakemake.params.overlap_threshold
target_CRS = snakemake.params.target_CRS


# Load DISES data
dises_shape_table = gpd.read_parquet(Path(DISES_input_folder) / "DISES_shape_table.parquet")

# Rename DISES columns for clarity
dises_shape_table = dises_shape_table.add_suffix('_dises')

for year in years:
    csb_dises_concat = pd.DataFrame()
    
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
        csb_shape_table = csb_shape_table.to_crs(target_CRS)
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
        csb_dises['parcel_assigned_dises'] = csb_dises['overlap_area_dises'].notna().astype(str)
        csb_dises['parcel_assigned_dises'] = csb_dises['parcel_assigned_dises'].replace({'True': 'Y', 'False': 'N'})
        
        # Calculate overlap area between Regrow and assigned DISES fields (in acres) and its share as % of CSB field area
        mask_overplap = csb_dises['parcel_assigned_dises'] == 'Y'
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
            (csb_dises['field_crop_23_dises_mapped'].isna()) | (csb_dises['CDL2023'].isna()),
            np.nan,
            (
                (csb_dises['field_crop_23_dises_mapped'] == csb_dises['CDL2023'])
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
        
        # Collect all merged regrow_dises together
        if csb_dises_concat.empty:
            csb_dises_concat = csb_dises.copy()
        else:
            csb_dises_concat = pd.concat([csb_dises_concat, csb_dises], axis=0, ignore_index=True)
    
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

    # Add representative field attribute to CSB fields
    def assign_representative_field(csb_dises, overlap_threshold):

        mask_field_assigned = csb_dises['parcel_assigned_dises'] == 'Y'
        mask_survey_responded = csb_dises['survey_responded_dises'] == 'Y'
        mask_overlap_area = csb_dises["overlap_area_share_dises"] >= overlap_threshold
        
        csb_dises_matching = csb_dises[(mask_field_assigned) & (mask_survey_responded) & (mask_overlap_area)].reset_index(drop=True)
        
        level_1 = (
            # A: size + crop match
            csb_dises_matching["match_quality_dises"] == "A"
        )

        level_2 = (
            # B_area: size match only, crop missing
            (
                (csb_dises_matching["match_quality_dises"] == "B_area") &
                (csb_dises_matching["field_crop_23_dises"].isna())
            ) |
            # B_crop: crop match only, size missing
            (
                (csb_dises_matching["match_quality_dises"] == "B_crop") &
                (csb_dises_matching["field_size_dises"].isna())
            )
        )

        level_3 = (
            # Size match only, crop mismatch
            (
                (csb_dises_matching["match_quality_dises"] == "B_area") &
                (csb_dises_matching["field_crop_23_dises"].notna())
            )
        )

        level_4 = (
            # RF size missing, crop mismatch
            (
                (csb_dises_matching["match_quality_dises"] == "F") &
                (csb_dises_matching["field_size_dises"].isna()) &
                (csb_dises_matching["field_crop_23_dises"].notna())
            ) |
            # Both RF size and crop missing, only single parcel
            (
                (csb_dises_matching["field_crop_23_dises"].isna()) &
                (csb_dises_matching["field_size_dises"].isna()) &
                (csb_dises_matching["n_parcels_dises"] == 1)
            )
        )

        # Collapse the 4 levels into one column. np.select takes the first matching
        # condition in the list, so Level 1 wins whenever a row happens to qualify for more than one level.
        csb_dises_matching["RF_assignment_dises"] = np.select(
            [level_1, level_2, level_3, level_4],
            ["Level 1", "Level 2", "Level 3", "Level 4"],
            default=None
        )
        RF_rank_dises = csb_dises_matching["RF_assignment_dises"].map(
            {"Level 1": 1, "Level 2": 2, "Level 3": 3, "Level 4": 4}
        )

        csb_dises_matching["area_diff"] = (csb_dises_matching["CSBACRES"] - csb_dises_matching["field_size_dises"]).abs()
        csb_dises_matching["RF_rank_dises"] = RF_rank_dises
        csb_dises_matching = csb_dises_matching.sort_values(
            by=["comp_id_dises", "RF_rank_dises", "area_diff", "overlap_area_share_dises"],
            ascending=[True, True, True, False]
        )

        csb_dises_matching = csb_dises_matching.drop_duplicates(subset=["comp_id_dises"], keep="first")
        csb_dises_matching.reset_index(drop=True, inplace=True)

        col_to_keep = csb_dises_matching.filter(regex="CSBID|RF_assignment_dises").columns
        return(csb_dises_matching[col_to_keep])


    representative_field_indicator = assign_representative_field(csb_dises_concat, overlap_threshold)
    del csb_dises_concat


    for state in states:
        output_path_table = os.path.join(CSB_DISES_output_folder, f"{state}_CSB{year}_dises_table.parquet")
        csb_dises = pd.read_parquet(output_path_table)

        csb_dises = csb_dises.merge(representative_field_indicator, on='CSBID', how='left')

        csb_dises.to_parquet(output_path_table, compression="zstd")
        print(f"CSB{year} and DISES with RF indicators for {state} are saved")
