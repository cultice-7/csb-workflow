import geopandas as gpd
import pandas as pd
import numpy as np
import os
from shapely import intersection
from pathlib import Path

# Import parameters from Snakemake
regrow_input_folder = snakemake.params.regrow_input_dir
DISES_input_folder = snakemake.params.DISES_input_dir
regrow_DISES_output_folder = snakemake.params.regrow_DISES_output_dir
states = snakemake.params.states
buffer_margin = snakemake.params.buffer_margin
area_match_lower_bound = snakemake.params["area_match_coefs"][0]
area_match_upper_bound = snakemake.params["area_match_coefs"][1]
target_CRS = snakemake.params.target_CRS


# Load DISES data
dises_shape_table = gpd.read_parquet(Path(DISES_input_folder) / "DISES_shape_table.parquet")

# Rename all DISES columns for clarity
dises_shape_table = dises_shape_table.add_suffix('_dises')

regrow_dises_concat = pd.DataFrame()

for state in states:
    
    regrow_geometry_input_path = os.path.join(regrow_input_folder, f"{state}_regrow_fieldID_geometry.parquet")
    regrow_table_input_path = os.path.join(regrow_input_folder, f"{state}_regrow_table.parquet")
    output_path_geospatial = os.path.join(regrow_DISES_output_folder, f"{state}_regrow_dises_spatial.parquet")
    output_path_table = os.path.join(regrow_DISES_output_folder, f"{state}_regrow_dises_table.parquet")
    
    # Load Regrow data
    regrow_shape = gpd.read_parquet(regrow_geometry_input_path)
    regrow_table = pd.read_parquet(regrow_table_input_path)
    
    # Create regrow shape_table file by merging regrow_shape and regrow_table
    regrow_shape_table = regrow_shape.merge(regrow_table, on="field_id", how="left")

    # Setting active geometry column
    regrow_shape_table = regrow_shape_table.set_geometry('geometry')
    dises_shape_table = dises_shape_table.set_geometry('geometry_dises')

    # Reproject both to an equal-area CRS (NAD83/CONUS Albers)
    regrow_shape_table = regrow_shape_table.to_crs(target_CRS)
    dises_shape_table = dises_shape_table.to_crs(target_CRS)

    # Preserve original geometry before buffering
    regrow_shape_table['original_geometry_regrow'] = regrow_shape_table.geometry
    dises_shape_table['original_geometry_dises'] = dises_shape_table.geometry

    # Create a buffered copy for spatial matching
    buffered = regrow_shape_table.copy()
    buffered['geometry'] = buffered.geometry.buffer(buffer_margin)
    buffered = buffered[buffered.is_valid & ~buffered.is_empty]

    # Perform intersection
    intersections = gpd.overlay(buffered, dises_shape_table, how='intersection')

    # Calculate overlap area (in acres)
    intersections['overlap_area_dises'] = intersections.geometry.area / 4046.8564224
    
    # Keep only the largest overlap per Regrow polygon
    largest_overlap = intersections.sort_values('overlap_area_dises', ascending=False).drop_duplicates('field_id')

    # Merge attributes back to original Regrow data
    columns_to_merge = largest_overlap.columns.difference(['geometry'])
    regrow_dises = regrow_shape_table.merge(largest_overlap[columns_to_merge], on='field_id', how='left', suffixes=('', '_temp'))

    # Drop temporary columns
    cols_to_drop = [col for col in regrow_dises.columns if col.endswith('_temp')]
    regrow_dises.drop(columns=cols_to_drop, inplace=True)

    # Add Regrow-DISES assignment column
    regrow_dises['field_assigned_dises'] = regrow_dises['overlap_area_dises'].notna().astype(str)
    regrow_dises['field_assigned_dises'] = regrow_dises['field_assigned_dises'].replace({'True': 'Y', 'False': 'N'})
    
    # Calculate overlap area between Regrow and assigned DISES fields (in acres) and its share as % of Regrow field area
    mask_overplap = regrow_dises['field_assigned_dises'] == 'Y'
    regrow_dises.loc[mask_overplap, 'overlap_area_dises'] = (
        intersection(regrow_dises.loc[mask_overplap, 'original_geometry_regrow'], regrow_dises.loc[mask_overplap, 'original_geometry_dises']).area) / 4046.8564224
    regrow_dises.loc[mask_overplap, 'overlap_area_share_dises'] = (
        intersection(regrow_dises.loc[mask_overplap, 'original_geometry_regrow'], regrow_dises.loc[mask_overplap, 'original_geometry_dises']).area
        ) / regrow_dises.loc[mask_overplap, 'original_geometry_regrow'].area
    
    # Restore original geometry
    regrow_dises.drop(columns='original_geometry_dises', inplace=True)
    regrow_dises = gpd.GeoDataFrame(regrow_dises, geometry=regrow_dises['original_geometry_regrow'], crs=regrow_shape_table.crs)
    regrow_dises.drop(columns='original_geometry_regrow', inplace=True)

    # Crop mapping to match crop categories between Regrow and DISES
    # Define mapping to match "field_crop_23_dises" with Regrow crop categories
    crop_mapping = {
    1: 1,
    2: 5,
    3: np.nan}
    
    # Create a temporary column in regrow_dises
    regrow_dises["field_crop_23_dises_mapped"] = regrow_dises["field_crop_23_dises"].replace(crop_mapping)
    
    # Add field match conditions (1, 0, or NaN)
    regrow_dises['crop_match_dises'] = np.where(
        (regrow_dises['field_crop_23_dises_mapped'].isna()) | (regrow_dises['crop_23_1'].isna()),
        np.nan,
        (
            (regrow_dises['field_crop_23_dises_mapped'] == regrow_dises['crop_23_1'].astype(float)) |
            (regrow_dises['field_crop_23_dises_mapped'] == regrow_dises['crop_23_2'].astype(float))
        ).astype(int)
    )

    regrow_dises['area_match_dises'] = np.where(
        (regrow_dises['field_size_dises'].isna()) | (regrow_dises['area_acre'].isna()),
        np.nan,
        (
            (regrow_dises['area_acre'] >= area_match_lower_bound * regrow_dises['field_size_dises']) &
            (regrow_dises['area_acre'] <= area_match_upper_bound * regrow_dises['field_size_dises'])
        ).astype(int)
    )
    
    # Drop temporary column in regrow_dises
    regrow_dises.drop(columns="field_crop_23_dises_mapped", inplace=True)

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
    
    # Collect all merged regrow_dises together
    if regrow_dises_concat.empty:
        regrow_dises_concat = regrow_dises.copy()
    else:
        regrow_dises_concat = pd.concat([regrow_dises_concat, regrow_dises], axis=0, ignore_index=True)
        
    # Keep only the necessary columns
    cols_to_keep = [c for c in regrow_dises.columns if c.endswith("_dises") or c in ["field_id", "geometry"]]
    regrow_dises = regrow_dises[cols_to_keep]
    
    # Preview result
    print(regrow_dises.head())

    # Save spatial joined Regrow_dises shape file
    #regrow_dises.to_parquet(output_path_geospatial, compression="zstd")

    # Save attribute table as parquet
    attribute_table = regrow_dises.drop(columns='geometry')
    attribute_table.to_parquet(output_path_table, compression="zstd")
    print(f"Regrow and DISES for {state} are merged and saved")


# Add representative field attribute to Regrow fields
def assign_representative_field(regrow_dises):
    
    OVERLAP_THRESHOLD = 0.5
    CROP_CONF_THRESHOLD = 75
    
    mask_field_assigned = regrow_dises['field_assigned_dises'] == 'Y'
    mask_survey_responded = regrow_dises['survey_responded_dises'] == 'Y'
    mask_overlap_area = regrow_dises["overlap_area_share_dises"] > OVERLAP_THRESHOLD
    
    regrow_dises_matching = regrow_dises[(mask_field_assigned) & (mask_survey_responded) & (mask_overlap_area)].reset_index(drop=True)
    
    crop_conf_col = regrow_dises_matching[['crop_conf_23_1', 'crop_conf_23_2']].min(axis=1)
    
    regrow_dises_matching["RF_level_1_dises"] = (
        # A: size + crop match, high confidence
        (regrow_dises_matching["match_quality_dises"] == "A") &
        (crop_conf_col > CROP_CONF_THRESHOLD)
    ).replace(False, np.nan)
    
    regrow_dises_matching["RF_level_2_dises"] = (
        # A: size + crop match, low confidence
        (
            (regrow_dises_matching["match_quality_dises"] == "A") &
            (crop_conf_col <= CROP_CONF_THRESHOLD)
        ) |
        # B_area: size match only, crop missing
        (
            (regrow_dises_matching["match_quality_dises"] == "B_area") &
            (regrow_dises_matching["field_crop_23_dises"].isna())
        ) |
        # B_crop: crop match only, size missing, high confidence
        (
            (regrow_dises_matching["match_quality_dises"] == "B_crop") &
            (regrow_dises_matching["field_size_dises"].isna()) &
            (crop_conf_col > CROP_CONF_THRESHOLD)
        )
    ).replace(False, np.nan)
    
    regrow_dises_matching["RF_level_3_dises"] = (
        # Crop match only, size missing, low confidence
        (
            (regrow_dises_matching["match_quality_dises"] == "B_crop") &
            (regrow_dises_matching["field_size_dises"].isna()) &
            (crop_conf_col <= CROP_CONF_THRESHOLD)
        ) |
        # Size match only, crop mismatch, low confidence
        (
            (regrow_dises_matching["match_quality_dises"] == "B_area") &
            (regrow_dises_matching["field_crop_23_dises"].notna()) &
            (crop_conf_col <= CROP_CONF_THRESHOLD)
        )
    ).replace(False, np.nan)
    
    regrow_dises_matching["RF_level_4_dises"] = (
        # RF size missing, crop mismatch, low confidence
        (
            (regrow_dises_matching["match_quality_dises"] == "F") &
            (regrow_dises_matching["field_size_dises"].isna()) &
            (regrow_dises_matching["field_crop_23_dises"].notna()) &
            (crop_conf_col <= CROP_CONF_THRESHOLD)
        ) |
        # Both RF size and crop missing, only single parcel
        (
            (regrow_dises_matching["field_crop_23_dises"].isna()) &
            (regrow_dises_matching["field_size_dises"].isna()) &
            (regrow_dises_matching["n_parcels_dises"] == 1)
        )
    ).replace(False, np.nan)

    regrow_dises_matching["area_diff"] = (regrow_dises_matching["area_acre"] - regrow_dises_matching["field_size_dises"]).abs()
    regrow_dises_matching = regrow_dises_matching.sort_values(
        by=(["comp_id_dises"] + regrow_dises_matching.filter(regex="RF_level_").columns.tolist() + ["area_diff", "overlap_area_share_dises"]),
        ascending=[True, False, False, False, False, True, False]
    )
    
    regrow_dises_matching = regrow_dises_matching.drop_duplicates(subset=["comp_id_dises"], keep="first")
    regrow_dises_matching.reset_index(drop=True, inplace=True)
    
    col_to_keep = regrow_dises_matching.filter(regex="field_id|RF_level_").columns
    return(regrow_dises_matching[col_to_keep])


representative_field_indicator = assign_representative_field(regrow_dises_concat)
del regrow_dises_concat

for state in states:
    output_path_table = os.path.join(regrow_DISES_output_folder, f"{state}_regrow_dises_table.parquet")
    regrow_dises = pd.read_parquet(output_path_table)
    
    regrow_dises = regrow_dises.merge(representative_field_indicator, on='field_id', how='left')
    
    regrow_dises.to_parquet(output_path_table, compression="zstd")
    print(f"Regrow and DISES with RF indicators for {state} are saved")