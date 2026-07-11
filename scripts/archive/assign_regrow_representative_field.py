import geopandas as gpd
import pandas as pd
import numpy as np
import os
from shapely import intersection
from pathlib import Path

# Import parameters from Snakemake
regrow_input_folder = snakemake.params.regrow_input_dir
regrow_DISES_input_folder  = snakemake.params.regrow_DISES_input_dir
regrow_DISES_output_folder = snakemake.params.regrow_DISES_output_dir
states = snakemake.params.states


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
        # B_size: size match only, crop missing
        (
            (regrow_dises_matching["match_quality_dises"] == "B_size") &
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
            (regrow_dises_matching["match_quality_dises"] == "B_size") &
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

# Main code
regrow_dises_concat = pd.DataFrame()
for state in states:
    regrow_table_input_path = os.path.join(regrow_input_folder, f"{state}_regrow_table.parquet")
    regrow_dises_input_path = os.path.join(regrow_DISES_output_folder, f"{state}_regrow_dises_table.parquet")
    
    regrow_table = pd.read_parquet(regrow_table_input_path)
    regrow_dises = pd.read_parquet(regrow_dises_input_path)
    
    regrow_dises = regrow_table.merge(regrow_dises, how='field_id', on='left')
    
    if regrow_dises_concat.empty():
        regrow_dises_concat = regrow_dises.copy()
    else:
        regrow_dises_concat = pd.concat([regrow_dises_concat, regrow_dises], axis=0, ignore_index=True)


representative_field_indicator = assign_representative_field(regrow_dises_concat)


for state in states:
    regrow_dises_input_path = os.path.join(regrow_DISES_output_folder, f"{state}_regrow_dises_table.parquet")
    regrow_dises_output_path = os.path.join(regrow_DISES_output_folder, f"{state}_regrow_dises_table.parquet")
    
    regrow_dises = pd.read_parquet(regrow_dises_input_path)
    
    regrow_dises = regrow_dises.merge(representative_field_indicator, on='field_id', how='left')
    
    regrow_dises.to_parquet(regrow_dises_output_path, compression="zstd")
    print(f"Regrow and DISES for {state} are merged and saved")