import geopandas as gpd
import pandas as pd
from pathlib import Path
import os

# Import parameters from Snakemake
regrow_input_folder = snakemake.params.regrow_input_dir
regrow_checks_output_folder = snakemake.params.regrow_checks_dir
states = snakemake.params.states
overlap_share_threshold = snakemake.params.overlap_share_threshold


### Clean Regrow: Regrow geometries contain overlaps which harms rasterization, we need to remove overlaps ###
def extract_duplicate_fields(state, overlap_share_threshold, regrow_input_folder, regrow_checks_output_folder):
        
    # Path to input and output files
    regrow_geometry_input_path = os.path.join(regrow_input_folder, f"{state}_regrow_fieldID_geometry.parquet")
    overlapping_fields_output_path = os.path.join(regrow_checks_output_folder, f"{state}_regrow_overlapping_fields.parquet")
    
    # Check whether output path exists, if not, create it
    if not Path(overlapping_fields_output_path).parent.exists():
        Path(overlapping_fields_output_path).parent.mkdir(parents=True, exist_ok=True)
    
    # Load regrow shape (geometry) data
    regrow_geometry = gpd.read_parquet(regrow_geometry_input_path)

    regrow_geometry["area"] = regrow_geometry.geometry.area

    # Overlay Regrow on itself. This will produce all overlapping polygons (including boundaries)
    regrow_overlaps = gpd.overlay(regrow_geometry, regrow_geometry, how="intersection")

    # Remove self-overlaps
    regrow_overlaps = regrow_overlaps[regrow_overlaps["field_id_1"] != regrow_overlaps["field_id_2"]]

    # Compute overlap ratio
    # Determine the smaller area in each pair
    regrow_overlaps["min_area"] = regrow_overlaps[["area_1", "area_2"]].min(axis=1)

    # Compute overlap fraction relative to smaller polygon
    regrow_overlaps["overlap_fraction"] = regrow_overlaps.geometry.area / regrow_overlaps["min_area"]

    # Keep only those with overlap >= 50%
    regrow_overlaps = regrow_overlaps[regrow_overlaps["overlap_fraction"] >= overlap_share_threshold]

    # Keep only unique pairs (A,B) same as (B,A)
    # Make a tuple of (smaller_area_parcel, larger_area_parcel)
    regrow_overlaps["pair"] = regrow_overlaps.apply(
        lambda row: (row["field_id_1"], row["field_id_2"]) if row["area_1"] <= row["area_2"] else (row["field_id_2"], row["field_id_1"]), axis=1)
    # Drop duplicates
    regrow_overlaps = regrow_overlaps.drop_duplicates(subset="pair").reset_index(drop=True)
    # Assign smaller field IDs in each pair to field_id_1 column and larger field IDs to field_id_2
    regrow_overlaps.loc[:, 'field_id_1'] = regrow_overlaps["pair"].apply(lambda x: x[0] if pd.notna(x) else None)
    regrow_overlaps.loc[:, 'field_id_2'] = regrow_overlaps["pair"].apply(lambda x: x[1] if pd.notna(x) else None)
    # Rename field_id_1 and field_id_2 accordingly
    regrow_overlaps.rename(columns={"field_id_1": "field_id_smaller", "field_id_2": "field_id_larger"}, inplace = True)
    # Keep only necessary columns
    regrow_overlaps = regrow_overlaps[['field_id_smaller', 'field_id_larger', 'overlap_fraction']]
    # Save a separate dataset with overlapping fields
    regrow_overlaps.to_parquet(overlapping_fields_output_path, compression="zstd")


# Main code
for state in states:
    extract_duplicate_fields(state, overlap_share_threshold, regrow_input_folder, regrow_checks_output_folder)