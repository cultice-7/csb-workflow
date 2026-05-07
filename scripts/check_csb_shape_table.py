import geopandas as gpd
import pandas as pd
from pathlib import Path
import os

# Import parameters from Snakemake
CSB_input_folder = snakemake.params.CSB_input_dir
CSB_checks_output_folder = snakemake.params.CSB_checks_dir
states = snakemake.params.states
CSB_years = snakemake.params.CSB_years
overlap_share_threshold = snakemake.params.overlap_share_threshold


### Clean CSB: CSB geometries contain overlaps which harms rasterization, we need to remove overlaps ###
def extract_duplicate_fields(state, CSB_year, overlap_share_threshold, CSB_input_folder, CSB_checks_output_folder):
    
    # Path to input and output files
    csb_geometry_input_path = os.path.join(CSB_input_folder, f"{state}_CSB{CSB_year}_CSBID_geometry.parquet")
    overlapping_fields_output_path = os.path.join(CSB_checks_output_folder, f"{state}_CSB{CSB_year}_overlapping_fields.parquet")
    
    # Check whether output path exists, if not, create it
    if not Path(overlapping_fields_output_path).parent.exists():
        Path(overlapping_fields_output_path).parent.mkdir(parents=True, exist_ok=True)
    
    # Load CSB shape (geometry) data
    csb_geometry = gpd.read_parquet(csb_geometry_input_path)

    csb_geometry["area"] = csb_geometry.geometry.area

    # Overlay CSBID on itself. This will produce all overlapping polygons (including boundaries)
    CSB_overlaps = gpd.overlay(csb_geometry, csb_geometry, how="intersection")

    # Remove self-overlaps
    CSB_overlaps = CSB_overlaps[CSB_overlaps["CSBID_1"] != CSB_overlaps["CSBID_2"]]

    # Compute overlap ratio
    # Determine the smaller area in each pair
    CSB_overlaps["min_area"] = CSB_overlaps[["area_1", "area_2"]].min(axis=1)

    # Compute overlap fraction relative to smaller polygon
    CSB_overlaps["overlap_fraction"] = CSB_overlaps.geometry.area / CSB_overlaps["min_area"]

    # Keep only those with overlap >= 50%
    CSB_overlaps = CSB_overlaps[CSB_overlaps["overlap_fraction"] >= overlap_share_threshold]

    # Keep only unique pairs (A,B) same as (B,A)
    # Make a tuple of (smaller_area_parcel, larger_area_parcel)
    CSB_overlaps["pair"] = CSB_overlaps.apply(
        lambda row: (row["CSBID_1"], row["CSBID_2"]) if row["area_1"] <= row["area_2"] else (row["CSBID_2"], row["CSBID_1"]), axis=1)
    # Drop duplicates
    CSB_overlaps = CSB_overlaps.drop_duplicates(subset="pair").reset_index(drop=True)
    # Assign smaller field IDs in each pair to CSBID_1 column and larger field IDs to CSBID_2
    CSB_overlaps.loc[:, 'CSBID_1'] = CSB_overlaps["pair"].apply(lambda x: x[0] if pd.notna(x) else None)
    CSB_overlaps.loc[:, 'CSBID_2'] = CSB_overlaps["pair"].apply(lambda x: x[1] if pd.notna(x) else None)
    # Rename CSBID_1 and CSBID_2 accordingly
    CSB_overlaps.rename(columns={"CSBID_1": "CSBID_smaller", "CSBID_2": "CSBID_larger"}, inplace = True)
    # Keep only necessary columns
    CSB_overlaps = CSB_overlaps[['CSBID_smaller', 'CSBID_larger', 'overlap_fraction']]
    # Save a separate dataset with overlapping fields
    CSB_overlaps.to_parquet(overlapping_fields_output_path, compression="zstd")
    
# Main code
for CSB_year in CSB_years:
    for state in states:
        extract_duplicate_fields(state, CSB_year, overlap_share_threshold, CSB_input_folder, CSB_checks_output_folder)