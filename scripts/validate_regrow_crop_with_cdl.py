import pandas as pd
import os
from pathlib import Path
from rasterstats import zonal_stats
from concurrent.futures import ProcessPoolExecutor
from shapely.geometry import mapping

# Import parameters from the Snakemake
regrow_input_folder = snakemake.params.regrow_input_dir
regrow_validation_input_folder = snakemake.params.regrow_validation_dir
regrow_validation_output_folder = snakemake.params.regrow_validation_dir
states = snakemake.params.states
years_range = snakemake.params.years 


# Define land use values columns, validity check columns, and summary columns
def summarize(regrow_table_cdl, year, filter_cdl=None):
    
    cdl_col, valid_col = f"cdl_crop_{year}", f"cdl_valid_{year}"
    regrow_crop_cols = regrow_table_cdl.filter(regex=f"^crop_{str(year)[2:]}_*").columns
    # Drop all raws with Rergow main crop value of 999 for a given year since 999 indicates all types of non-cropland alnd
    regrow_table_cdl = regrow_table_cdl[~(regrow_table_cdl[regrow_crop_cols].eq(999).any(axis=1))]
    df = regrow_table_cdl[regrow_table_cdl[regrow_crop_cols].isin(filter_cdl).any(axis=1)] if filter_cdl else regrow_table_cdl
    valid = int(df[valid_col].sum())
    total = int(df[valid_col].count())
    return {
        "year": year,
        "cdl_category": filter_cdl,
        "total_fields": total,
        "valid_fields": valid,
        "invalid_fields": total - valid,
        "percent_valid_fields": (valid / total) if total > 0 else 0
    }

# Main script
# Create an output directory
os.makedirs(Path(regrow_validation_output_folder), exist_ok = True)

for state in states:

    regrow_cdl_input_path = os.path.join(regrow_validation_input_folder, f"{state}_regrow_cdl_validation_table.parquet")
    regrow_table_input_path = os.path.join(regrow_input_folder, f"{state}_regrow_table.parquet")
    output_path_xlsx = os.path.join(regrow_validation_output_folder, f"{state}_regrow_cdl_summary_by_crop_category.xlsx")
    
    regrow_cdl = pd.read_parquet(regrow_cdl_input_path)
    regrow_table = pd.read_parquet(regrow_table_input_path)
    
    regrow_table_cdl = regrow_table.merge(regrow_cdl, on='field_id', how='left')
    
    configs = [
        ("all", None),
        ("cdl_1", [1]),
        ("cdl_5", [5]),
        ("cdl_24", [24]),
        ("cdl_1_5_22_23_24", [1, 5, 22, 23, 24])]

    summaries = {name: [] for name, _ in configs}

    for year in years_range:
        print(f"Calculating matching indicators for {state} and year {year}...")
        for name, filter_cdl in configs:
            kwargs = {"filter_cdl": filter_cdl} if filter_cdl is not None else {}
            summaries[name].append(summarize(regrow_table_cdl, year, **kwargs))

    summary_cdl_by_category = [s for name, _ in configs for s in summaries[name]]

    # Save summary tables for validation
    pd.DataFrame(summary_cdl_by_category).to_excel(output_path_xlsx, index = False)
    print(f"Summary tables for {state} are created and saved")

