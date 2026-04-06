import geopandas as gp
import pandas as pd
import os
from pathlib import Path
from rasterstats import zonal_stats
from concurrent.futures import ProcessPoolExecutor
from shapely.geometry import mapping

# Define process of zonal statistics and validity check
def process_polygon(args):
    geom, crop_values, raster_path = args
    try:
        stats = zonal_stats(geom, raster_path, stats="majority", nodata=0)
        majority = stats[0].get("majority", None) if stats else None
        valid = 1 if majority in crop_values else 0
        return majority, valid
    except Exception:
        return None, 0

# Define land use values columns, validity check columns, and summary columns
def summarize(gdf, year, filter_cdl=None):
    cdl_col, valid_col = f"cdl_{year}", f"valid_{year}"
    df = gdf[gdf[cdl_col].isin(filter_cdl)] if filter_cdl else gdf
    valid = int(df[valid_col].sum())
    total = int(df[valid_col].count())
    return {
        "year": year,
        "total_count": total,
        "valid_count": valid,
        "invalid_count": total - valid,
        "percent_valid": (valid / total) * 100 if total > 0 else 0
    }

# Define main processing loop
def main():

    # Import parameters from the Snakemake
    regrow_input_folder = snakemake.params.regrow_input_dir
    regrow_output_folder = snakemake.params.regrow_output_dir
    cdl_input_folder = snakemake.params.cdl_input_dir
    years_range = snakemake.params.years 
    
    # Create an output directory
    os.makedirs(Path(regrow_output_folder) / "validation", exist_ok = True)
    
    states = snakemake.params.states
    
    for state in states:
    
        input_path = os.path.join(regrow_input_folder, f"{state}_regrow_shape_table.parquet")
        output_path = os.path.join(regrow_output_folder, "validation", f"{state}_regrow_with_cdl_validation.parquet")
        output_path_xlsx_1 = os.path.join(regrow_output_folder, "validation", f"{state}_regrow_validity_summary_by_year.xlsx")
        output_path_xlsx_2 = os.path.join(regrow_output_folder, "validation", f"{state}_regrow_validity_summary_by_year_cdl_1_5.xlsx")
        
        gdf = gp.read_file(input_path)
        summary_all, summary_cdl_1_5 = [], []

        for year in years_range:
            raster_path = os.path.join(cdl_input_folder, f"{year}_30m_cdls_clipped.tif")
            if not os.path.exists(raster_path):
                print(f"Raster not found for year {year}")
                continue

            print(f"Processing year {year} with {len(gdf)} polygons...")
            crop_cols = [f"crop{str(year)[2:]}_{i}" for i in range(1, 4) if f"crop{str(year)[2:]}_{i}" in gdf.columns]
            args_list = [(mapping(row.geometry), [row[c] for c in crop_cols if pd.notnull(row[c])], raster_path)
                        for _, row in gdf.iterrows()]

            with ProcessPoolExecutor() as executor:
                results = list(executor.map(process_polygon, args_list))

            gdf[f"cdl_{year}"], gdf[f"valid_{year}"] = zip(*results)
            summary_all.append(summarize(gdf, year))
            summary_cdl_1_5.append(summarize(gdf, year, filter_cdl=[1, 5]))

        # Save Regrow polygons with validation and summary tables for validation
        gdf.to_parquet(output_path, compression="zstd")
        pd.DataFrame(summary_all).to_excel(output_path_xlsx_1, index = False)
        pd.DataFrame(summary_cdl_1_5).to_excel(output_path_xlsx_2, index = False)
        print("Saved GeoJSON and both summary tables for {state}")

if __name__ == "__main__":
    main()
