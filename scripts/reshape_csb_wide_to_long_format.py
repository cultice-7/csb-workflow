import os
import re
import numpy as np
import pandas as pd
import pyarrow as pa
import pyarrow.parquet as pq

# Import parameters from Snakemake
states = snakemake.params.states
CSB_years = snakemake.params.CSB_years
CSB_input_folder = snakemake.params.CSB_input_dir
CSB_output_folder = snakemake.params.CSB_output_dir

list_of_datasets = ["census_tract_table", "elevation_slope_table", "watershed_table",
                    "nearest_roads_table", "neighbor_field_mgmt_table", "weather_table",
                    "crop_prices_table_reduced", "soil_composition_table"]

dynamic_vars = ['CDL', 'CDL_corn_nbr', 'CDL_soybean_nbr', 'CDL_wheat_nbr', 'CDL_other_nbr',
                'ppt_centroid', 'tmin_centroid', 'tmean_centroid', 'tdmean_centroid', 'soltotal_centroid',
                'corn_price_elevator_5-nearest', 'soybeans_price_elevator_5-nearest']
static_vars = ['CSBACRES',
               'county_state_name', 'watershed_name',
               'elevation_mean', 'slope_mean',
               'dist_to_road',
               'claytotal_r_30cm_weighted', 'sandtotal_r_30cm_weighted', 'ph1to1h2o_r_30cm_weighted', 'om_r_30cm_weighted', 'resdept_r_dominant']

id_var = ['CSBID']


def classify_column(col, valid_yyyy, valid_yy):
    parts = col.split('_')
    candidates = []  # (position, type, raw_suffix)

    for i, p in enumerate(parts):
        if re.fullmatch(r'\d{6}', p):
            yr, mo = int(p[:4]), int(p[4:])
            if yr in valid_yyyy and 1 <= mo <= 12:
                candidates.append((i, 'YYYYMM', p))
        elif re.fullmatch(r'\d{4}', p):
            yr = int(p)
            if yr in valid_yyyy:
                candidates.append((i, 'YYYY', p))
        elif re.fullmatch(r'\d{2}', p):
            yy = int(p)
            if yy in valid_yy:
                candidates.append((i, 'YY', p))

    if len(candidates) == 0:
        return col, 'static', None, None      # no time-like token found -> static
    if len(candidates) > 1:
        return col, 'ambiguous', None, None     # more than one plausible token -> flag for manual check

    idx, ptype, suffix = candidates[0]
    stub_parts = parts[:idx] + parts[idx+1:]    # remove the time token, keep everything else in order
    stub = '_'.join(stub_parts)
    return stub, ptype, suffix, idx


def normalize_columns(df, valid_yyyy, valid_yy):
    rename_map = {}
    review = []

    for col in df.columns:
        stub, ptype, suffix, idx = classify_column(col, valid_yyyy, valid_yy)

        if ptype == 'YYYYMM':
            year, month = suffix[:4], suffix[4:]
            rename_map[col] = f"{stub}_{month}_{year}"   # month folded into stub, year as final suffix
        elif ptype == 'YY':
            yr = int(suffix)
            yyyy = 2000 + yr
            rename_map[col] = f"{stub}_{yyyy}"
        elif ptype == 'YYYY':
            rename_map[col] = f"{stub}_{suffix}"
        elif ptype == 'ambiguous':
            review.append(col)
        # static columns: leave unrenamed

    if review:
        print("Columns needing manual review (ambiguous time token):")
        print(review)

    return df.rename(columns=rename_map)


year_pattern = re.compile(r'^(.+)_(\d{4})$')
chunk_size = 20000

os.makedirs(CSB_output_folder, exist_ok=True)

for year in CSB_years:
    # CSB_years is a release-window label like "1724", covering 2017-2024
    MIN_YEAR, MAX_YEAR = 2000 + int(year[:2]), 2000 + int(year[2:])
    VALID_YY = {y % 100 for y in range(MIN_YEAR, MAX_YEAR + 1)}
    VALID_YYYY = set(range(MIN_YEAR, MAX_YEAR + 1))

    for state in states:

        # Read the base CSB table and merge in every supplement block
        csb_df = pd.read_parquet(os.path.join(CSB_input_folder, f"{state}_CSB{year}_table.parquet"))
        for dataset_name in list_of_datasets:
            df_temp = pd.read_parquet(os.path.join(CSB_input_folder, f"{state}_CSB{year}_{dataset_name}.parquet"))
            csb_df = csb_df.merge(df_temp, how="left", on="CSBID")

        dynamic_cols = csb_df.columns.str.startswith(tuple(dynamic_vars))
        static_cols = csb_df.columns.str.startswith(tuple(static_vars))

        csb_df = normalize_columns(csb_df, VALID_YYYY, VALID_YY)

        stub_cols = [c for c in csb_df.columns if year_pattern.match(c)]
        stubnames = sorted(set(year_pattern.match(c).group(1) for c in stub_cols))

        # Reshape from wide to long in row chunks, writing incrementally to avoid holding
        # the full long-format table (many more rows than the wide table) in memory at once
        output_path = os.path.join(CSB_output_folder, f"{state}_CSB{year}_supplement_long_table.parquet")
        writer = None

        for start in range(0, len(csb_df), chunk_size):
            chunk = csb_df.iloc[start:start + chunk_size]

            long_chunk = pd.wide_to_long(
                chunk,
                stubnames=stubnames,
                i=id_var,
                j='year',
                sep='_',
                suffix=r'\d{4}'
            ).reset_index()

            table = pa.Table.from_pandas(long_chunk)
            if writer is None:
                writer = pq.ParquetWriter(output_path, table.schema)
            writer.write_table(table)

        writer.close()
        print(f"CSB{year} supplement long-format table for {state} is created and saved")
