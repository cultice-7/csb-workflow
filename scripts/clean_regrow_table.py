import pandas as pd
import os
import numpy as np

# Import parameters from Snakemake
regrow_table_input_folder = snakemake.params.input_dir
regrow_table_output_folder = snakemake.params.output_dir
states = snakemake.params.states


def clean_regrow_table(state, regrow_table_input_folder, regrow_table_output_folder):
    
    regrow_table_input_path = os.path.join(regrow_table_input_folder, f"{state}_Monitor_data_2014-2025.parquet")
    regrow_table_output_path = os.path.join(regrow_table_output_folder, f"{state}_regrow_monitor_wide_coded.parquet")
    
    df = pd.read_parquet(regrow_table_input_path)
    
    print("State name:", state)
    # Drop completely duplicate rows
    df.drop_duplicates(keep = 'first', inplace = True)
     
    # Drop rows with completely missing monitor data
    cols_to_check_nan = df.columns.difference(['field_id'])
    df = df.dropna(subset=cols_to_check_nan, how='all')

    df.reset_index(drop = True, inplace = True)

    # Rename columns
    df.rename(columns={'boundary_id': 'field_id',
                    'PostHarvest_Intensity': 'PHtill', 
                    'PostHarvest_till_start': 'PHtill_start',
                    'PostHarvest_till_end': 'PHtill_end',
                    'PostHarvest_till_residue_conf': 'PHtill_conf',
                    'Cover_Crop_Class': 'cover',
                    'Cover_Crop_Start_Date': 'cover_start',
                    'Cover_Crop_End_Date': 'cover_end',
                    'Cover_Crop_Confidence': 'cover_conf',
                    'PrePlant_Intensity': 'PPtill',
                    'PrePlant_till_start': 'PPtill_start',
                    'PrePlant_till_end': 'PPtill_end',
                    'PrePlant_till_residue_conf': 'PPtill_conf',
                    'crop': 'crop', 
                    'crop_plant_date': 'crop_start',
                    'crop_harvest_date': 'crop_end',
                    'crop_confidence': 'crop_conf',
                    }, inplace=True)

    # Change the format of date columns to datetime
    cols_to_convert = ['cycle_start', 'cycle_end', 'PHtill_start', 'PHtill_end',
                    'cover_start', 'cover_end', 'PPtill_start', 'PPtill_end', 'crop_start', 'crop_end']
    df[cols_to_convert] = df[cols_to_convert].apply(pd.to_datetime, errors = 'coerce')


    # Aggregate multiple activities within one cycle together
    # Columns determining cultivation cycle
    key_cols = ['field_id', 'crop', 'crop_start', 'crop_end', 'crop_conf']

    # List pre-plant till, post-harvest till and cover crop columns (aggregate their values if there are multiple acitivities within a single cycle)
    activity_cols = ["PPtill", "PHtill", "cover"]

    # List corresponding date and confidence columns for activity_cols (aggregate their values if there are multiple acitivities within a single cycle)
    numeric_ops = ["start", "end", "conf"]
    numeric_cols = [f"{col}_{op}" for col in activity_cols for op in numeric_ops]

    # List corresponding cycle_start/end columns (aggregate their values if there are multiple acitivities within a signle cycle)
    numeric_cols += ["cycle_start", "cycle_end"]

    # Split rows with duplicate main crops into a separate dataset
    dup_mask = df.duplicated(subset=key_cols, keep=False)
    df_duplicates = df[dup_mask].copy()
    df_non_duplicates = df[~dup_mask].copy()

    # Functions to aggregate land management activities with multiple values within a single cultivation cycle
    def aggregate_activity_values(values, start_dates, end_dates):
        """
        - If all start/end dates identical and values identical -> keep value
        - If all start/end identical but values differ -> NaN
        - If dates differ -> join all activity values
        - Preserve NaNs
        """
        mask = values.notna()
        values = values[mask]
        start_dates = start_dates[mask]
        end_dates = end_dates[mask]

        if len(values) == 0:
            return np.nan

        if (start_dates.nunique() == 1) and (end_dates.nunique() == 1):
            if values.nunique() == 1:
                return values.iloc[0]
            else:
                return np.nan

        return ' & '.join(values.astype(str))

    def aggregate_group(g):
        result = {}

        # Activity columns
        for col in activity_cols:
            # Compute start/end date and conf value of aggregated activities
            result[f"{col}_start"] = g[f"{col}_start"].min(skipna=True)
            result[f"{col}_end"] = g[f"{col}_end"].max(skipna=True)
            result[f"{col}_conf"] = g[f"{col}_conf"].max(skipna=True)

            # Aggregate activity values according to dates
            result[col] = aggregate_activity_values(g[col], g[f"{col}_start"], g[f"{col}_end"])

        # Compute start/end dates of cultivation cycles containing aggregated activities
        result["cycle_start"] = g["cycle_start"].min(skipna=True)
        result["cycle_end"] = g["cycle_end"].max(skipna=True)

        return pd.Series(result)

    # Process entries with multiple land management activities within a single cultivation cycle
    df_processed = df_duplicates.groupby(key_cols).apply(aggregate_group).reset_index()

    # Merge back with non-duplicate rows
    df = pd.concat([df_non_duplicates, df_processed], ignore_index=True)


    # Extract starting year of each cycle and compute cycle length
    df['start_year'] = df['cycle_start'].dt.year
    df['end_year'] = df['cycle_end'].dt.year
    df['cycle_length'] = (df['cycle_end'] - df['cycle_start']).dt.days

    # Sort dataset by field_id, start_year, plant_date (so that earlier cultivation cycles within same year come first)
    df.sort_values(by = ['field_id', 'start_year', 'cycle_start'], inplace = True)
    df.reset_index(drop = True, inplace = True)

    # Compute the number of cultivation cycles for each field in each year
    df['cycle_count'] = df.groupby(['field_id', 'end_year']).cumcount() + 1

    # Keep observations with not more than 3 cultivation cycles per calendar year
    df = df[df['cycle_count'] <= 3]
    df.reset_index(drop = True, inplace = True)

    # Create calendar year_cycle identifier
    df['year_cycle'] = df['end_year'].astype(str).str[2:] + '_' + df['cycle_count'].astype(str)

    # Reshape post harvest tillage from long to wide
    df_PHtill_wide = df.pivot(index = 'field_id', columns = 'year_cycle', values = ['PHtill', 'PHtill_start', 'PHtill_end', 'PHtill_conf'])
    df_PHtill_wide.reset_index(inplace = True)
    df_PHtill_wide.columns = ['field_id'] + ['{}_{}'.format(col[0], col[1]) for col in df_PHtill_wide.columns[1:]]

    # Reshape cover crop from long to wide
    df_cover_wide = df.pivot(index = 'field_id', columns = 'year_cycle', values = ['cover', 'cover_start', 'cover_end', 'cover_conf'])
    df_cover_wide.reset_index(inplace = True)
    df_cover_wide.columns = ['field_id'] + ['{}_{}'.format(col[0], col[1]) for col in df_cover_wide.columns[1:]]

    # Reshape pre-plant tillage from long to wide
    df_PPtill_wide = df.pivot(index = 'field_id', columns = 'year_cycle', values = ['PPtill', 'PPtill_start', 'PPtill_end', 'PPtill_conf']) 
    df_PPtill_wide.reset_index(inplace = True)
    df_PPtill_wide.columns = ['field_id'] + ['{}_{}'.format(col[0], col[1]) for col in df_PPtill_wide.columns[1:]]

    # Reshape main crop from long to wide
    df_crop_wide = df.pivot(index = 'field_id', columns = 'year_cycle', values = ['crop', 'crop_start', 'crop_end', 'crop_conf'])
    df_crop_wide.reset_index(inplace = True)
    df_crop_wide.columns = ['field_id'] + ['{}_{}'.format(col[0], col[1]) for col in df_crop_wide.columns[1:]]

    # Recode main crop names to numerical codes 
    main_crop_map = {
        "corn": 1, "cotton": 2, "rice": 3, "sorghum": 4, "soybean": 5, "sunflower": 6,
        "peanut": 10, "tobacco": 11, "sweet_corn": 12, "pop_or_orn_corn": 13, "mint": 14,
        "barley": 21, "wheat_durum": 22, "wheat_spring": 23, "wheat_winter": 24, "rye": 27, "oat": 28,
        "millet": 29, "speltz": 30, "canola": 31, "flaxseed": 32, "flax": 32, "safflower": 33,
        "mustard": 35, "alfalfa": 36, "hay": 37, "ryegrass": 37, "camelina": 38, "buckwheat": 39,
        "sugar_beet": 41, "dry_bean": 42, "potato": 43, "other": 44, "sugarcane": 45, "sweet_potato": 46,
        "vegetable": 47, "watermelon": 48, "onion": 49, "cucumber": 50, "lentil": 52, "pea": 53, "tomato": 54, "herb": 57,
        "clover": 58, "sod_grass": 59, "switchgrass": 60, "fallow": 61,
        "cherry": 66, "peach": 67, "apple": 68, "grape": 69, "pecan": 74, "walnut": 76,
        "forest_deciduous": 141, "forest_evergreen": 142, "evergreen": 142, "shrub": 152,  "grass_perennial": 176,
        "pistachio": 204, "triticale": 205, "carrot": 206, "canteloupe": 209, "broccoli": 214, "pepper": 216, "greens": 219, "strawberry": 221,
        "squash": 222, "vetch": 224, "pumpkin": 229, "blueberry": 242, "cabbage": 243, "cauliflower": 244, "radish": 246, 
        "turnip": 247, "cranberry": 250, "non_cropland": 999, "berry": 999
    }  

    df_crop_wide.replace(main_crop_map, inplace = True)
    # Select columns containing "crop_1" or "crop_2"
    cols_to_convert = df_crop_wide.columns[df_crop_wide.columns.str.startswith(("crop_1", "crop_2"))]
    # Convert to numeric and then to nullable Int16
    df_crop_wide[cols_to_convert] = (df_crop_wide[cols_to_convert].apply(pd.to_numeric, errors="coerce").astype("Int16"))

    # Recode tillage names to numerical codes 
    tillage_map = {
        "TILLAGE_INTENSITY_NO_TILLAGE": 1, 
        "TILLAGE_INTENSITY_REDUCED_TILLAGE": 2, 
        "TILLAGE_INTENSITY_CONVENTIONAL_TILLAGE": 3,
        "TILLAGE_INTENSITY_NO_TILLAGE & TILLAGE_INTENSITY_NO_TILLAGE": 1,
        "TILLAGE_INTENSITY_NO_TILLAGE & TILLAGE_INTENSITY_REDUCED_TILLAGE": 2, "TILLAGE_INTENSITY_REDUCED_TILLAGE & TILLAGE_INTENSITY_NO_TILLAGE": 2,
        "TILLAGE_INTENSITY_NO_TILLAGE & TILLAGE_INTENSITY_CONVENTIONAL_TILLAGE": 3, "TILLAGE_INTENSITY_CONVENTIONAL_TILLAGE & TILLAGE_INTENSITY_NO_TILLAGE": 3,
        "TILLAGE_INTENSITY_REDUCED_TILLAGE & TILLAGE_INTENSITY_REDUCED_TILLAGE": 4,
        "TILLAGE_INTENSITY_REDUCED_TILLAGE & TILLAGE_INTENSITY_CONVENTIONAL_TILLAGE": 5, "TILLAGE_INTENSITY_CONVENTIONAL_TILLAGE & TILLAGE_INTENSITY_REDUCED_TILLAGE": 5,
        "TILLAGE_INTENSITY_CONVENTIONAL_TILLAGE & TILLAGE_INTENSITY_CONVENTIONAL_TILLAGE": 6,
        "TILLAGE_INTENSITY_NO_TILLAGE_DATA": np.nan
    }

    df_PHtill_wide.replace(tillage_map, inplace = True)
    # Select columns containing "PHtill_1" or "PHtill_2"
    cols_to_convert = df_PHtill_wide.columns[df_PHtill_wide.columns.str.startswith(("PHtill_1", "PHtill_2"))]
    # Convert to numeric and then to nullable Int16
    df_PHtill_wide[cols_to_convert] = (df_PHtill_wide[cols_to_convert].apply(pd.to_numeric, errors="coerce").astype("Int16"))


    df_PPtill_wide.replace(tillage_map, inplace = True)
    # Select columns containing "PPtill_1" or "PPtill_2"
    cols_to_convert = df_PPtill_wide.columns[df_PPtill_wide.columns.str.startswith(("PPtill_1", "PPtill_2"))]
    # Convert to numeric and then to nullable Int16
    df_PPtill_wide[cols_to_convert] = (df_PPtill_wide[cols_to_convert].apply(pd.to_numeric, errors="coerce").astype("Int16"))


    # Recode cover crop names to numerical codes 
    cover_crop_map = {
        "GREEN_COVER_CLASS_COVER_CROP_NOT_TRACKED": 1, 
        "GREEN_COVER_CLASS_POTENTIAL_COVER_CROP": 2, 
        "GREEN_COVER_CLASS_COVER_CROP": 3,
        "GREEN_COVER_CLASS_NO_DATA": np.nan,
        "GREEN_COVER_CLASS_NOT_APPLICABLE": 999
    }

    df_cover_wide.replace(cover_crop_map, inplace = True)
    # Select columns containing "cover_1" or "cover_2"
    cols_to_convert = df_cover_wide.columns[df_cover_wide.columns.str.startswith(("cover_1", "cover_2"))]
    # Convert to numeric and then to nullable Int16
    df_cover_wide[cols_to_convert] = (df_cover_wide[cols_to_convert].apply(pd.to_numeric, errors="coerce").astype("Int16"))


    # Join all tables into one
    regrow_cleaned_table = df_PHtill_wide.merge(df_cover_wide, on = 'field_id', how = 'left')
    del df_PHtill_wide, df_cover_wide # remove already merged datasets from memory
    regrow_cleaned_table = regrow_cleaned_table.merge(df_PPtill_wide, on = 'field_id', how = 'left')
    del df_PPtill_wide # remove already merged datasets from memory
    regrow_cleaned_table = regrow_cleaned_table.merge(df_crop_wide, on = 'field_id', how = 'left')
    del df_crop_wide # remove already merged datasets from memory

    # Drop all-NaN columns
    regrow_cleaned_table =  regrow_cleaned_table.dropna(axis=1, how='all')

    # Save to csv
    os.makedirs(regrow_table_output_folder, exist_ok = True)
    regrow_cleaned_table.to_parquet(regrow_table_output_path, compression="zstd")
    
    return regrow_cleaned_table.columns, regrow_table_output_path


# Main script
# After cleaning, unify column list in all cleaned regrow table files to ensure consistency 
all_cols = set() 
regrow_cleaned_table_paths = []

for state in states:
    col_list, output_path = clean_regrow_table(state, regrow_table_input_folder, regrow_table_output_folder)
    all_cols.update(col_list)
    regrow_cleaned_table_paths.append(output_path)
   
 # Sort column names to ensure consistency
id_col = "field_id"
all_cols = [id_col] + sorted(c for c in all_cols if c != id_col)

# Ensure each clean_regrow_table has all columns
for regrow_cleaned_table_path in regrow_cleaned_table_paths:
    
    regrow_cleaned_table = pd.read_parquet(regrow_cleaned_table_path)
    
    missing_cols = set(all_cols) - set(regrow_cleaned_table.columns)
    
    # Add missing columns with a dtype matching what the column would have if it were present:
    # bare "{activity}_{year}_{cycle}" code columns -> nullable Int16
    # "{activity}_conf_{year}_{cycle}" -> float64
    # "{activity}_start/end_{year}_{cycle}" -> datetime64
    activity_names = ("PHtill", "PPtill", "cover", "crop")
    for col in missing_cols:
        parts = col.split("_")
        if len(parts) == 3 and parts[0] in activity_names and len(parts[1]) == 2:
            # e.g. "PHtill_19_2"
            regrow_cleaned_table[col] = pd.Series(pd.NA, index=regrow_cleaned_table.index, dtype="Int16")
        elif len(parts) == 4 and parts[0] in activity_names and parts[1] == "conf" and len(parts[2]) == 2:
            # e.g. "PHtill_conf_19_2"
            regrow_cleaned_table[col] = np.nan
        elif len(parts) == 4 and parts[0] in activity_names and parts[1] in ("start", "end") and len(parts[2]) == 2:
            # e.g. "PHtill_start_19_2" / "PHtill_end_19_2"
            regrow_cleaned_table[col] = pd.Series(pd.NaT, index=regrow_cleaned_table.index, dtype="datetime64[ns]")
        else:
            # Fallback for any unexpected column name
            regrow_cleaned_table[col] = np.nan
    
    # Reorder columns
    regrow_cleaned_table = regrow_cleaned_table[all_cols]
    
    regrow_cleaned_table.to_parquet(regrow_cleaned_table_path, compression="zstd")
