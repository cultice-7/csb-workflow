import pandas as pd
import os

# Input and output paths
input_folder_dises_table = snakemake.params.input_dir
output_folder_dises_table = snakemake.params.output_dir

input_path_dises_table = os.path.join(input_folder_dises_table, "combined_data_clean.csv")
output_path_dises_table = os.path.join(output_folder_dises_table, "DISES_table_cleaned.csv")

# Read DISES table
dises_table = pd.read_csv(input_path_dises_table)

# Delete entires where where all values are missing:
dises_table_cleaned = dises_table.dropna(how="all", axis=0).reset_index(drop=True)

# Delete entries with missing comprehensive ID
dises_table_cleaned = dises_table_cleaned[dises_table_cleaned['Comprehensive_ID'].notna()].reset_index(drop=True)

# Rename column Comprehensive_ID to comp_id, field_crop to field_crop_23, field_tillage to field_till_23, field_CC to field_cover_23  
dises_table_cleaned.rename(columns={'Comprehensive_ID': 'comp_id',
                         'field_crop': 'field_crop_23',
                         'field_tillage': 'field_till_23',
                         'field_CC': 'field_cover_23'}, inplace = True)

# Compute farmer types using FI answers from the survey
# Productivism index
dises_table_cleaned['productivism_index'] = dises_table_cleaned[['FI_1', 'FI_4', 'FI_7', 'FI_9', 'FI_10']].mean(axis=1)
# Conservationism index
dises_table_cleaned['conservationism_index'] = dises_table_cleaned[['FI_2', 'FI_3', 'FI_5', 'FI_6', 'FI_8']].mean(axis=1)
# Civic index
dises_table_cleaned['civic_index'] = dises_table_cleaned[['FI_11', 'FI_12', 'FI_13', 'FI_14', 'FI_15', 'FI_16']].mean(axis=1)

# Save dises data table for future join with dises shape
os.makedirs(input_folder_dises_table, exist_ok = True)
dises_table_cleaned.to_csv(output_path_dises_table, index = False) 