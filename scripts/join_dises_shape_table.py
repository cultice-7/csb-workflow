import geopandas as gpd
import pandas as pd
import os


# Input and output paths
input_folder_dises = snakemake.params.input_dir
output_folder_dises = snakemake.params.output_dir

input_path_dises_shape = os.path.join(input_folder_dises, "DISES_shape_cleaned.parquet")
input_path_dises_table = os.path.join(input_folder_dises, "DISES_table_cleaned.csv")
output_path_dises_shape_table = os.path.join(output_folder_dises, "DISES_shape_table.parquet")

# Read DISES shape file
dises_shape = gpd.read_parquet(input_path_dises_shape)

# Load dises data table
dises_table = pd.read_csv(input_path_dises_table)

# Add indication that this comp_id exist in data table, to be shown in joined shape file
dises_table['survey_responded'] = "Y"

"""
# Replace field_crop values to match CDL classification (We'll do the same with Regrow)
dises_table['field_crop_23'] = dises_table['field_crop_23'].replace(2,5)

# Replace field_CC values to match Regrow
dises_table['field_cover_23'] = dises_table['field_cover_23'].replace(1,3) #1-Yes to 3-confirmed CC
dises_table['field_cover_23'] = dises_table['field_cover_23'].replace(2,1) #2-No to 1-No cC
"""

# Join shape file and table using comprehensive ID
dises_shape_table = dises_shape.merge(dises_table, on='comp_id', how='left')

# Replace nan values in the "survey_responded" column witn "N" value
dises_shape_table["survey_responded"] = dises_shape_table["survey_responded"].fillna("N")

# Save joined DISES shape and table files to parquet
dises_shape_table.to_parquet(output_path_dises_shape_table, compression="zstd")