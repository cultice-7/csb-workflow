import geopandas as gpd
import os

# Import variables from snakemake
input_folder_dises_shape = snakemake.params.input_dir
output_folder_dises_shape = snakemake.params.output_dir
target_CRS = snakemake.params.target_CRS

# Input and output paths
input_path_dises_shape = os.path.join(input_folder_dises_shape, "DISES_All_Parcels_11.12.25.shp")
output_path_dises_shape = os.path.join(output_folder_dises_shape, "DISES_shape_cleaned.parquet")

# Read DISES shape file
dises_shape = gpd.read_file(input_path_dises_shape)

# Delete comprehensive ID = 0
dises_shape_cleaned = dises_shape[dises_shape['comprehens'] != 0].copy()

# Rename comprehensive ID column from Comprehe_1 to Comp_ID
dises_shape_cleaned.rename(columns={'comprehens': 'comp_id'}, inplace=True)

# Keep only comp_id and geometry
dises_shape_cleaned = dises_shape_cleaned[['comp_id', 'geometry']]

# Consolidate multiple tax parcels into one for each owner
dises_shape_cleaned_consolidated = dises_shape_cleaned.dissolve(by='comp_id')

# Reset index so comp_id becomes a column again
dises_shape_cleaned_consolidated.reset_index(inplace=True)

# Reproject to NAD83/Conus Albers
dises_shape_cleaned_consolidated = dises_shape_cleaned_consolidated.to_crs(target_CRS)

#  Calculate area in US survey acres (1 acre = 4046.8564224 sq m)
dises_shape_cleaned_consolidated["area_acre"] = dises_shape_cleaned_consolidated.geometry.area/4046.8564224

# Save shape file to parquet
dises_shape_cleaned_consolidated.to_parquet(output_path_dises_shape, compression="zstd")