#---# This code adds information about land management activities on neighboring fields
import geopandas as gpd
import pandas as pd
from shapely.ops import nearest_points
from shapely.geometry import box
import os
import numpy as np

# Import parameters from Snakemake
regrow_input_folder = snakemake.params.regrow_input_dir
regrow_output_folder = snakemake.params.regrow_output_dir
states = snakemake.params.states


for state in states:
    
    regrow_shape_input_path = os.path.join(regrow_input_folder, f"{state}_regrow_fieldID_geometry.parquet")
    regrow_table_input_path = os.path.join(regrow_input_folder, f"{state}_regrow_table.parquet")
    output_path_spatial = os.path.join(regrow_output_folder, f"{state}_regrow_neighbor_field_mgmt_spatial.parquet")
    output_path_table = os.path.join(regrow_output_folder, f"{state}_regrow_neighbor_field_mgmt_table.parquet")

    # Load regrow shape and table data
    regrow_shape = gpd.read_parquet(regrow_shape_input_path)
    regrow_table = pd.read_parquet(regrow_table_input_path)
    
    # Create regrow shape_table file by merging regrow_shape and regrow_table
    regrow_shape_table = regrow_shape.merge(regrow_table, on="field_id", how="left")
    
    # Create regrow_supplement_5 dataset 
    regrow_supplement_5 = regrow_shape_table.copy()
    # Keep only the columns necessary for the spatial join
    cols_to_keep = ['field_id', 'geometry']
    regrow_supplement_5 = regrow_supplement_5[cols_to_keep]
    
    
    # Adding activities on neighboring fields
    print(f"Adding activities on neighboring fields for {state}...")
    
    # Keep land management activities from regrow_shape_table
    exact_cols = ["field_id", "area_acre", "geometry"]
    prefix_cols = ("PHtill_1", "PHtill_2", "PPtill_1", "PPtill_2", "cover_1", "cover_2", "crop_1", "crop_2")
    columns_to_keep = [c for c in regrow_shape_table.columns if c in exact_cols or c.startswith(prefix_cols)]
    regrow_shape_table = regrow_shape_table[columns_to_keep]

    # Identify neighboring fields
    neighbor_gdf = gpd.sjoin(regrow_shape_table, regrow_shape_table, how="left", predicate="intersects", lsuffix="self", rsuffix="nbr")
    
    # Remove self-matches (each field intersects with itself)
    neighbor_gdf = neighbor_gdf[neighbor_gdf["field_id_self"] != neighbor_gdf["field_id_nbr"]]
    
    # Extract columns for specific land management activities (PHtill, cover crop and PPtill) in neighboring fields
    nbr_activity_dict = {}
    for col in neighbor_gdf.columns:
        col_parts = col.split("_")
        if col.startswith("PHtill") and col.endswith("_nbr"):
            year = col_parts[1]
            nbr_activity_dict.setdefault(f"PHtill_nbr_{year}", []).append(col)
        elif col.startswith("cover") and col.endswith("_nbr"):
            year = col_parts[1]
            nbr_activity_dict.setdefault(f"cover_nbr_{year}", []).append(col)
        elif col.startswith("PPtill") and col.endswith("_nbr"):
            year = col_parts[1]
            nbr_activity_dict.setdefault(f"PPtill_nbr_{year}", []).append(col)
    
    # Create a new geodf with field IDs of fields that have at least one neighboring field
    neighbor_temp = neighbor_gdf.groupby("field_id_self")
    neighbor_till_cover = neighbor_temp.size().reset_index()[["field_id_self"]]
    
    # For a given field, compute the mean values of tillage intesity and cover crop practiced on all neighboring fields
    for new_col, cols in nbr_activity_dict.items():
        # sum of values across all rows and columns
        total_sum = neighbor_temp[cols].sum().sum(axis=1)
        total_count = neighbor_temp[cols].count().sum(axis=1)
        mean_values = total_sum.div(total_count.replace(0, np.nan))
        neighbor_till_cover[new_col] = neighbor_till_cover["field_id_self"].map(mean_values)
    
    # For a given field, compute the weighted values of main crops practiced on all neighboring fields
    # Weights are field areas. Crops are corn, soybean, wheat and others
    crop_cols = [c for c in neighbor_gdf.columns if c.startswith("crop_") and c.endswith("_nbr")]
    
    crop_groups = { 
    "corn": [1],
    "soybean": [5],
    "wheat": [22, 23, 24]
    }
    
    # Precompute set of all group codes
    all_group_codes = set(sum(crop_groups.values(), []))

    # List to store results
    results = []

    # Loop only over crop columns
    for col in crop_cols:
        # Valid values mask
        mask_valid = neighbor_gdf[col].notna()
        
        # Denominator per field (sum of area over valid values)
        denom = neighbor_gdf.loc[mask_valid].groupby("field_id_self")["area_acre_nbr"].sum()
        
        # Initialize dataframe for this column
        col_df = pd.DataFrame({"field_id_self": denom.index})
        
        # Compute numerator per group
        for group_name, codes in crop_groups.items():
            mask_group = mask_valid & neighbor_gdf[col].isin(codes)
            numer = neighbor_gdf.loc[mask_group].groupby("field_id_self")["area_acre_nbr"].sum()
            numer = numer.reindex(denom.index, fill_value=0)  # align with denom
            col_df[f"crop_{group_name}_nbr_{col.split('_')[1]}_{col.split('_')[2]}"] = np.where(denom != 0, numer / denom, np.nan)

        # Compute "other"
        mask_other = mask_valid & ~neighbor_gdf[col].isin(all_group_codes)
        numer = neighbor_gdf.loc[mask_other].groupby("field_id_self")["area_acre_nbr"].sum()
        numer = numer.reindex(denom.index, fill_value=0)
        col_df[f"crop_other_nbr_{col.split('_')[1]}_{col.split('_')[2]}"] = np.where(denom != 0, numer / denom, np.nan)

        results.append(col_df)

    neighbor_crop = neighbor_temp.size().reset_index()[["field_id_self"]]

    for df in results[0:]:
        neighbor_crop = pd.merge(neighbor_crop, df, on="field_id_self", how="outer")
    
    # Join land management practices on neighboring fields back to the main dataset
    regrow_supplement_5 = regrow_supplement_5.merge(neighbor_till_cover, left_on="field_id", right_on="field_id_self", how='left')
    regrow_supplement_5 = regrow_supplement_5.merge(neighbor_crop, left_on="field_id", right_on="field_id_self", how='left')
    regrow_supplement_5 = regrow_supplement_5.drop(columns=['field_id_self_x', 'field_id_self_y'])
    
    print(f"Adding activities on neighboring fields for {state} is complete.")
    
    # Convert float64 columns to float32 to save memory
    float64_cols = regrow_supplement_5.select_dtypes(include=["float64"]).columns
    regrow_supplement_5[float64_cols] = regrow_supplement_5[float64_cols].astype("float32")

    #---# Save files w/ and w/o geometry
    #regrow_supplement_5.to_parquet(output_path_spatial, compression="zstd")
    attribute_table = regrow_supplement_5.drop(columns='geometry')
    attribute_table.to_parquet(output_path_table, index = False, compression="zstd")
    
    