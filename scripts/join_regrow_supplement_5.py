#---# This code adds information about land management activities on neighboring fields
import geopandas as gpd
import pandas as pd
from shapely.ops import nearest_points
from shapely.geometry import box
import os
import numpy as np

# Input and output folders for Regrow and road
input_folder_Regrow = "data/edited/Regrow/"
output_folder_Regrow = "data/edited/Regrow/"

states = snakemake.params.states

for state in states:
    
    input_path_Regrow = os.path.join(input_folder_Regrow, f"{state}_regrow_shape_table.geojson")
    output_path_geojson = os.path.join(output_folder_Regrow, f"{state}_regrow_supplement_5_spatial.geojson")
    output_path_csv = os.path.join(output_folder_Regrow, f"{state}_regrow_supplement_5_table.csv")

    # Load regrow joined datasets
    regrow_shape = gpd.read_file(input_path_Regrow)
    
    # Reproject geometry to the same CRS (NAD83/CONUS Albers)
    regrow_shape = regrow_shape.to_crs(epsg=5070)
    
    # Create regrow_supplement_5 dataset 
    regrow_supplement_5 = regrow_shape.copy()
    # Keep only the columns necessary for the spatial join
    cols_to_keep = ['field_id', 'geometry']
    regrow_supplement_5 = regrow_supplement_5[cols_to_keep]
    
    
    # Adding activities on neighboring fields
    print(f"Adding activities on neighboring fields for {state}...")
    
    # Keep land management activities from regrow_shape
    columns_to_keep = [c for c in regrow_shape.columns 
                       if (c =="field_id") or (c == "area_acre") or (c == "geometry") or c.startswith("PHtill_1") or c.startswith("PHtill_2") 
                       or c.startswith("PPtill_1") or c.startswith("PPtill_2") or c.startswith("cover_1") or c.startswith("cover_2") or c.startswith("crop_1") or c.startswith("crop_2")]
    regrow_shape = regrow_shape[columns_to_keep]

    # Identify neighboring fields
    neighbor_gdf = gpd.sjoin(regrow_shape, regrow_shape, how="left", predicate="intersects", lsuffix="self", rsuffix="nbr")
    
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

    #---# Save geojson and csv files
    #regrow_supplement_5.to_file(output_path_geojson, driver="GeoJSON")
    attribute_table = regrow_supplement_5.drop(columns='geometry')
    attribute_table.to_csv(output_path_csv, index = False)
    
    