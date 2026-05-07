import geopandas as gpd
import pandas as pd
from shapely.ops import nearest_points
from shapely.geometry import box
import os
import numpy as np

# Import parameters from Snakemake
CSB_input_folder = snakemake.params.CSB_input_dir
CSB_output_folder = snakemake.params.CSB_output_dir
states = snakemake.params.states
CSB_years = snakemake.params.CSB_years


for year in CSB_years:
    for state in states:
        
        CSB_shape_input_path = os.path.join(CSB_input_folder, f"{state}_CSB{year}_CSBID_geometry.parquet")
        CSB_table_input_path = os.path.join(CSB_input_folder, f"{state}_CSB{year}_table.parquet")
        output_path_spatial = os.path.join(CSB_output_folder, f"{state}_CSB{year}_supplement_5_spatial.parquet")
        output_path_table = os.path.join(CSB_output_folder, f"{state}_CSB{year}_supplement_5_table.parquet")
        
        # Load CSB data
        csb_shape = gpd.read_parquet(CSB_shape_input_path)
        csb_table = pd.read_parquet(CSB_table_input_path)
        
        # Create CSB shape_table file by merging csb_shape and csb_table
        csb_shape_table = csb_shape.merge(csb_table, on="CSBID", how="left")
        
        # Create csb_supplement_5 dataset 
        CSB_supplement_5 = csb_shape_table.copy()
        # Keep only the columns necessary for the spatial join
        cols_to_keep = ['CSBID', 'geometry']
        CSB_supplement_5 = CSB_supplement_5[cols_to_keep]
        
        
        print(f"Adding activities on neighboring fields for {state}...")
        
        # Keep land management activities from CSB_clipped
        exact_cols = ["CSBID", "CSBACRES", "geometry"]
        prefix_cols = ("CDL")
        columns_to_keep = [c for c in csb_shape_table.columns if c in exact_cols or c.startswith(prefix_cols)]
        csb_shape_table = csb_shape_table[columns_to_keep]
        
        # Identify neighboring fields
        neighbor_gdf = gpd.sjoin(csb_shape_table, csb_shape_table, how="left", predicate="intersects", lsuffix="self", rsuffix="nbr")
        
        # Remove self-matches (each field intersects with itself)
        neighbor_gdf = neighbor_gdf[neighbor_gdf["CSBID_self"] != neighbor_gdf["CSBID_nbr"]]
        
        # Group spatially merged df by CSBID
        neighbor_temp = neighbor_gdf.groupby("CSBID_self")
            
        # For a given field, compute the weighted values of main crops practiced on all neighboring fields
        # Weights are field areas. Crops are corn, soybean, wheat and others
        crop_cols = [c for c in neighbor_gdf.columns if c.startswith("CDL") and c.endswith("_nbr")]
        
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
            denom = neighbor_gdf.loc[mask_valid].groupby("CSBID_self")["CSBACRES_nbr"].sum()
            
            # Initialize dataframe for this column
            col_df = pd.DataFrame({"CSBID_self": denom.index})
            
            # Compute numerator per group
            for group_name, codes in crop_groups.items():
                mask_group = mask_valid & neighbor_gdf[col].isin(codes)
                numer = neighbor_gdf.loc[mask_group].groupby("CSBID_self")["CSBACRES_nbr"].sum()
                numer = numer.reindex(denom.index, fill_value=0)  # align with denom
                col_df[f"CDL_{group_name}_nbr_{col[3:7]}"] = np.where(denom != 0, numer / denom, np.nan)

            # Compute "other"
            mask_other = mask_valid & ~neighbor_gdf[col].isin(all_group_codes)
            numer = neighbor_gdf.loc[mask_other].groupby("CSBID_self")["CSBACRES_nbr"].sum()
            numer = numer.reindex(denom.index, fill_value=0)
            col_df[f"CDL_other_nbr_{col[3:7]}"] = np.where(denom != 0, numer / denom, np.nan)

            results.append(col_df)

        neighbor_crop = neighbor_temp.size().reset_index()[["CSBID_self"]]

        for df in results[0:]:
            neighbor_crop = pd.merge(neighbor_crop, df, on="CSBID_self", how="outer")
        
        # Join land management practices on neighboring fields back to the main dataset
        CSB_supplement_5 = CSB_supplement_5.merge(neighbor_crop, left_on="CSBID", right_on="CSBID_self", how='left')
        CSB_supplement_5 = CSB_supplement_5.drop(columns=['CSBID_self'])
        
        print(f"Adding activities on neighboring fields for {state} is complete.")
        
        # Convert float64 columns to float32 to save memory
        float64_cols = CSB_supplement_5.select_dtypes(include=["float64"]).columns
        CSB_supplement_5[float64_cols] = CSB_supplement_5[float64_cols].astype("float32")
        
        #---# Save files w/ and w/o geometry
        #CSB_supplement_5.to_parquet(output_path_spatial, compression="zstd")
        attribute_table = CSB_supplement_5.drop(columns='geometry')
        attribute_table.to_parquet(output_path_table, index = False, compression="zstd")
    
    