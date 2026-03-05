import geopandas as gpd
import pandas as pd
from shapely.ops import nearest_points
from shapely.geometry import box
import os
import numpy as np

# Input and output folders for CSB and road
input_folder_CSB = "data/edited/CSB/"
output_folder_CSB = "data/edited/CSB/"

states = snakemake.params.states

for state in states:
    
    input_path_CSB = os.path.join(input_folder_CSB, f"{state}_CSB1724_clipped.gpkg")
    output_path_geojson = os.path.join(output_folder_CSB, f"{state}_CSB1724_supplement_5_spatial.geojson")
    output_path_csv = os.path.join(output_folder_CSB, f"{state}_CSB1724_supplement_5_table.csv")

    # Load CSB_clipped datasets
    CSB_shape = gpd.read_file(input_path_CSB)
    
    # Reproject geometry to the same CRS (NAD83/CONUS Albers)
    CSB_shape = CSB_shape.to_crs(epsg=5070)
    
    # Create csb_supplement_5 dataset 
    CSB_supplement_5 = CSB_shape.copy()
    # Keep only the columns necessary for the spatial join
    cols_to_keep = ['CSBID', 'geometry']
    CSB_supplement_5 = CSB_supplement_5[cols_to_keep]
    
    
    print(f"Adding activities on neighboring fields for {state}...")
    
    # Keep land management activities from CSB_clipped
    columns_to_keep = [c for c in CSB_shape.columns 
                       if (c == "CSBID") or (c == "CSBACRES") or (c == "geometry") or c.startswith("CDL")]
    CSB_shape = CSB_shape[columns_to_keep]
    
    # Identify neighboring fields
    neighbor_gdf = gpd.sjoin(CSB_shape, CSB_shape, how="left", predicate="intersects", lsuffix="self", rsuffix="nbr")
    
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
    
    #---# Save geojson and csv files
    #CSB_supplement_5.to_file(output_path_geojson, driver="GeoJSON")
    attribute_table = CSB_supplement_5.drop(columns='geometry')
    attribute_table.to_csv(output_path_csv, index = False)
    
    