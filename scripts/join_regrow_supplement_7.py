import geopandas as gpd
import pandas as pd
import numpy as np
import os
from pathlib import Path
from sklearn.neighbors import KDTree
import shutil

# Import parameters from Snakemake
regrow_input_folder = snakemake.params.regrow_input_dir
regrow_output_folder = snakemake.params.regrow_output_dir
crop_price_input_folder = snakemake.params.crop_price_input_dir
states = snakemake.params.states
crops = snakemake.params.crops
number_of_neighbors = snakemake.params.number_of_neighbors
K = snakemake.params.K
target_CRS = snakemake.params.target_CRS


def combine_clean_regrow_files(regrow_input_folder, state):
    
    # Path to Regrow files
    regrow_shape_input_path = os.path.join(regrow_input_folder, f"{state}_regrow_fieldID_geometry.parquet")
    regrow_supplement_1_table_input_path = os.path.join(regrow_input_folder, f"{state}_regrow_supplement_1_table.parquet")
    
    # Load Regrow datasets
    regrow_shape = gpd.read_parquet(regrow_shape_input_path)
    regrow_supplement_1_table = pd.read_parquet(regrow_supplement_1_table_input_path)
    
    # Create regrow regrow_supplement_1 file by merging regrow_shape and regrow_supplement_1
    regrow_supplement_1 = regrow_shape.merge(regrow_supplement_1_table, on="field_id", how="left")
    
    # Replace St to Saint in the county_name column
    regrow_supplement_1['county_name'] = regrow_supplement_1['county_name'].str.replace(r'\bSt\.?\b', 'Saint', regex=True)
    # Select columns for the analysis
    selected_columns = ['field_id', 'county_name', 'geometry']
    regrow_supplement_1 = regrow_supplement_1[selected_columns]
    
    return regrow_supplement_1
    
    
def aggregate_crop_price_elevator_level(regrow_supplement_1, crop, state, number_of_neighbors, K, target_CRS, crop_price_input_folder, temp_dir):
    
    #---# Elevator-level crop price
    # Read monthly average elevator price
    df_crop_avg = pd.read_csv(Path(crop_price_input_folder) / f"{crop}_monthly_average_elevator_price.csv")
    
    # Read elevator location
    elevator_location = gpd.read_file(Path(crop_price_input_folder) / f"{crop}_elevator_location.geojson")
    elevator_location = elevator_location[~elevator_location.geometry.isna()]

    
    # Compute price of the nearest elevator 
    gdf_parcel = regrow_supplement_1.to_crs(target_CRS)
    elevator_location = elevator_location.to_crs(target_CRS)

    # Compute coordinates of parcel centroids and coondinates of elevators
    centroids = gdf_parcel.geometry.centroid
    parcel_coords = np.column_stack((centroids.x, centroids.y))
    elevator_coords = np.array([[geom.x, geom.y] for geom in elevator_location.geometry])

    # Create a KDTree to speed up searching for the nearest elevators
    tree = KDTree(elevator_coords)
    
    K = K # number of nearest elevators considered in the nearest elevator search
    N = number_of_neighbors  # number of nearest elevators considered in the N-nearest elevator average

    distances, indices = tree.query(parcel_coords, k = N)
    
    # Keep only month columns
    month_cols = df_crop_avg.columns[1:]

    # Dictionary: ticker -> DataFrame of months (1 row per ticker)
    value_dict = {row['ticker']: row[month_cols] for _, row in df_crop_avg.iterrows()}
    
    all_nearest_elevators = []

    for n in range(N):
        # get nth nearest elevator index for each parcel
        obj_indices_n = indices[:, n]
        tickers_n = elevator_location.iloc[obj_indices_n]['ticker'].values

        # get values for these tickers
        nearest_elevator_n = pd.DataFrame(
            [value_dict[t] for t in tickers_n],
            columns=month_cols,
            index=gdf_parcel.index
        )

        all_nearest_elevators.append(nearest_elevator_n)
        
    # Start with first nearest
    gdf_parcel_price = all_nearest_elevators[0].copy()

    # Fill missing values with second and then third nearest
    for k in range(1, K):
        vals = all_nearest_elevators[k]
        gdf_parcel_price = gdf_parcel_price.fillna(vals)

    # Remove heavy variables from memory
    del value_dict
    
    # Add parcel_id back
    gdf_parcel_price.insert(0, 'field_id', gdf_parcel['field_id'])

    # Rename columns
    gdf_parcel_price.columns = [col.replace(f'{crop}_price_elevator_', f'{crop}_price_elevator_nearest_') if col.startswith(f'{crop}_price_elevator_') else col for col in gdf_parcel_price.columns]
    
    # Save price file
    temp_file_path = os.path.join(temp_dir, f"{state}_{crop}_price_nearest.parquet")
    gdf_parcel_price.to_parquet(temp_file_path , index=False, compression="zstd")
    
    
    # Compute weighted average price of N-nearest elevator 
    # Stack distances to all N-nearest elevators together
    data_stack = np.stack([df.values for df in all_nearest_elevators], axis=1)
    
    # Convert distances to weights (closer = bigger weight)
    weights = 1 / distances
    
    # Expand weights to match shape (num_parcels x N x num_months)
    weights_expanded = np.expand_dims(weights, axis=2)
    weights_expanded = np.repeat(weights_expanded, data_stack.shape[2], axis=2)
    
    # Mask NaNs
    mask = ~np.isnan(data_stack)
    # Apply mask to weights (zero weight for NaNs)
    weights_expanded = weights_expanded * mask
    del mask
    
    # Normalize weights along N axis so sum = 1 for available values
    weight_sums = weights_expanded.sum(axis=1, keepdims=True)
    normalized_weights = np.divide(weights_expanded, weight_sums, where=weight_sums != 0, out=np.zeros_like(weights_expanded, dtype=float))
    
    weighted_mean = np.nansum(data_stack * normalized_weights, axis=1)
    
    weighted_prices = pd.DataFrame(weighted_mean, columns=all_nearest_elevators[0].columns)
    weighted_prices.replace(0, np.nan, inplace=True)
    weighted_prices.insert(0, 'field_id', gdf_parcel['field_id'])
    weighted_prices.columns = [col.replace(f'{crop}_price_elevator_', f'{crop}_price_elevator_{N}-nearest_') if col.startswith(f'{crop}_price_elevator_') else col for col in weighted_prices.columns]
    
    # Save price file
    temp_file_path = os.path.join(temp_dir, f"{state}_{crop}_price_{N}-nearest.parquet")
    weighted_prices.to_parquet(temp_file_path, index=False, compression="zstd")

    print(f'Elevator-level {crop} price in {state} is computed and saved.')


def aggregate_crop_price_county_level(regrow_supplement_1, crop, state, number_of_neighbors, target_CRS, crop_price_input_folder, temp_dir):

    if crop != 'wheat':
        # Read monthly average county price
        df_crop_county_avg = pd.read_csv(Path(crop_price_input_folder) / f"{crop}_monthly_average_county_price.csv")
        
        # Read county price index location
        index_county_location = gpd.read_file(Path(crop_price_input_folder) / f"{crop}_index_county_location.geojson")
        index_county_location = index_county_location[~index_county_location.geometry.isna()]

        # Convert geometries of both datasets to target CRS
        gdf_parcel = regrow_supplement_1.to_crs(target_CRS)
        index_county_location = index_county_location.to_crs(target_CRS)
        
        # Compute coordinates of parcel centroids and coondinates of county centroids
        centroids = gdf_parcel.geometry.centroid
        parcel_coords = np.column_stack((centroids.x, centroids.y))
        county_coords = np.array([[geom.x, geom.y] for geom in index_county_location.geometry])

        # Create a KDTree to speed up searching for the nearest elevators
        tree = KDTree(county_coords)
        
        # Number of nearest counties to average
        N = number_of_neighbors

        distances, indices = tree.query(parcel_coords, k = N)
        
        county_name_to_index = dict(zip(index_county_location['county'], index_county_location.index))

        final_indices = []
        for i, row in gdf_parcel.iterrows():
            nearest_idx = indices[i].copy()
            
            parcel_county = row['county_name']
            own_county_idx = county_name_to_index.get(parcel_county)
            
            if own_county_idx is not None:
                # Remove own county if it is already in nearest list
                nearest_idx = nearest_idx[nearest_idx != own_county_idx]
                
                # Put own county first
                nearest_idx = np.insert(nearest_idx, 0, own_county_idx)
                
                # Keep only N counties in a list of nearest counties
                nearest_idx = nearest_idx[:N]
            
            final_indices.append(nearest_idx)

        final_indices = np.vstack(final_indices)
        
        # Keep only month columns
        month_cols = df_crop_county_avg.columns[1:]

        # Dictionary: ticker -> DataFrame of months (1 row per ticker)
        value_dict = {row['ticker']: row[month_cols] for _, row in df_crop_county_avg.iterrows()}
        
        all_nearest_counties = []

        for n in range(N):
            # get nth nearest county index for each parcel
            obj_indices_n = final_indices[:, n]
            tickers_n = index_county_location.iloc[obj_indices_n]['ticker'].values

            # get values for these tickers
            nearest_county_n = pd.DataFrame(
                [value_dict[t] for t in tickers_n],
                columns=month_cols,
                index=gdf_parcel.index
            )

            all_nearest_counties.append(nearest_county_n)
            
        # Start with first nearest
        county_price_index = all_nearest_counties[0].copy()

        # Fill missing values with second and then third nearest
        for i in range(N):
            vals = all_nearest_counties[i]
            county_price_index = county_price_index.fillna(vals)

        # Add parcel_id back
        county_price_index.insert(0, 'field_id', gdf_parcel['field_id'])
        
        # Save price file
        temp_file_path = os.path.join(temp_dir, f"{state}_{crop}_price_county_index.parquet")
        county_price_index.to_parquet(temp_file_path, index=False, compression="zstd")
        
        # Remove heavy variables from memory 
        del value_dict, all_nearest_counties

        print(f'County-level {crop} price in {state} is computed and saved.')
  
              
# Merge the prices for all crops for each state and save the merged files
def merge_all_crop_prices_by_state(state, regrow_output_folder, temp_dir):
    all_crop_prices = pd.DataFrame()
    # Collect filenames of price data for a given state
    temp_files = sorted(temp_dir.glob(f"{state}_*.parquet"))
    
    for temp_file in temp_files:
        df = pd.read_parquet(temp_file)
        # Combine all matching files into one DataFrame
        if all_crop_prices.empty:
            all_crop_prices = df.copy()
        else:
            all_crop_prices = all_crop_prices.merge(df, on = 'field_id', how = 'outer')
        temp_file.unlink()

    # Convert float64 columns to float32 to save memory
    float64_cols = all_crop_prices.select_dtypes(include=["float64"]).columns
    all_crop_prices[float64_cols] = all_crop_prices[float64_cols].astype("float32")
    
    # Save the merged price dataset for a given state
    output_path = os.path.join(regrow_output_folder, f"{state}_regrow_supplement_7_table.parquet")
    all_crop_prices.to_parquet(output_path, index=False, compression="zstd")
    print(f"Saved for {state}")


# Main script
# Create a temporary folder to run this code
temp_dir =  Path(crop_price_input_folder) / "temp"
# Delete temporary folder if it already exists
if temp_dir.exists():
    shutil.rmtree(temp_dir)
# Create parent folder if it does not exist
temp_dir.mkdir(parents=True, exist_ok=False)

for state in states:
    regrow_supplement_1 = combine_clean_regrow_files(regrow_input_folder, state)
     
    for crop in crops: 
         aggregate_crop_price_elevator_level(regrow_supplement_1, crop, state, number_of_neighbors, K, target_CRS, crop_price_input_folder, temp_dir)
         aggregate_crop_price_county_level(regrow_supplement_1, crop, state, number_of_neighbors, target_CRS, crop_price_input_folder, temp_dir)
    
    merge_all_crop_prices_by_state(state, regrow_output_folder, temp_dir)
         
# Delete temporary folder
shutil.rmtree(temp_dir)