import geopandas as gpd
import pandas as pd
import numpy as np
import os
from pathlib import Path
from sklearn.neighbors import KDTree
import shutil


# Input and output folders for CSB
input_folder_CSB = "data/edited/CSB/"
output_folder_CSB = "data/edited/CSB/"

states = snakemake.params.states
crops = snakemake.params.crops
number_of_neighbors = snakemake.params.number_of_neighbors

# Create a temporary folder to run this code
temp_dir = Path("data/edited/Grain Price/temp/")
# Delete temporary folder if it already exists
if temp_dir.exists():
    shutil.rmtree(temp_dir)
# Create parent folder if it doesn't exist
temp_dir.mkdir(parents=True, exist_ok=False)

# Create a list of all months from 2010 to 2025
dates = pd.date_range(start="2010-01-01", end="2025-12-01", freq="MS")

for state_CSB in states:
    
    # Path to CSB file
    input_path_CSB = os.path.join(input_folder_CSB, f"{state_CSB}_CSB1724_supplement_1_spatial.geojson")
    # Load CSB joined datasets
    CSB_supplement_1 = gpd.read_file(input_path_CSB)
    # Replace St to Saint in the county_name column
    CSB_supplement_1['county_name'] = CSB_supplement_1['county_name'].str.replace(r'\bSt\.?\b', 'Saint', regex=True)
    # Select columns for the analysis
    selected_columns = ['CSBID', 'county_name', 'geometry']
    CSB_supplement_1 = CSB_supplement_1[selected_columns]


    for crop in crops:
        # Create empty dataset with Date (month) column
        df_crop_elevator = pd.DataFrame({"Date": dates})
        df_crop_county = pd.DataFrame({"Date": dates})
        
        #---# Elevator-level crop price
        for state in states:
            
            df_elevator_price = pd.read_excel(f"data/edited/Grain Price/{state}_{crop}_elevator_level.xlsx", sheet_name="Values", header=[0, 1])
            df_elevator_price[('History', 'Date')] = pd.to_datetime(df_elevator_price[('History', 'Date')])

            new_cols = {}
            for col in df_elevator_price.columns.get_level_values(0)[1:]:
                new_cols[(col, 'Average')] = (df_elevator_price[(col, 'Low')] + df_elevator_price[(col, 'High')]) / 2
            avg_df_elevator_price = pd.DataFrame(new_cols, index=df_elevator_price.index)
            df_elevator_price = pd.concat([df_elevator_price, avg_df_elevator_price], axis=1)

            df_elevator_price.columns = [f"{ticker}_{price}" for ticker, price in df_elevator_price.columns]

            df_elevator_price.rename(columns={"History_Date": "Date"}, inplace = True)
            
            df_crop_elevator = df_crop_elevator.merge(df_elevator_price, on = "Date", how = "left")

        # Drop duplicate column names (keep the first occurence)
        df_crop_elevator = df_crop_elevator.loc[:, ~df_crop_elevator.columns.duplicated()]
        # Select only Date and columns with average elevator prices
        df_crop_avg = df_crop_elevator.loc[:, df_crop_elevator.columns.str.endswith(('Date', 'Average'))]
        # Remove "Average" from column names
        df_crop_avg.columns = df_crop_avg.columns.str.replace('_Average', '', regex=False)
        # Select values from January 2014 onwards
        df_crop_avg = df_crop_avg[df_crop_avg['Date'] >= '2014-01-01']
        df_crop_avg.reset_index(drop = True, inplace = True)
        # Drop the elevators with few observations: fewer than 12 observations between 2015 and 2025
        mask = df_crop_avg['Date'].between('2015-01-01', '2024-12-31')
        cols_to_keep = (df_crop_avg.loc[mask].dropna(axis = 1, thresh = 12).columns)
        df_crop_avg = df_crop_avg.loc[:, cols_to_keep]

        # Set first column as index
        df_crop_avg = df_crop_avg.set_index(df_crop_avg.columns[0])
        # Transpose
        df_crop_avg = df_crop_avg.T
        # Rename columns using YYYYMM format
        df_crop_avg.columns = [
            f"{crop}_price_elevator_{col.strftime('%Y%m')}"
            for col in df_crop_avg.columns]
        # Reset index so original column names become first column
        df_crop_avg = df_crop_avg.reset_index().rename(columns={"index": "ticker"})
        
        elevator_location = gpd.read_file(f"data/edited/Grain price/{crop}_elevator_location.geojson")
        elevator_location = elevator_location[~elevator_location.geometry.isna()]
        

        # Compute price of the nearest elevator 
        gdf_parcel = CSB_supplement_1.to_crs(epsg=5070)
        elevator_location = elevator_location.to_crs(epsg=5070)

        # Compute coordinates of parcel centroids and coondinates of elevators
        centroids = gdf_parcel.geometry.centroid
        parcel_coords = np.column_stack((centroids.x, centroids.y))
        elevator_coords = np.array([[geom.x, geom.y] for geom in elevator_location.geometry])

        # Create a KDTree to speed up searching for the nearest elevators
        tree = KDTree(elevator_coords)
        
        K = 3 # number of nearest elevators considered in the nearest elevator search
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
                index=gdf_parcel.index,
                dtype=np.float32
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
        gdf_parcel_price.insert(0, 'CSBID', gdf_parcel['CSBID'])

        # Rename columns
        gdf_parcel_price.columns = [col.replace(f'{crop}_price_elevator_', f'{crop}_price_elevator_nearest_') if col.startswith(f'{crop}_price_elevator_') else col for col in gdf_parcel_price.columns]
        
        # Save price file
        temp_file_path = os.path.join(temp_dir, f"{state_CSB}_{crop}_price_nearest.csv")
        gdf_parcel_price.to_csv(temp_file_path, index=False)
        
        # Remove heavy variables from memory
        del gdf_parcel_price
        
        
        # Compute weighted average price of N-nearest elevator 
        # Stack distances to all N-nearest elevators together
        data_stack = np.stack([df.values for df in all_nearest_elevators], axis=1)
        col_names = all_nearest_elevators[0].columns
        
        # Remove heavy variables from memory
        del all_nearest_elevators
        
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
        normalized_weights = np.divide(weights_expanded, weight_sums, where=weight_sums!=0)
        
        weighted_mean = np.nansum(data_stack * normalized_weights, axis=1)
        
        weighted_prices = pd.DataFrame(weighted_mean, columns=col_names)
        weighted_prices.replace(0, np.nan, inplace=True)
        weighted_prices.insert(0, 'CSBID', gdf_parcel['CSBID'])
        weighted_prices.columns = [col.replace(f'{crop}_price_elevator_', f'{crop}_price_elevator_{N}-nearest_') if col.startswith(f'{crop}_price_elevator_') else col for col in weighted_prices.columns]
                
        # Save price file
        temp_file_path = os.path.join(temp_dir, f"{state_CSB}_{crop}_price_{N}-nearest.csv")
        weighted_prices.to_csv(temp_file_path, index=False)
        
        # Remove heavy variables from memory
        del data_stack, weights_expanded, normalized_weights, weighted_mean, weighted_prices

        print(f'Elevator-level {crop} price in {state_CSB} is computed and saved.')
    

        #---# County-level crop price
        if crop != 'wheat':
            for state in states:
                
                df_county_price = pd.read_excel(f"data/edited/Grain Price/{state}_{crop}_county_level.xlsx", sheet_name="Values", header=[0, 1])
                df_county_price[('History', 'Date')] = pd.to_datetime(df_county_price[('History', 'Date')])

                new_cols = {}
                for col in df_county_price.columns.get_level_values(0)[1:]:
                    new_cols[(col, 'Average')] = (df_county_price[(col, 'Low')] + df_county_price[(col, 'High')]) / 2
                avg_df_county_price = pd.DataFrame(new_cols, index=df_county_price.index)
                df_county_price = pd.concat([df_county_price, avg_df_county_price], axis=1)

                df_county_price.columns = [f"{ticker}_{price}" for ticker, price in df_county_price.columns]

                df_county_price.rename(columns={"History_Date": "Date"}, inplace = True)
                
                df_crop_county = df_crop_county.merge(df_county_price, on = "Date", how = "left")

            # Drop duplicate column names (keep the first occurence)
            df_crop_county = df_crop_county.loc[:, ~df_crop_county.columns.duplicated()]
            # Select only Date and columns with average elevator prices
            df_crop_county_avg = df_crop_county.loc[:, df_crop_county.columns.str.endswith(('Date', 'Average'))]
            # Remove "Average" from column names
            df_crop_county_avg.columns = df_crop_county_avg.columns.str.replace('_Average', '', regex=False)
            # Select values from January 2014 onwards
            df_crop_county_avg = df_crop_county_avg[df_crop_county_avg['Date'] >= '2014-01-01']
            df_crop_county_avg.reset_index(drop = True, inplace = True)
            # Drop the elevators with few observations: fewer than 12 observations between 2015 and 2025
            mask = df_crop_county_avg['Date'].between('2015-01-01', '2024-12-31')
            cols_to_keep = (df_crop_county_avg.loc[mask].dropna(axis = 1, thresh = 12).columns)
            df_crop_county_avg = df_crop_county_avg.loc[:, cols_to_keep]
            
            # Set first column as index
            df_crop_county_avg = df_crop_county_avg.set_index(df_crop_county_avg.columns[0])
            # Transpose
            df_crop_county_avg = df_crop_county_avg.T
            # Rename columns using YYYYMM format
            df_crop_county_avg.columns = [
                f"{crop}_price_county_{col.strftime('%Y%m')}"
                for col in df_crop_county_avg.columns]
            # Reset index so original column names become first column
            df_crop_county_avg = df_crop_county_avg.reset_index().rename(columns={"index": "ticker"})
            
            index_county_location = gpd.read_file(f"data/edited/Grain price/{crop}_index_county_location.geojson")
            index_county_location = index_county_location[~index_county_location.geometry.isna()]
            
        
            gdf_parcel = CSB_supplement_1.to_crs(epsg=5070)
            index_county_location = index_county_location.to_crs(epsg=5070)
            
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
                    index=gdf_parcel.index,
                    dtype=np.float32
                )

                all_nearest_counties.append(nearest_county_n)
                
            # Start with first nearest
            county_price_index = all_nearest_counties[0].copy()

            # Fill missing values with second and then third nearest
            for i in range(N):
                vals = all_nearest_counties[i]
                county_price_index = county_price_index.fillna(vals)

            # Add parcel_id back
            county_price_index.insert(0, 'CSBID', gdf_parcel['CSBID'])
            
            # Save price file
            temp_file_path = os.path.join(temp_dir, f"{state_CSB}_{crop}_price_county_index.csv")
            county_price_index.to_csv(temp_file_path, index=False)
            
            # Remove heavy variables from memory 
            del value_dict, all_nearest_counties

            print(f'County-level {crop} price in {state_CSB} is computed and saved.')
                

# Merge the prices for all crops for each state and save the merged files
for state_CSB in states:
    all_crop_prices = pd.DataFrame()
    # Collect filenames of price data for a given state
    csv_files = sorted(temp_dir.glob(f"{state_CSB}_*.csv"))
    
    for csv_file in csv_files:
        df = pd.read_csv(csv_file)
        # Combine all matching CSVs into one DataFrame
        if all_crop_prices.empty:
            all_crop_prices = df.copy()
        else:
            all_crop_prices = all_crop_prices.merge(df, on = 'CSBID', how = 'outer')
        csv_file.unlink()

    # Save the merged price dataset for a given state
    all_crop_prices.to_csv(f'data/edited/CSB/{state_CSB}_CSB1724_supplement_7_table.csv', index=False)
    print(f"Saved for {state_CSB}")

# Delete temporary folder
shutil.rmtree(temp_dir)