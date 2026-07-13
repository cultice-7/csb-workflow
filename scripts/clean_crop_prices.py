import geopandas as gpd
import pandas as pd
import numpy as np
from geopy.geocoders import Nominatim
from shapely.geometry import Point
import time
from pathlib import Path

# Import variables from snakemake
states = snakemake.params.states
crops = snakemake.params.crops
price_levels = snakemake.params.price_levels
drop_threshold = snakemake.params.drop_threshold
input_dir = Path(snakemake.params.input_dir)
output_dir = Path(snakemake.params.output_dir)


def merge_crop_price_by_state(crop, states, level, input_dir):
    # Create a list of all months from 2010 to 2025
    dates = pd.date_range(start="2010-01-01", end="2025-12-01", freq="MS")
    # Create empty dataset with Date (month) column
    df_crop_price = pd.DataFrame({"Date": dates})
    
    # Merge price data across all available states
    for state in states:
        
        df_state_crop_price = pd.read_excel(input_dir / f"{state}_{crop}_{level}_level.xlsx", sheet_name="Values", header=[0, 1])
        df_state_crop_price[('History', 'Date')] = pd.to_datetime(df_state_crop_price[('History', 'Date')])

        new_cols = {}
        for col in df_state_crop_price.columns.get_level_values(0)[1:]:
            new_cols[(col, 'Average')] = (df_state_crop_price[(col, 'Low')] + df_state_crop_price[(col, 'High')]) / 2
        avg_df_state_crop_price = pd.DataFrame(new_cols, index=df_state_crop_price.index)
        df_state_crop_price = pd.concat([df_state_crop_price, avg_df_state_crop_price], axis=1)

        df_state_crop_price.columns = [f"{ticker}_{price}" for ticker, price in df_state_crop_price.columns]

        df_state_crop_price.rename(columns={"History_Date": "Date"}, inplace = True)
        
        df_crop_price = df_crop_price.merge(df_state_crop_price, on = "Date", how = "left")

    return df_crop_price


def clean_transform_save_merged_price_data(df_crop_price, crop, level, drop_threshold, output_dir):
    # Drop duplicate column names (keep the first occurence)
    df_crop_price = df_crop_price.loc[:, ~df_crop_price.columns.duplicated()]
    # Select only Date and columns with average elevator prices
    df_crop_price_avg = df_crop_price.loc[:, df_crop_price.columns.str.endswith(('Date', 'Average'))]
    # Remove "Average" from column names
    df_crop_price_avg.columns = df_crop_price_avg.columns.str.replace('_Average', '', regex=False)
    # Select values from January 2014 onwards
    df_crop_price_avg = df_crop_price_avg[df_crop_price_avg['Date'] >= '2014-01-01']
    df_crop_price_avg.reset_index(drop = True, inplace = True)
    # Drop the elevators with few observations: fewer than "drop_threshold" observations between 2015 and 2025
    mask = df_crop_price_avg['Date'].between('2015-01-01', '2025-12-31')
    cols_to_keep = (df_crop_price_avg.loc[mask].dropna(axis = 1, thresh = drop_threshold).columns)
    df_crop_price_avg = df_crop_price_avg.loc[:, cols_to_keep]

    # Set first column as index
    df_crop_price_avg = df_crop_price_avg.set_index(df_crop_price_avg.columns[0])
    # Transpose
    df_crop_price_avg = df_crop_price_avg.T
    # Rename columns using YYYYMM format
    df_crop_price_avg.columns = [
        f"{crop}_price_{level}_{col.strftime('%Y%m')}"
        for col in df_crop_price_avg.columns]
    # Reset index so original column names become first column
    df_crop_price_avg = df_crop_price_avg.reset_index().rename(columns={"index": "ticker"})
    
    # Save monthly average price data
    df_crop_price_avg.to_csv(output_dir / f"{crop}_monthly_average_{level}_price.csv", index=False)
    print(f"{crop} {level} monthly average price is saved")
    
    return df_crop_price_avg



def geocode_location_elevator(location):
    geolocator = Nominatim(user_agent="geo_centroid")
    try:
        time.sleep(1)
        loc = geolocator.geocode(location, addressdetails = True)
        point = Point(loc.longitude, loc.latitude)
        county = loc.raw['address'].get('county')
        state = loc.raw['address'].get('state')
        return point, county, state
    except Exception:
        try:
            time.sleep(1)
            loc = geolocator.geocode(location, addressdetails = True)
            point = Point(loc.longitude, loc.latitude)
            county = loc.raw['address'].get('county')
            state = loc.raw['address'].get('state')
            return point, county, state
        except Exception:
            return None, None, None


def geocode_location_county(location):
    geolocator = Nominatim(user_agent="geo_centroid")
    try: 
        time.sleep(1) 
        loc = geolocator.geocode(location)
        point = Point(loc.longitude, loc.latitude)
        return point
    except Exception:
        try:
            time.sleep(1) 
            loc = geolocator.geocode(location)
            point = Point(loc.longitude, loc.latitude)
            return point
        except Exception:
            return None


def determine_save_elevator_location(df_crop_price, df_crop_price_avg, crop, input_dir, output_dir):
    names = []
    tickers = []
    for col in df_crop_price.columns:
        if col.endswith("_Name"):
            ticker = col.split("_")[0] 
            value = df_crop_price[col].dropna().unique()[0] if df_crop_price[col].notna().any() else None
            tickers.append(ticker)
            names.append(value)

    # Create ticker + elevator location dataframe
    elevator_location = pd.DataFrame({
        "ticker": tickers,
        "elevator description": names
    })
    elevator_location = elevator_location[elevator_location['ticker'].isin(df_crop_price_avg['ticker'])]
    elevator_location.reset_index(drop = True, inplace = True)

    # Extract elevator location from elevator description
    elevator_location['elevator location'] = elevator_location['elevator description'].str.split(";").str[1].str.strip()
    # Replace "St" and "St." with "Saint" in elevator location
    elevator_location['elevator location'] = elevator_location['elevator location'].str.replace(r'\bSt\.?\b', 'Saint', regex=True)

    # Replace incorrect elevator locations
    elevator_location_corrected = pd.read_excel(input_dir / f"corrected_elevator_location_corn_soybeans_wheat.xlsx")
    elevator_location['elevator location'] = elevator_location['ticker'].map(elevator_location_corrected.set_index('ticker')['corrected location']).fillna(elevator_location['elevator location'])
    
    # Apply geocode_location to find the centroids of towns where elevators are located
    elevator_location[['geometry', 'county', 'state']] = elevator_location['elevator location'].apply(lambda x: pd.Series(geocode_location_elevator(x)))
    # Extract crop type from elevator description
    elevator_location['crop type'] = elevator_location['elevator description'].str.split(";").str[2].str.split().str[0].str.lower()
    # Create a GeoDataFrame with elevator locations
    elevator_location = gpd.GeoDataFrame(elevator_location, geometry = "geometry", crs = "EPSG:4326")
    
    # Verify that the crop type listed for elevators match the crop type we are currently working with
    crop_type_check = (elevator_location['crop type'].notna()) & (elevator_location['crop type'] != crop)
    if crop_type_check.any():
        print("Warning: crop type mismatch detected for {crop}")
        print(elevator_location.loc[crop_type_check, 'crop type'])
        
    # Save location of each elevator
    elevator_location.to_file(output_dir / f"{crop}_elevator_location.geojson", driver="GeoJSON")
    elevator_location.to_excel(output_dir / f"{crop}_elevator_location.xlsx", index=False)
    print(f"{crop} elevator location is saved")


def determine_save_county_index_location(df_crop_price, df_crop_price_avg, crop, output_dir):
    indices = []
    tickers = []
    for col in df_crop_price.columns:
        if col.endswith("_Name"):
            ticker = col.split("_")[0]
            value = df_crop_price[col].dropna().unique()[0] if df_crop_price[col].notna().any() else None
            tickers.append(ticker)
            indices.append(value)

    # Create ticker + county location dataframe
    index_county_location = pd.DataFrame({
        "ticker": tickers,
        "index description": indices
    })
    index_county_location = index_county_location[index_county_location['ticker'].isin(df_crop_price_avg['ticker'])]
    index_county_location.reset_index(drop = True, inplace = True)
    
    # Extract county name from index description
    index_county_location['county'] = index_county_location['index description'].str.split(",").str[0].str.strip()
    # Replace "St" and "St." with "Saint" in county name
    index_county_location['county'] = index_county_location['county'].str.replace(r'\bSt\.?\b', 'Saint', regex=True)
    # Extract state name from index description
    index_county_location['state'] = index_county_location['index description'].str.split(",").str[1].str.split().str[0]
    # Apply geocode_location function to find the centroids of counties
    index_county_location[['geometry']] = (index_county_location['county'] + ', ' + index_county_location['state']).apply(lambda x: pd.Series(geocode_location_county(x)))
    # Extract crop type from price index description
    index_county_location['crop type'] = index_county_location['index description'].str.split(",").str[1].str.split().str[1].str.lower()
    # Replace "soybean" with "soybeans" in crop type
    index_county_location['crop type'] = index_county_location['crop type'].str.replace('soybean', 'soybeans', regex=False) 
    # Create a GeoDataFrame with counties location
    index_county_location = gpd.GeoDataFrame(index_county_location, geometry = "geometry", crs = "EPSG:4326")
    
    # Verify that the crop type listed for price index match the crop type we are currently working with
    crop_type_check = (index_county_location['crop type'].notna()) & (index_county_location['crop type'] != crop)
    if crop_type_check.any():
        print("Warning: crop type mismatch detected for {crop}")
        print(index_county_location.loc[crop_type_check, 'crop type'])
    
    # Save location of each crop price index county
    index_county_location.to_file(output_dir / f"{crop}_index_county_location.geojson", driver="GeoJSON")
    index_county_location.to_excel(output_dir / f"{crop}_index_county_location.xlsx", index=False)
    print(f"{crop} index county location is saved")
    
    

# Main script
for crop in crops:
    for level in price_levels:
        
        try:
            df_crop_price_merged = merge_crop_price_by_state(crop, states, level, input_dir)
            df_crop_price_merged_avg = clean_transform_save_merged_price_data(df_crop_price_merged, crop, level, drop_threshold, output_dir)
            
            if level == "elevator":
                determine_save_elevator_location(df_crop_price_merged, df_crop_price_merged_avg, crop, input_dir, output_dir)
            elif level == "county":
                determine_save_county_index_location(df_crop_price_merged, df_crop_price_merged_avg, crop, output_dir)
        
            print(f"All steps for {crop} at the {level} level are completed.")
        
        except Exception as e:
            print(f"Skipping {crop} and {level} because of: {e}")
            continue