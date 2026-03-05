import geopandas as gpd
import pandas as pd
import numpy as np
from geopy.geocoders import Nominatim
from shapely.geometry import Point
import time

states = snakemake.params.states
crops = snakemake.params.crops


# Create a list of all months from 2010 to 2025
dates = pd.date_range(start="2010-01-01", end="2025-12-01", freq="MS")

for crop in crops:  
    # Create empty dataset with Date (month) column
    df_crop_elevator = pd.DataFrame({"Date": dates})
    df_crop_county = pd.DataFrame({"Date": dates})
    
    # Elevator-level prices
    for state in states:
        
        df_elevator_price = pd.read_excel(f"data/Grain Price/{state}_{crop}_elevator_level.xlsx", sheet_name="Values", header=[0, 1])
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
    # Drop the elevators with few observations: fewer than 12 observations between 2015 and 2024
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
    
    names = []
    tickers = []
    for col in df_crop_elevator.columns:
        if col.endswith("_Name"):
            ticker = col.split("_")[0] 
            value = df_crop_elevator[col].dropna().unique()[0] if df_crop_elevator[col].notna().any() else None
            tickers.append(ticker)
            names.append(value)

    # Create ticker + elevator location dataframe
    elevator_location = pd.DataFrame({
        "ticker": tickers,
        "elevator description": names
    })
    elevator_location = elevator_location[elevator_location['ticker'].isin(df_crop_avg['ticker'])]
    elevator_location.reset_index(drop = True, inplace = True)
    

    def geocode_location(location):
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

    # Extract elevator location from elevator description
    elevator_location['elevator location'] = elevator_location['elevator description'].str.split(";").str[1].str.strip()
    # Replace "St" and "St." with "Saint" in elevator location
    elevator_location['elevator location'] = elevator_location['elevator location'].str.replace(r'\bSt\.?\b', 'Saint', regex=True)
    # Replace incorrect elevator locations
    elevator_location_corrected = pd.read_excel(f"data/edited/Grain Price/corrected_elevator_location_corn_soybeans_wheat.xlsx")
    elevator_location['elevator location'] = elevator_location['ticker'].map(elevator_location_corrected.set_index('ticker')['corrected location']).fillna(elevator_location['elevator location'])
    # Apply geocode_location to find the centroids of towns where elevators are located
    elevator_location[['geometry', 'county', 'state']] = elevator_location['elevator location'].apply(lambda x: pd.Series(geocode_location(x)))
    # Extract crop type from elevator description
    elevator_location['crop type'] = elevator_location['elevator description'].str.split(";").str[2].str.split().str[0].str.lower()
    # Create a GeoDataFrame with elevator locations
    elevator_location = gpd.GeoDataFrame(elevator_location, geometry = "geometry", crs = "EPSG:4326")
    
    # Verify that the crop type listed for elevators match the crop type we are currently working with
    crop_type_check = (elevator_location['crop type'].notna()) & (elevator_location['crop type'] != crop)
    if crop_type_check.any():
        print("Warning: crop type mismatch detected")
        print(elevator_location.loc[crop_type_check, 'crop type'])
        
    # Save location of each elevator
    elevator_location.to_file(f'data/edited/Grain price/{crop}_elevator_location.geojson', driver="GeoJSON")
    elevator_location.to_excel(f"data/edited/Grain price/{crop}_elevator_location.xlsx", index=False)
    print(f"{crop}_elevator_location is saved")


    # County-level prices
    if crop != 'wheat':
        for state in states:
            
            df_county_price = pd.read_excel(f"data/Grain Price/{state}_{crop}_county_level.xlsx", sheet_name="Values", header=[0, 1])
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
        # Select only Date and columns with average county price indices
        df_crop_county_avg = df_crop_county.loc[:, df_crop_county.columns.str.endswith(('Date', 'Average'))]
        # Remove "Average" from column names
        df_crop_county_avg.columns = df_crop_county_avg.columns.str.replace('_Average', '', regex=False)
        # Select values from January 2014 onwards
        df_crop_county_avg = df_crop_county_avg[df_crop_county_avg['Date'] >= '2014-01-01']
        df_crop_county_avg.reset_index(drop = True, inplace = True)
        # Drop the counties with few observations: fewer than 12 observations between 2015 and 2024
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
        
        indices = []
        tickers = []
        for col in df_crop_county.columns:
            if col.endswith("_Name"):
                ticker = col.split("_")[0]
                value = df_crop_county[col].dropna().unique()[0] if df_crop_county[col].notna().any() else None
                tickers.append(ticker)
                indices.append(value)

        # Create ticker + county location dataframe
        index_county_location = pd.DataFrame({
            "ticker": tickers,
            "index description": indices
        })
        index_county_location = index_county_location[index_county_location['ticker'].isin(df_crop_county_avg['ticker'])]
        index_county_location.reset_index(drop = True, inplace = True)
        
        def geocode_location(location):
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

        # Extract county name from index description
        index_county_location['county'] = index_county_location['index description'].str.split(",").str[0].str.strip()
        # Replace "St" and "St." with "Saint" in county name
        index_county_location['county'] = index_county_location['county'].str.replace(r'\bSt\.?\b', 'Saint', regex=True)
        # Extract state name from index description
        index_county_location['state'] = index_county_location['index description'].str.split(",").str[1].str.split().str[0]
        # Apply geocode_location function to find the centroids of counties
        index_county_location[['geometry']] = (index_county_location['county'] + ', ' + index_county_location['state']).apply(lambda x: pd.Series(geocode_location(x)))
        # Extract crop type from price index description
        index_county_location['crop type'] = index_county_location['index description'].str.split(",").str[1].str.split().str[1].str.lower()
        # Replace "soybean" with "soybeans" in crop type
        index_county_location['crop type'] = index_county_location['crop type'].str.replace('soybean', 'soybeans', regex=False) 
        # Create a GeoDataFrame with counties location
        index_county_location = gpd.GeoDataFrame(index_county_location, geometry = "geometry", crs = "EPSG:4326")
        
        # Verify that the crop type listed for price index match the crop type we are currently working with
        crop_type_check = (index_county_location['crop type'].notna()) & (index_county_location['crop type'] != crop)
        if crop_type_check.any():
            print("Warning: crop type mismatch detected")
            print(index_county_location.loc[crop_type_check, 'crop type'])
        
        # Save location of each crop price index county
        index_county_location.to_file(f'data/edited/Grain price/{crop}_index_county_location.geojson', driver="GeoJSON")
        index_county_location.to_excel(f"data/edited/Grain price/{crop}_index_county_location.xlsx", index=False)
        print(f"{crop}_index_county_location is saved")