import geopandas as gpd
import pandas as pd
from shapely.ops import nearest_points
import os

#---# Load road dataset
roads = gpd.read_file("data/Census/road/prisecroads.shp")

# Reproject vector datasets
roads = roads.to_crs(epsg=5070)

# Input and output folders for Regrow and road
input_folder_Regrow = "data/edited/Regrow/"
output_folder_Regrow = "data/edited/Regrow/"
output_folder_road = "data/edited/road/"

states = snakemake.params.states

for state in states:
    
    input_path_Regrow = os.path.join(input_folder_Regrow, f"{state}_regrow_shape_table.geojson")
    output_path_geojson = os.path.join(output_folder_Regrow, f"{state}_regrow_supplement_4_spatial.geojson")
    output_path_csv = os.path.join(output_folder_Regrow, f"{state}_regrow_supplement_4_table.csv")
    output_path_road_geojson = os.path.join(output_folder_road, f"{state}_regrow_points_on_road.geojson")
    output_path_road_csv = os.path.join(output_folder_road, f"{state}_regrow_points_on_road_table.csv")
    
    #---# Load regrow_shape datasets
    regrow_shape = gpd.read_file(input_path_Regrow)
    
    # Reproject geometry to the same CRS (NAD83/CONUS Albers)
    regrow_shape = regrow_shape.to_crs(epsg=5070)
    
    # Keep only the columns necessary for the spatial join
    cols_to_keep = ['field_id', 'geometry']
    regrow_shape = regrow_shape[cols_to_keep]

    #---# Nearest distance to road
    try:
        print(f"Calculating nearest road distances for {state}...")

        # Distance to nearest road
        regrow_supplement_4 = gpd.sjoin_nearest(regrow_shape, roads, how = "left", distance_col = "dist_to_road")
        
        # Store nearest point coordinates
        nearest_points_data = []
        for idx, row in regrow_supplement_4.iterrows():
            field_id = row["field_id"]
            field_geom = row.geometry
            road_geom = roads.loc[row.index_right].geometry

            # Get Shapely nearest points
            point_on_field, point_on_road = nearest_points(field_geom, road_geom)

            # Store results in a list of dicts
            nearest_points_data.append({
                'field_id': field_id,
                'parcel_x': point_on_field.x,
                'parcel_y': point_on_field.y,
                'road_x': point_on_road.x,
                'road_y': point_on_road.y,
            })
        
        print(f"Calculating nearest road distances for {state} is complete")
        
        # S1100 stands for a primary road, S1200 stands for a secondary road
        priority = {"S1100": 1, "S1200": 2}
        regrow_supplement_4["road_priority"] = regrow_supplement_4["MTFCC"].map(priority)
        # Rename column RTTYP -> road type
        regrow_supplement_4.rename(columns={"RTTYP": "road_type"}, inplace=True)
        
        # If more than one road is linked, keep only the one with the highest priority
        regrow_supplement_4 = regrow_supplement_4.sort_values(["field_id", "dist_to_road", "road_priority"], ascending=[True, True, True]).drop_duplicates('field_id')
        # Drop columns with road information
        columns_to_drop = ['index_right', 'LINEARID', 'FULLNAME', 'MTFCC']
        regrow_supplement_4.drop(columns = columns_to_drop, inplace = True)

    except Exception as e:
        print(f"Error calculating distances: {e}")
        raise

    # Check whether there are any missing values
    cols_to_check = ['field_id', 'dist_to_road']
    if regrow_supplement_4[cols_to_check].isna().any().any():
        print(regrow_supplement_4[regrow_supplement_4[cols_to_check].isna().any(axis=1)])
    
    #---# Save geojson and csv files
    #regrow_supplement_4.to_file(output_path_geojson, driver="GeoJSON")
    attribute_table = regrow_supplement_4.drop(columns='geometry')
    attribute_table.to_csv(output_path_csv, index = False)

    #---# Save nearest points table
    nearest_points_df = pd.DataFrame(nearest_points_data)
    nearest_points_df.to_csv(output_path_road_csv, index = False)

    #---# Generate nearest points on roads feature
    try: 
        print("Generating nearest points on road...")
        point_on_road_gdf = gpd.GeoDataFrame(
            nearest_points_df,
            geometry = gpd.points_from_xy(nearest_points_df['road_x'], nearest_points_df['road_y']),
            crs = roads.crs
        )
        point_on_road_gdf.to_file(output_path_road_geojson, driver="GeoJSON")
        print("Point feature saved successfully.")
    except Exception as e:
        print(f"Error saving feature: {e}")
        raise