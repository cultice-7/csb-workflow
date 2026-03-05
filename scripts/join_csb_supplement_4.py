import geopandas as gpd
import pandas as pd
from shapely.ops import nearest_points
import os

#---# Load road datasets
roads = gpd.read_file("data/Census/road/prisecroads.shp")

# Reproject vector datasets
roads = roads.to_crs(epsg=5070)

# Input and output folders for CSB and road
input_folder_CSB = "data/edited/CSB/"
output_folder_CSB = "data/edited/CSB/"
output_folder_road = "data/edited/road/"

states = snakemake.params.states

for state in states:
    
    input_path_CSB = os.path.join(input_folder_CSB, f"{state}_CSB1724_clipped.gpkg")
    output_path_geojson = os.path.join(output_folder_CSB, f"{state}_CSB1724_supplement_4_spatial.geojson")
    output_path_csv = os.path.join(output_folder_CSB, f"{state}_CSB1724_supplement_4_table.csv")
    output_path_road_geojson = os.path.join(output_folder_road, f"{state}_CSB1724_points_on_road.geojson")
    output_path_road_csv = os.path.join(output_folder_road, f"{state}_CSB1724_points_on_road_table.csv")

    CSB_shape = gpd.read_file(input_path_CSB)
    
    # Reproject geometry to the same CRS (NAD83/CONUS Albers)
    CSB_shape = CSB_shape.to_crs(epsg=5070)
    
    # Keep only the columns necessary for the spatial join
    cols_to_keep = ['CSBID', 'geometry']
    CSB_shape = CSB_shape[cols_to_keep]


    #---# Nearest distance to road
    try:
        print(f"Calculating nearest road distances for {state}...")
        
        # Distance to nearest road
        CSB_supplement_4 = gpd.sjoin_nearest(CSB_shape, roads, how = "left", distance_col="dist_to_road")
        
        # Store nearest point coordinates
        nearest_points_data = []
        for idx, row in CSB_supplement_4.iterrows():
            field_id = row["CSBID"]
            field_geom = row.geometry
            road_geom = roads.loc[row.index_right].geometry

            # Get Shapely nearest points
            point_on_field, point_on_road = nearest_points(field_geom, road_geom)

            # Store results in a list of dicts
            nearest_points_data.append({
                'CSBID': field_id,
                'parcel_x': point_on_field.x,
                'parcel_y': point_on_field.y,
                'road_x': point_on_road.x,
                'road_y': point_on_road.y,
            })
        
        print(f"Calculating nearest road distances for {state} is complete")
        
        # S1100 stands for a primary road, S1200 stands for a secondary road
        priority = {"S1100": 1, "S1200": 2}
        CSB_supplement_4["road_priority"] = CSB_supplement_4["MTFCC"].map(priority)
        # Rename column RTTYP -> road type
        CSB_supplement_4.rename(columns={"RTTYP": "road_type"}, inplace=True)
        
        # If more than one road is linked, keep only the one with the highest priority
        CSB_supplement_4 = CSB_supplement_4.sort_values(["CSBID", "dist_to_road", "road_priority"], ascending=[True, True, True]).drop_duplicates('CSBID')
        
        # Drop columns with road information
        columns_to_drop = ['index_right', 'LINEARID', 'FULLNAME', 'MTFCC']
        CSB_supplement_4.drop(columns = columns_to_drop, inplace = True)
        
    except Exception as e:
        print(f"Error calculating distances: {e}")
        raise
    
    # Check whether there are any missing values
    cols_to_check = ['CSBID', 'dist_to_road']
    if CSB_supplement_4[cols_to_check].isna().any().any():
        print(CSB_supplement_4[CSB_supplement_4[cols_to_check].isna().any(axis=1)])

    #---# Save geojson and csv files
    #CSB_supplement_4.to_file(output_path_geojson, driver="GeoJSON")
    attribute_table = CSB_supplement_4.drop(columns='geometry')
    attribute_table.to_csv(output_path_csv, index=False)

    #---# Save nearest points table
    nearest_points_df = pd.DataFrame(nearest_points_data)
    nearest_points_df.to_csv(output_path_road_csv, index=False)

    #---# Generate nearest points on roads feature
    try: 
        print("Generating nearest points on road...")
        point_on_road_gdf = gpd.GeoDataFrame(
            nearest_points_df,
            geometry=gpd.points_from_xy(nearest_points_df['road_x'], nearest_points_df['road_y']),
            crs=roads.crs
        )
        point_on_road_gdf.to_file(output_path_road_geojson, driver="GeoJSON")
        print("Point feature saved successfully.")
    except Exception as e:
        print(f"Error saving feature: {e}")
        raise