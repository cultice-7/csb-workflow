import geopandas as gpd
import pandas as pd
from shapely.ops import nearest_points
import os
from pathlib import Path

# Import parameters from Snakemake
CSB_input_folder = snakemake.params.CSB_input_dir
CSB_output_folder = snakemake.params.CSB_output_dir
roads_input_folder = snakemake.params.roads_input_dir
roads_output_folder = snakemake.params.roads_output_dir
states = snakemake.params.states
CSB_years = snakemake.params.CSB_years
target_CRS = snakemake.params.target_CRS


#---# Load road datasets
roads = gpd.read_file(Path(roads_input_folder) / "prisecroads.shp")

# Reproject vector datasets
roads = roads.to_crs(target_CRS)

for year in CSB_years:
    for state in states:
        
        CSB_input_path_CSB = os.path.join(CSB_input_folder, f"{state}_CSB{year}_CSBID_geometry.parquet")
        output_path_spatial = os.path.join(CSB_output_folder, f"{state}_CSB{year}_nearest_roads_spatial.parquet")
        output_path_table = os.path.join(CSB_output_folder, f"{state}_CSB{year}_nearest_roads_table.parquet")
        road_output_path_geojson = os.path.join(roads_output_folder, f"{state}_CSB{year}_points_on_road.geojson")
        road_output_path_csv = os.path.join(roads_output_folder, f"{state}_CSB{year}_points_on_road_table.csv")
        
        Path(roads_output_folder).mkdir(parents=True, exist_ok=True)
        
        # Load CSB_shape datasets
        CSB_shape = gpd.read_parquet(CSB_input_path_CSB)
        
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
        
        # Convert float64 columns to float32 to save memory
        float64_cols = CSB_supplement_4.select_dtypes(include=["float64"]).columns
        CSB_supplement_4[float64_cols] = CSB_supplement_4[float64_cols].astype("float32")

        #---# Save files w/ and w/o geometry
        #CSB_supplement_4.to_parquet(output_path_spatial, compression="zstd")
        attribute_table = CSB_supplement_4.drop(columns='geometry')
        attribute_table.to_parquet(output_path_table, index=False, compression="zstd")

        #---# Save nearest points table
        nearest_points_df = pd.DataFrame(nearest_points_data)
        nearest_points_df.to_csv(road_output_path_csv, index=False)

        #---# Generate nearest points on roads feature
        try: 
            print("Generating nearest points on road...")
            point_on_road_gdf = gpd.GeoDataFrame(
                nearest_points_df,
                geometry=gpd.points_from_xy(nearest_points_df['road_x'], nearest_points_df['road_y']),
                crs=roads.crs
            )
            point_on_road_gdf.to_file(road_output_path_geojson, driver="GeoJSON")
            print("Point feature saved successfully.")
        except Exception as e:
            print(f"Error saving feature: {e}")
            raise