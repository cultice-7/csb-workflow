import geopandas as gp
import pandas as pd
from shapely.ops import nearest_points
from shapely.geometry import box

#---# Load required datasets
main_shape = gp.read_file("data/edited/Regrow/regrow_dises_supplemented_1.geojson")
roads = gp.read_file("data/Census/road/prisecroads.shp")

# Reproject vector datasets
target_crs = "EPSG:5070"
roads = roads.to_crs(target_crs)

#---# Clip roads to bounding box of main_shape with buffer
buffer_dist = 50000  # buffer 50,000 meters (50 kilometers)
minx, miny, maxx, maxy = main_shape.total_bounds
bbox = box(minx - buffer_dist, miny - buffer_dist, maxx + buffer_dist, maxy + buffer_dist)
roads_clipped = roads[roads.intersects(bbox)]

#---# Join supplementary data to main dataset
# Nearest distance to road
try:
    print("Calculating nearest road distances...")
    road_union = roads_clipped.geometry.union_all()
    nearest_points_data = []

    for idx, row in main_shape.iterrows():
        polygon_geom = row.geometry
        boundary_id = row['boundary_id']
        point_on_polygon, point_on_road = nearest_points(polygon_geom, road_union)

        # Distance to nearest road
        main_shape.at[idx, 'dist_to_road'] = point_on_polygon.distance(point_on_road)

        # Store nearest point coordinates
        nearest_points_data.append({
            'boundary_id': boundary_id,
            'polygon_x': point_on_polygon.x,
            'polygon_y': point_on_polygon.y,
            'road_x': point_on_road.x,
            'road_y': point_on_road.y
        })
except Exception as e:
    print(f"Error calculating distances: {e}")
    raise

#---# Save geojson and csv files
main_shape.to_file("data/edited/Regrow/regrow_dises_supplemented_2.geojson", driver="GeoJSON")
attribute_table = main_shape.drop(columns='geometry')
attribute_table.to_csv("data/edited/Regrow/regrow_dises_supplemented_2_table.csv", index=False)

#---# Save nearest points table
nearest_points_df = pd.DataFrame(nearest_points_data)
nearest_points_df.to_csv("data/edited/Regrow/regrow_points_on_road_table.csv", index=False)

#---# Generate nearest points on roads feature
try: 
    print("Generating nearest points on road...")
    point_on_road_gdf = gp.GeoDataFrame(
        nearest_points_df,
        geometry=gp.points_from_xy(nearest_points_df['road_x'], nearest_points_df['road_y']),
        crs=target_crs
    )
    point_on_road_gdf.to_file("data/edited/Regrow/regrow_points_on_road.geojson")
    print("Point feature saved successfully.")
except Exception as e:
    print(f"Error saving feature: {e}")
    raise