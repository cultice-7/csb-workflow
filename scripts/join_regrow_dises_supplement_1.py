import geopandas as gp
from rasterstats import zonal_stats

#---# Load required datasets
# Load main dataset
main_shape = gp.read_file("data/edited/Regrow/regrow_dises_spatialjoin.geojson")

# Projected raster paths
elevation_proj_path = "data/Geo/elevation/elevation_reproj.tif"
slope_proj_path = "data/Geo/elevation/slope_reproj.tif"

#---# Join supplementary data to main dataset
# Zonal statistics for elevation
try:
    print("Calculating mean elevation...")
    elevation_stats = zonal_stats(main_shape, elevation_proj_path, stats="mean", geojson_out=True)
    elevation_means = [feature['properties']['mean'] for feature in elevation_stats]
    main_shape['elevation_mean'] = elevation_means
    print("Mean elevation added to attribute table.")
except Exception as e:
    print(f"Error processing elevation: {e}")
    raise

# Zonal statistics for slope
try:
    print("Calculating mean slope...")
    slope_stats = zonal_stats(main_shape, slope_proj_path, stats="mean", geojson_out=True)
    slope_means = [feature['properties']['mean'] for feature in slope_stats]
    main_shape['slope_mean'] = slope_means
    print("Mean slope added to attribute table.")
except Exception as e:
    print(f"Error processing slope: {e}")
    raise

#---# Save geojson and csv files
main_shape.to_file("data/edited/Regrow/regrow_dises_supplemented_1.geojson", driver="GeoJSON")
attribute_table = main_shape.drop(columns='geometry')
attribute_table.to_csv("data/edited/Regrow/regrow_dises_supplemented_1_table.csv", index=False)

