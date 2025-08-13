import geopandas as gp
import rasterio
from rasterio.warp import calculate_default_transform, reproject, Resampling
from rasterstats import zonal_stats
from shapely.ops import nearest_points

#---# Load required datasets
# Load main dataset
regrow_shape = gp.read_file("data/edited/Regrow/regrow_dises_spatialjoin.geojson")

# Load supplementary vector data
roads = gp.read_file("data/Census/road/prisecroads.shp")

# Paths for supplementary raster data
elevation_path = "data/Geo/elevation/elevation.tif"
slope_path = "data/Geo/elevation/slope.tif"

#---# Prepare dataset
# Paths for reprojected raster data
elevation_proj_path = "data/Geo/elevation/elevation_epsg5070.tif"
slope_proj_path = "data/Geo/elevation/slope_epsg5070.tif"

# Reproject raster files (bilinear resampling for continous data)
def reproject_raster(input_path, output_path, dst_crs="EPSG:5070"):
    with rasterio.open(input_path) as src:
        transform, width, height = calculate_default_transform(
            src.crs, dst_crs, src.width, src.height, *src.bounds
        )
        kwargs = src.meta.copy()
        kwargs.update({
            'crs': dst_crs,
            'transform': transform,
            'width': width,
            'height': height
        })

        with rasterio.open(output_path, 'w', **kwargs) as dst:
            for i in range(1, src.count + 1):
                reproject(
                    source=rasterio.band(src, i),
                    destination=rasterio.band(dst, i),
                    src_transform=src.transform,
                    src_crs=src.crs,
                    dst_transform=transform,
                    dst_crs=dst_crs,
                    resampling=Resampling.bilinear
                )

reproject_raster(elevation_path, elevation_proj_path)
reproject_raster(slope_path, slope_proj_path)

# Reproject vector files
target_crs = "EPSG:5070"
regrow_shape = regrow_shape.to_crs(target_crs)
roads = roads.to_crs(target_crs)

#---# Join supplementary data to main dataset
# Zonal statistics for elevation
elevation_stats = zonal_stats(regrow_shape, elevation_proj_path, stats="mean", geojson_out=True)
elevation_means = [feature['properties']['mean'] for feature in elevation_stats]
regrow_shape['elevation_mean'] = elevation_means

# Zonal statistics for slope
slope_stats = zonal_stats(regrow_shape, slope_proj_path, stats="mean", geojson_out=True)
slope_means = [feature['properties']['mean'] for feature in slope_stats]
regrow_shape['slope_mean'] = slope_means

# Distance to nearest road
road_union = roads.geometry.union_all()
regrow_shape['dist_to_road'] = regrow_shape.geometry.apply(
    lambda geom: geom.distance(nearest_points(geom, road_union)[1])
)

# Save geojson and csv files
regrow_shape.to_file("data/edited/Regrow/regrow_dises_supplemented.geojson", driver="GeoJSON")
attribute_table = regrow_shape.drop(columns='geometry')
attribute_table.to_csv("data/edited/Regrow/regrow_dises_supplemented_table.csv", index=False)
