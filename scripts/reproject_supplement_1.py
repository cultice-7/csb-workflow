import rasterio
from rasterio.warp import calculate_default_transform, reproject, Resampling

# Original raster paths
elevation_path = "data/Geo/elevation/elevation.tif"
slope_path = "data/Geo/elevation/slope.tif"

# Projected raster paths
elevation_proj_path = "data/Geo/elevation/elevation_reproj.tif"
slope_proj_path = "data/Geo/elevation/slope_reproj.tif"

# Reproject raster files (bilinear resampling for continuous data)
def reproject_raster(input_path, output_path, dst_crs="EPSG:5070", resampling_method=Resampling.bilinear):
    print(f"Reprojecting {input_path} to {output_path} ...")
    try:
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
                                               src_crs=src.crs,
                        dst_transform=transform,
                        dst_crs=dst_crs,
                        resampling=resampling_method
                    )
        print(f"Finished reprojecting {input_path}.")
    except Exception as e:
        print(f"Error reprojecting {input_path}: {e}")

# Run reprojection
reproject_raster(elevation_path, elevation_proj_path)
reproject_raster(slope_path, slope_proj_path)
