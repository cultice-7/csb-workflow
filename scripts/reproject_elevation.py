import rasterio
from rasterio.warp import calculate_default_transform, reproject, Resampling
from osgeo import gdal

# Import parameters from Snakemake
dst_crs = snakemake.params.target_CRS
resampling_method = snakemake.params.resampling_method

def reproject_raster_gdal(input_path, output_path, dst_crs, resampling_method):
    """
    Reproject a raster using GDAL Warp

    Parameters
    ----------
    input_path : str
        Path to input raster
    output_path : str
        Path to output raster
    dst_crs : str
        Target CRS
    resampling_method : str
        GDAL resampling method
    """
    print(f"Reprojecting {input_path} to {output_path} ...") 
    try:
        warp_options = gdal.WarpOptions(
            dstSRS=dst_crs,
            resampleAlg=resampling_method,
            format="GTiff",
            multithread=True
        )

        gdal.Warp(
            output_path,
            input_path,
            options=warp_options
        )

        print(f"Finished reprojecting {input_path}.")

    except Exception as e:
        print(f"Error reprojecting {input_path}: {e}")

# Original raster paths
elevation_path = "data/Geo/elevation/elevation.tif"
# Projected raster paths
elevation_reproj_path = "data/Geo/elevation/elevation_reprojected.tif"

# Run reprojection
reproject_raster_gdal(
    input_path=elevation_path,
    output_path=elevation_reproj_path,
    dst_crs=dst_crs,
    resampling_method=resampling_method
)



"""
# Alternative reprojection approach using Rasterio

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
                        src_transform=src.transform,
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
"""