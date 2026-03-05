# Import packages
import requests
from pathlib import Path
import rasterio
from shapely.geometry import box
import geopandas as gpd
import math

#---# Download function
def download_raw_html(
    html: str,
    raw_dir: Path,
    redownload: bool = True,
    size_tolerance: float = 0.05,
    retries: int =3,
    delay: int = 5
) -> Path:
    
    # process the html to get out the raw filename
    raw_filename = html.split("/")[-1] # split the html and get out the last element

    # path for downloaded file
    raw_path = Path(raw_dir) / raw_filename # stack the path using Path object

    # ensure parent directory exists
    if not raw_path.parent.exists():
        raw_path.parent.mkdir()

    # if file already exists, print a message and delete based on redownload parameter
    if raw_path.exists():
        print(f"{raw_path.name} already exists")
        if redownload:
            print(f"deleting and redownloading {raw_path.name}...")
            raw_path.unlink()

    # send get request; set stream to true to ensure no large file issues; wait timeout = 10 secs
    response = requests.get(html, stream=True, timeout=10)
    response_size = int(response.headers["Content-Length"])

    # check if response is okay
    if response.ok and response.status_code == 200:


        # in error catching block...
        try:
            # begin writing file in chunks
            with open(raw_path, mode="wb") as file:
                for chunk in response.iter_content(chunk_size=10000 * 1024):
                    file.write(chunk)

        # if an error throws, print that an issue occurred and delete the broken file
        except Exception as e:
            raw_path.unlink(missing_ok=True)
            raise ValueError(f"Error downloading {raw_path.name}; deleting file and closing connection.")
        # after any condition, close the open http link
        finally:
            response.close()
    else:
        raise ValueError(
            f"""
            Issue with opening download link:
            Status code {response.status_code}
            """
        )

    # last check: make sure file size in response is close enough to filesize on disk
    raw_size = raw_path.stat().st_size
    if response_size and abs(raw_size - response_size) / response_size > size_tolerance:
        raw_path.unlink()
        raise RuntimeError("downloaded raw file not within set filesize tolerance; deleting raw file.")

    return raw_path

#---# Merge function
def merge_rasters_streaming(tif_files: list[Path], output_path: Path) -> None:
    # Read metadata from first file
    with rasterio.open(tif_files[0]) as src0:
        meta = src0.meta.copy()
        res_x = src0.transform.a
        res_y = -src0.transform.e

    # Compute mosaic bounding box (union of all tiles)
    minx = min(rasterio.open(f).bounds.left   for f in tif_files)
    miny = min(rasterio.open(f).bounds.bottom for f in tif_files)
    maxx = max(rasterio.open(f).bounds.right  for f in tif_files)
    maxy = max(rasterio.open(f).bounds.top    for f in tif_files)
    
    # Compute new mosaic transform and dimensions
    width  = math.ceil((maxx - minx) / res_x)
    height = math.ceil((maxy - miny) / res_y)
    transform = rasterio.Affine(res_x, 0, minx, 0, -res_y, maxy)
    
    meta.update({
        "driver": "GTiff",
        "height": height,
        "width": width,
        "transform": transform,
        "count": 1
    })

    # Create destination raster
    with rasterio.open(output_path, "w", **meta) as dst:
        # For each tile: read block-by-block and write into the mosaic
        for fp in tif_files:
            with rasterio.open(fp) as src:
                left, bottom, right, top = src.bounds
                # Compute placement of src tile in mosaic coordinates
                window = rasterio.windows.from_bounds(
                    left, bottom, right, top, transform=dst.transform
                )

                data = src.read(1)
                dst.write(data, 1, window=window)    

    print(f"Merged raster saved to {output_path}")

#---# Intersection between state-boundary polygon and DEM tile polygon
def tile_intersects_states(lat, lon, state_bound):
    tile = box(lon-1, lat-1, lon+1, lat+1)
    return state_bound.intersects(tile).any()    


#---# Main execution for snakemake
if __name__ == "__main__":
    y_range = range(snakemake.params["y_range"][0], snakemake.params["y_range"][1])
    x_range = range(snakemake.params["x_range"][0], snakemake.params["x_range"][1])
    raw_dir = snakemake.params.raw_dir
    output_dir = snakemake.params.output_dir
    html = snakemake.params.html
    states = snakemake.params.states
    
    df_state_bound = gpd.read_file("data/Census/state_bound/cb_2023_us_state_500k.shp")
    df_state_bound = df_state_bound[df_state_bound["STUSPS"].isin(states)]

    downloaded_files = []
    for y in y_range:
        for x in x_range:
            if tile_intersects_states(y, -x, df_state_bound):
                y_str = str(y).zfill(2)
                x_str = str(x).zfill(3)
                url = f"{html}/n{y_str}w{x_str}/USGS_1_n{y_str}w{x_str}.tif"
                try:
                    tif_path = download_raw_html(url, raw_dir)
                    downloaded_files.append(tif_path)
                except Exception as e:
                    print(f"Skipping {url}: {e}")

    if downloaded_files:
        merged_path = Path(output_dir) / "elevation.tif"
        merge_rasters_streaming(downloaded_files, merged_path)
