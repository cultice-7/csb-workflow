# Import packages
import requests
from pathlib import Path
import rasterio
from rasterio.merge import merge

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
def merge_rasters(tif_files: list[Path], output_path: Path) -> None:
    src_files = [rasterio.open(str(fp)) for fp in tif_files]
    mosaic, out_trans = merge(src_files)

    out_meta = src_files[0].meta.copy()
    out_meta.update({
        "driver": "GTiff",
        "height": mosaic.shape[1],
        "width": mosaic.shape[2],
        "transform": out_trans,
        "count": mosaic.shape[0]
    })

    with rasterio.open(output_path, "w", **out_meta) as dest:
        dest.write(mosaic)

    print(f"Merged raster saved to {output_path}")

#---# Main execution for snakemake
if __name__ == "__main__":
    y_range = range(snakemake.params["y_range"][0], snakemake.params["y_range"][1])
    x_range = range(snakemake.params["x_range"][0], snakemake.params["x_range"][1])
    raw_dir = snakemake.params.raw_dir
    output_dir = snakemake.params.output_dir
    html = snakemake.params.html

    downloaded_files = []
    for y in y_range:
        for x in x_range:
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
        merge_rasters(downloaded_files, merged_path)
