# Import packages
import geopandas as gp
import pandas as pd
import requests 
from pathlib import Path
import zipfile


#---# Setting up the request to the website
def download_raw_html(
    html: str,
    raw_dir: Path,
    redownload: bool = True,
    size_tolerance: float = 0.05
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

    # send get request; set stream to true to ensure no large file issues
    response = requests.get(html, stream=True)
    response_size = int(response.headers.get("Content-Length", 0))

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

#---# Zip file extraction
def extract_zip_to_raw(zip_path: Path, raw_dir: Path) -> None:
    # Create a subdirectory named after the ZIP file (without extension)
    subfolder_name = zip_path.stem
    extract_path = raw_dir / subfolder_name
    extract_path.mkdir(parents=True, exist_ok=True)

    # Extract contents into the subfolder
    with zipfile.ZipFile(zip_path, 'r') as zip_ref:
        zip_ref.extractall(extract_path)


#---# Main execution for snakemake
if __name__ == "__main__":

    # unpack snakemake
    html = snakemake.params.html
    raw_dir = Path(snakemake.params.raw_dir)
    output_dir = Path(snakemake.params.output_dir)
    year_range = range(snakemake.params["year_range"][0], snakemake.params["year_range"][1])

    # download files
    downloaded_files = []

    for year in year_range:
        year_str = str(year)
        url = f"{html}/{year_str}_30m_cdls.zip"
        try:
            zip_path = download_raw_html(url, raw_dir)
            downloaded_files.append(zip_path)
        except Exception as e:
            print(f"Skipping {url}: {e}")

    # extract zip
    for zip_path in downloaded_files:
        extract_zip_to_raw(zip_path, raw_dir) 
