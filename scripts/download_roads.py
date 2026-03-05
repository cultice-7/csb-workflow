# Import packages
import geopandas as gp
import pandas as pd
import requests 
from pathlib import Path
import zipfile
import urllib3

# Suppress insecure request warning (due to disabling SSL verification, Census.gov only!)
urllib3.disable_warnings(urllib3.exceptions.InsecureRequestWarning)


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
    response = requests.get(html, stream=True, verify=False) #Disable SSL verification (Do this for census.gov only!)
    response_size = int(response.headers.get("Content-Length", 0))

    # check if response is okay and connection remains open
    response_is_ok = response.ok
    response_successful = (response.status_code == 200)
    response_connect_open = (response.headers["Connection"] == 'keep-alive')

    if response_is_ok and response_successful and response_connect_open:
        
        # in error catching block...
        try:
            # begin writing file in chunks
            with open(raw_path, mode="wb") as file:
                for chunk in response.iter_content(chunk_size= 10000 * 1024): # set chunk size to ensure download happens in pieces
                    file.write(chunk)

        # if an error throws, print that an issue occurred and delete the broken file
        except Exception as e:
            raw_path.unlink(missing_ok=True)
            raise ValueError(f"error in downloading {raw_path.name}; deleting file and closing connection")
        # after any condition, close the open http link
        finally:
            response.close()
    else:
        raise ValueError(
            f"""
            Issue with opening download link to html:
            Response ok? {response_is_ok}
            Successful exit code? {response_successful}
            Able to keep alive connection? {response_connect_open}
            """
        )

    # last check: make sure file size in response is close enough to filesize on disk
    if response_size > 0:
        raw_size = raw_path.stat().st_size
        if abs(raw_size - response_size)/response_size > size_tolerance:
            raw_path.unlink()
            raise RuntimeError("downloaded raw file not within set filesize tolerance; deleting raw file")
    return raw_path

#---# Zip file extraction
def extract_zip_to_raw(zip_path: Path, raw_dir: Path) -> None:
    # Extract all contents directly into raw_dir
    with zipfile.ZipFile(zip_path, 'r') as zip_ref:
        zip_ref.extractall(raw_dir)

#---# Merge function
def merge_vectors(shp_files: list[Path], output_path: Path) -> None:
    # Read all shapefiles and concatenate them
    gdfs = [gp.read_file(shp) for shp in shp_files]
    merged = pd.concat(gdfs, ignore_index=True)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    merged.to_file(output_path)
    print(f"Merged vector saved to {output_path}")


#---# Main execution for snakemake
if __name__ == "__main__":

    # unpack snakemake
    html = snakemake.params.html
    raw_dir = Path(snakemake.params.raw_dir)
    output_dir = Path(snakemake.params.output_dir)
    prefixes = [int(p) for p in snakemake.params.prefixes]

    # download files
    downloaded_files = []

    for prefix in prefixes:
        prefix_st = str(prefix).zfill(2)
        url = f"{html}/tl_2023_{prefix_st}_prisecroads.zip"
        try:
            zip_path = download_raw_html(url, raw_dir)
            downloaded_files.append(zip_path)
        except Exception as e:
            print(f"Skipping {url}: {e}")

    # extract zip
    for zip_path in downloaded_files:
        extract_zip_to_raw(zip_path, raw_dir) 

    # collect all shapefiles from raw_dir
    shp_files = list(raw_dir.glob("*.shp"))

    # define output path for merged shapefile
    merged_output_path = Path(snakemake.output.roads)

    # merge file and save the merged file in output_dir
    merge_vectors(shp_files, merged_output_path)

