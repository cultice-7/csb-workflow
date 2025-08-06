# Import packages
import geopandas as gp
import pandas as pd
import requests
from pathlib import Path
import zipfile
import shutil
import tempfile

#---# Debugging placeholders
# year = str(2023)
# year7 = str(2023 - 7)
# html = f"https://www.nass.usda.gov/Research_and_Science/Crop-Sequence-Boundaries/datasets/NationalCSB_{year7}-{year}_rev23.zip"

#---# Setting up the request to the website
def download_raw_html(
    html: str,
    raw_dir: Path,
    redownload: bool = True,
    size_tolerance: float = 0.05
) -> None:

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
    response_size = int(response.headers["Content-Length"])

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
    raw_size = raw_path.stat().st_size
    if (raw_size - response_size)/response_size > size_tolerance:
        raw_path.unlink()
        raise RuntimeError("downloaded raw file not within set filesize tolerance; deleting raw file")

#---# Zip file extraction
def extract_zip_to_temp(zip_path: Path, output_dir: Path) -> None:
    temp_extract_dir = Path(tempfile.mkdtemp(dir=output_dir))

    # Extract all contents
    with zipfile.ZipFile(zip_path, 'r') as zip_ref:
        zip_ref.extractall(temp_extract_dir)

    # Locate the expected folder inside temp_extract
    extracted_root = next(temp_extract_dir.iterdir())

    if not extracted_root.is_dir():
        raise FileNotFoundError("Expected extracted folder not found.")

    return extracted_root

#---# Move all contents from extracted_root to output_dir
def move_extracted_contents(extracted_root: Path, output_dir: Path):
    for item in extracted_root.iterdir():
        target = Path(output_dir) / item.name
        if target.exists():
            if target.is_dir():
                shutil.rmtree(target)
            else:
                target.unlink()
        shutil.move(str(item), str(target))

#---# Main execution for snakemake
if __name__ == "__main__":

    # unpack snakemake
    html = snakemake.params.html
    raw_dir = snakemake.params.raw_dir
    output_dir = snakemake.params.output_dir

    # send html to download function
    download_raw_html(html, raw_dir)

    # define zip path
    zip_path = Path(raw_dir) / html.split("/")[-1]

    # extract zip to temp
    extracted_root = extract_zip_to_temp(zip_path, output_dir)

    # move contents to output directory
    move_extracted_contents(extracted_root, output_dir)

    # delete temp folder
    shutil.rmtree(extracted_root.parent, ignore_errors=True)