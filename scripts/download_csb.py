# Import packages
import geopandas as gp
import pandas as pd
import requests
from pathlib import Path

#---# Debugging placeholders
# year = str(2023)
# year7 = str(2023 - 7)
# html = f"https://www.nass.usda.gov/Research_and_Science/Crop-Sequence-Boundaries/datasets/NationalCSB_{year7}-{year}_rev23.zip"

#---# Setting up the request to the website
def download_raw_html(
    html: str,
    redownload: bool | None = True,
    size_tolerance: float | None = 0.05
) -> None:

    # process the html to get out the raw filename
    raw_filename = html.split("/")[-1] # split the html and get out the last element
    raw_path = Path() / "data" / raw_filename # stack the path using Path object

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
    response_size = response.headers["Content-Length"]

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

if __name__ == "__main__":

    # unpack snakemake
    html = snakemake.params["html"]

    # send html to download function
    download_raw_html(html)