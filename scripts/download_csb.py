# Import packages
import requests
from pathlib import Path
import zipfile
import shutil
import tempfile


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
        raw_path.parent.mkdir(parents=True, exist_ok=True)

    # if file already exists, print a message and delete based on redownload parameter
    if raw_path.exists():
        print(f"{raw_path.name} already exists")
        if redownload:
            print(f"deleting and redownloading {raw_path.name}...")
            raw_path.unlink()

    # send get request; set stream to true to ensure no large file issues
    response = requests.get(html, stream=True)
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
    raw_size = raw_path.stat().st_size
    if (raw_size - response_size)/response_size > size_tolerance:
        raw_path.unlink()
        raise RuntimeError("downloaded raw file not within set filesize tolerance; deleting raw file")

#---# Zip file extraction
def extract_zip_to_temp(zip_path: Path, output_dir: Path) -> Path:
    temp_extract_dir = Path(tempfile.mkdtemp(dir=output_dir))

    # Extract all contents
    with zipfile.ZipFile(zip_path, 'r') as zip_ref:
        zip_ref.extractall(temp_extract_dir)
        
    zip_path.unlink(missing_ok=True)

    return temp_extract_dir

#---# Move all contents from extracted_root to output_dir
def move_extracted_contents(temp_extract_dir: Path, output_dir: Path):
    output_dir.mkdir(parents=True, exist_ok=True)
    
    items = list(temp_extract_dir.iterdir())
    
    # If there is a single wrapper folder → unpack it
    if len(items) == 1 and items[0].is_dir():
        items = list(items[0].iterdir())

    # Move everything
    for item in items:
        target = output_dir / item.name

        # Remove existing target if needed
        if target.exists():
            if target.is_dir():
                shutil.rmtree(target)
            else:
                target.unlink()

        shutil.move(str(item), str(target))

    # Clean up temp directory
    shutil.rmtree(temp_extract_dir, ignore_errors=True)

#---# Main execution for Snakemake
if __name__ == "__main__":

    # Import parameters from Snakemake
    raw_dir = snakemake.params.raw_dir
    output_dir = snakemake.params.output_dir
    base_html = snakemake.params.base_html
    CSB_years = snakemake.params.CSB_years
    
    for CSB_year in CSB_years:

        # Design an html for a specific CSB_year
        html = f"{base_html}/NationalCSB_20{CSB_year[:2]}-20{CSB_year[2:]}_rev23.zip"
        
        # send html to download function
        download_raw_html(html, Path(raw_dir))

        # define zip path
        zip_path = Path(raw_dir) / html.split("/")[-1]

        # extract zip to temp
        temp_extract_dir = extract_zip_to_temp(zip_path, Path(output_dir))

        # move contents to output directory
        move_extracted_contents(temp_extract_dir, Path(output_dir))