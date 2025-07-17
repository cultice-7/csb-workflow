#untested

import os
import requests
import zipfile
from io import BytesIO

# Base URL and output folder
base_url = "https://www2.census.gov/geo/tiger/TIGER2023/ROADS/"
output_folder = "census_roads_data"
os.makedirs(output_folder, exist_ok=True)

# State FIPS prefixes
prefixes = ["18", "26", "39"]

# Loop through county suffixes 001 to 199
for prefix in prefixes:
    for i in range(1, 200):  
        suffix = str(i).zfill(3)
        filename = f"tl_2023_{prefix}{suffix}_roads.zip"
        url = base_url + filename

        try:
            response = requests.get(url, timeout=10)
            if response.status_code == 200:
                with zipfile.ZipFile(BytesIO(response.content)) as z:
                    z.extractall(output_folder)
        except Exception:
            pass 

print("All matching ZIP files downloaded and extracted.")
