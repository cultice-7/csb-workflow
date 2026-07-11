# Which Proprietary Data to Download and Where to Place It

Several raw datasets are proprietary or otherwise not downloadable by the pipeline itself, and must be obtained manually from the project's protected fileshare (see [How to Get Access to the Fileshare](05-how-to-get-access-to-fileshare.md)) and placed into specific subfolders of `data/` (relative to the repository root) before running the corresponding Snakemake rules. Create these folders yourself if they don't already exist:

- **Regrow** — create `data/Regrow/2014-2024/` and `data/Regrow/2025/` folders. Place the legacy (2014–2024) field boundary GeoJSONs and Monitor_data exports in `data/Regrow/2014-2024/`, and the 2025 boundary GeoJSONs and 2025 Monitor data exports in `data/Regrow/2025/`.

- **DISES** — create `data/DISES/` folder. Place the survey export (`combined_data_clean.csv`) and the parcel shapefile (`DISES_All_Parcels_11.12.25.shp`) and all of its sidecar files (`.dbf`, `.shx`, `.prj`, etc.) directly in this folder.

- **gSSURGO** — create `data/gSSURGO/` folder. Inside it, place the national mukey grid folder (`FY2026_gSSURGO_mukey_grid/`, containing `MURASTER_30m.tif` and its sidecar files) and one `gSSURGO_{state}/` subfolder per study state, each containing that state's `gSSURGO_{state}.gdb` geodatabase.

- **Barchart elevator price data** — create `data/Grain Price/` folder. Place the raw per-state, per-crop Excel exports here (e.g. `{state}_corn_elevator_level.xlsx`, `{state}_corn_county_level.xlsx`, `{state}_soybeans_elevator_level.xlsx`, `{state}_soybeans_county_level.xlsx`, `{state}_wheat_elevator_level.xlsx` — wheat has no county-level file), plus any manually-curated elevator-location correction spreadsheet referenced by `clean_grain_price.py`.
