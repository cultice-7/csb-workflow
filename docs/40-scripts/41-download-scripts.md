# Download Scripts

**`download_census_tract.py` / `download_state_bound.py` / `download_county_bound.py`**
- *Overview*: Download single Census TIGER 2023 boundary zip files (tract, state, county respectively) and extract them.
- *Key steps*: build the download URL from config, stream-download with a file-size tolerance check, extract the zip.
- *Assumptions*: SSL verification is deliberately disabled for census.gov downloads only. Header access for `Content-Length` now safely defaults to `0` if absent (see [Known Issues](../50-miscellaneous/52-known-issues.md) for the remaining division-by-zero edge case on the size-tolerance check in these three scripts specifically).

**`download_roads.py`**
- *Overview*: Downloads TIGER 2023 PRISECROADS shapefiles per state FIPS code (including neighboring states, for edge coverage) and merges them into one national `prisecroads.shp`.
- *Key steps*: loop over configured state FIPS codes, download/extract each per-state zip, concatenate all shapefiles, write merged output.
- *Assumptions*: assumes consistent schema across all per-state shapefiles; no deduplication of roads that cross state borders.

**`download_watershed.py`**
- *Overview*: Downloads NHD High-Resolution shapefile bundles per (full-name) state, extracts HUC8/HUC10/HUC12 boundary layers, and merges each level across all states.
- *Key steps*: download per state, extract, recursively locate `WBDHU8/10/12.shp`, merge each level separately with duplicate-row removal.
- *Assumptions*: relies on exact, unique NHD shapefile naming across all extracted folders.

**`download_weather.py`**
- *Overview*: Downloads monthly PRISM climate raster zips (one per variable × year × month, from config) and saves extracted `.tif` files to a flat per-variable directory.
- *Key steps*: nested loop over variable/year/month, download/extract to a temp folder, copy resulting `.tif`s into the final output directory.
- *Assumptions*: one `.tif` per zip.

**`download_csb.py`**
- *Overview*: Downloads the national CSB geodatabase zip for each configured CSB year-window (currently `"1724"` → 2017–2024) plus its metadata HTML.
- *Key steps*: build URL from the year-window string, download, extract, unwrap any single wrapper folder.
- *Assumptions*: year-window strings are exactly 4 characters (split as first-2/last-2 digits).

**`download_cdl.py`**
- *Overview*: Downloads annual USDA CDL national raster zips for every year in the configured range and extracts each to its own per-year subfolder.

**`download_elevation.py`**
- *Overview*: Downloads 1-degree USGS 3DEP elevation tiles intersecting the study states, then mosaics them into a single `elevation.tif`.
- *Key steps*: build a lat/lon tile grid, test each tile's bounding box against the state boundaries, download matching tiles, stream-merge into one mosaic raster.
- *Assumptions*: tile-intersection test uses a slightly generous bounding box (not the exact tile footprint), so it may pull in a small number of unnecessary tiles — harmless, just slightly conservative.
