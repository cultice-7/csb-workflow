# CSB Supplementary Data Merge Scripts

All CSB supplement scripts loop over both `CSB_years` and `states`, producing one output file per state×year combination (named `{state}_CSB{year}_{name}_table.parquet`, where `{name}` is the supplement name). The primary key throughout is `CSBID`; field area is `CSBACRES`.

## rasterize_CSB_to_gSSURGO_grid.py

Rasterizes CSB field polygons onto the gSSURGO 30m mukey raster grid, one state×year combination at a time. Before rasterizing, it removes the smaller partner from each overlapping-field pair (identified from a pre-built overlapping-fields parquet) and sorts remaining geometries by area in ascending order so that, where polygon edges still touch a pixel center, the larger field wins. Each `CSBID` is mapped to a compact integer PID; the `CSBID`↔`pid` correspondence is saved as a separate parquet file. Geometries are reprojected to match the gSSURGO raster CRS, then burned onto the grid using the pixel-center rule (`all_touched=False`), with non-field pixels filled with −1. The output is one GeoTIFF per state×year in which each pixel holds the PID of the field whose polygon contains that pixel's center — used by `join_csb_soil_composition.py` to aggregate soil attributes across each field's footprint. Follows the same procedure as `rasterize_regrow_to_gSSURGO_grid.py`.

## join_csb_census_tract.py

Census tract join. Loads CSB field geometries and TIGER 2023 census tract boundaries, reprojects both to a common projected CRS, then computes the geographic **centroid** of each field polygon (a single representative point at the center of the polygon). Each centroid is matched to the census tract polygon it falls inside via a spatial within-join — this centroid approach is significantly faster than a full polygon-polygon area overlap. Fields whose centroid falls outside all tract boundaries (which can happen at polygon edges due to precision) are recovered via a nearest-tract fallback: the tract geographically closest to the unmatched centroid is assigned. Adds state/county/census tract identifiers and tract land/water area to each field. Follows the same procedure as `join_regrow_census_tract.py`.

## join_csb_elevation_slope.py

Elevation and slope join. Loads state-clipped, equal-area-projected elevation and slope rasters (outputs of `clip_elevation_slope.py`) and computes **zonal statistics** (mean) for both across each field polygon: every raster pixel whose center falls inside the polygon contributes to that field's mean. Attaches `elevation_mean` (meters) and `slope_mean` (degrees) to each field. Follows the same procedure as `join_regrow_elevation_slope.py`.

## join_csb_watershed.py

Watershed join. Loads USGS watershed boundary shapefiles at three nested scales — subbasin (HUC8), watershed (HUC10), and subwatershed (HUC12) — and processes each scale independently. For each scale, the boundaries are reprojected to the project CRS, each field's geographic centroid is computed, and a spatial within-join assigns each centroid to the watershed unit polygon it falls inside. There is no nearest-unit fallback: fields whose centroid falls outside all boundaries for a given scale are left with missing values for that scale. Each scale adds its own ID, name, and (for HUC10/12) type columns to the field table. Follows the same procedure as `join_regrow_watershed.py`.

## join_csb_nearest_roads.py

Nearest road distance join. Loads TIGER PRISECROADS (primary + secondary roads for the study area), reprojects to the project CRS, and finds the nearest road segment to each CSB field polygon using `gpd.sjoin_nearest`. This function returns the minimum distance from the field polygon to any point on the road geometry, plus the matched road's MTFCC class code (S1100 = primary, S1200 = secondary). When more than one road is equidistant, the one with higher road-class priority (primary over secondary) is kept. Adds `dist_to_road` (distance in CRS units, approximately meters in EPSG:5070) and `road_type`. Follows the same procedure as `join_regrow_nearest_roads.py`.

## join_csb_neighbor_field_mgmt.py

Neighboring-field crop composition. Derives land-management context from the CSB dataset itself — no external data source. Neighboring field pairs are identified by a spatial self-join of the full CSB field polygon layer (intersects predicate), with self-matches removed. For each field and each annual CDL crop column, computes the area-weighted share of neighboring field acreage in each crop group (corn / soybean / wheat / other), where each neighbor's contribution is weighted by its field area in acres (`CSBACRES`). Fields with no intersecting neighbors are left with missing values. Unlike the Regrow counterpart (`join_regrow_neighbor_field_mgmt.py`), no tillage or cover-crop neighbor averages are computed, since CSB has no tillage or cover-crop columns.

## join_csb_weather.py

Weather join. Loops over all PRISM AN Monthly raster files for 6 weather variables (precipitation, min/max temperature, mean dew point, min/max vapor pressure deficit), covering 2014–2025 at 800m resolution. For each field and each monthly raster, the PRISM value is extracted at the field's geographic centroid using `rasterio.Dataset.sample` — a single point sample at the center of the polygon, not a spatial average across the full polygon. Column names include `_centroid_` to reflect this. One column per variable per year-month is added to each field record. Follows the same procedure as `join_regrow_weather.py`.

## join_csb_crop_prices.py + cut_csb_crop_prices.py

Crop price join. Computes two types of monthly field-level price series.

**Elevator-level prices** (corn, soybeans, wheat): Loads geocoded elevator locations and their monthly average cash prices. Builds a KD-tree over elevator coordinates in projected space to efficiently find, for each field centroid, the N nearest elevators. Computes: (1) a gap-filled series — starting from the nearest elevator's monthly price, filling any missing months from the 2nd nearest, then 3rd, up to K; (2) a distance-weighted average price across all N nearest elevators, where each elevator's weight is its inverse Euclidean distance to the field centroid, applied month-by-month over non-missing observations only.

**County-level prices** (corn and soybeans; no county series for wheat): Same KD-tree approach but over county centroid locations. The field's own county is always placed first in the neighbor list — guaranteeing it receives the highest weight when available — and the remaining slots are filled by the nearest counties.

All per-crop temporary files are merged into a single `crop_prices_table.parquet` per state×year. `cut_csb_crop_prices.py` then drops the gap-filled single-nearest-elevator series, keeping only the N-nearest distance-weighted averages and the county-level series in the final `crop_prices_table_reduced.parquet`. Follows the same procedure as `join_regrow_crop_prices.py` + `cut_regrow_crop_prices.py`.

## join_csb_soil_composition.py

Soil join. Two-stage process.

**Stage 1 — mukey-level soil variable construction**: Reads six gSSURGO tabular layers from the state's `.gdb` file (mapunit, component, chorizon, corestrictions, muaggatt, legend/sacatalog) and produces one row per soil map unit (`mukey`). **Categorical attributes** (`drainagecl`, `cropprodindex`, `resdept_r`) are taken from the map unit's **dominant soil component** — the component with `majcompflag == "Yes"` and the highest `comppct_r` (percent composition) when multiple components qualify. **Continuous horizon attributes** (sand total, clay total, pH, organic matter) are depth-weighted across all soil horizons within each component to a target depth of 0–`soil_depth_cm` cm (default 30 cm), then component-weighted across all components within the mukey using `comppct_r` as the weight. Slope (`slopegraddcp`) is taken directly from the pre-aggregated muaggatt table at mukey resolution. Stage 1 follows the same procedure as `join_regrow_soil_composition.py`.

**Stage 2 — field-level aggregation**: Reads the gSSURGO mukey raster and the pre-built CSB field-ID raster (from `rasterize_CSB_to_gSSURGO_grid.py`) simultaneously. For every field-raster pixel belonging to a field, the mukey value from the corresponding mukey-raster pixel is recorded — producing for each field the complete list of mukeys present under its polygon footprint. Numeric soil variables are averaged and categorical variables take the mode across those mukeys. Fields excluded from rasterization (smaller overlapping-field partners) have their soil values copied from the corresponding larger field, keyed on `CSBID`.

## join_csb_ag_census.py

USDA Census of Agriculture join (county-level). Reads USDA NASS QuickStats CSVs for the 2017 and 2022 Census of Agriculture from four source files per state: land use (farmland/cropland/pastureland/woodland acreage, owned vs. rented), field crops (corn grain/silage, soybean, winter wheat harvested acreage), conservation practices (no-till, reduced-till, conventional-till, cover-crop acreage), and producer characteristics (count, sex, age, race, occupation, farming experience, decision-making role, residence). Specific data items are filtered from each file by matching on the `Data Item` and `Domain` columns, stacked into one long table, and pivoted to wide format — one row per county, one column per variable×census year. The resulting county-level table is matched to the list of counties actually present in that state's CSB census_tract output (`{state}_CSB{year}_census_tract_table.parquet`), keyed on `county_state_name`. The output is county-level — it is not merged onto individual field records. Follows the same procedure as `join_regrow_ag_census.py`.
