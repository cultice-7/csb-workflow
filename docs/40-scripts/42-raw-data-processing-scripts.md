# Raw Data Processing Scripts

**`reproject_elevation.py`**: Reprojects the merged national elevation mosaic to the project's equal-area CRS (EPSG:5070) via GDAL Warp.

**`calculate_slope.py`**: Computes slope in degrees from the reprojected elevation raster via GDAL DEM processing. Slope is deliberately computed *after* reprojection to an equal-area CRS, so the gradient calculation uses planar (not geographic lat/lon) units, avoiding the distortion that would result from differentiating directly over a geographic grid.

**`clip_elevation_slope.py`**: Clips the national reprojected elevation and slope rasters to each state boundary, with an outward buffer (`state_buffer_margin` in config) to retain context just past the state line.

**`clip_gSSURGO_mukey_rasters.py`**: Clips the national gSSURGO 30m mukey raster to each state boundary (same buffer parameter) to produce the per-state soil map-unit grids used later in supplement 8.

**`clip_cdl_rasters.py`**: Clips national CDL rasters per state per year (same buffer parameter), producing the `{state}_{year}_30m_cdls_clipped.tif` files used for Regrow's CDL validation. Missing-file years are skipped silently (by design, since not every year may be downloaded yet); any other processing error is re-raised rather than silently skipped.

**`clean_crop_prices.py`**: Cleans raw Barchart Excel elevator/county price spreadsheets into monthly average price series and geocodes elevator/county-index locations (via OpenStreetMap Nominatim, with a manually-curated correction file for known bad addresses) for use in supplement 7's nearest-neighbor price join. Drops price series (elevators or county indices) with fewer than a configured minimum number of non-null monthly observations.

**`clip_reproject_weather_rasters.py`**: Reprojects and clips monthly PRISM rasters to each state boundary, producing per-state/variable/month clipped GeoTIFFs used by supplement 6.

**`clip_csb_shape.py`**: Extracts each state's portion of the national CSB geodatabase via an attribute filter on `STATEFIPS`, saving one parquet per state per CSB year-window.

**`split_csb_shape_table.py`**: Splits each state's clipped CSB parquet into geometry-only, attributes-only, and combined shape+table outputs, reprojected to the target CRS.

**`check_csb_shape_table.py` / `check_regrow_shape_table.py`**: Detects pairs of heavily-overlapping (near-duplicate) field polygons within the same dataset via a self-overlay, flags the smaller of each overlapping pair, and saves a small lookup table. This matters specifically for rasterization: if two polygons cover nearly the same area, naive pixel-by-pixel rasterization would arbitrarily assign each pixel to whichever polygon happens to be processed last, producing inconsistent field assignment. The downstream rasterization scripts use this lookup to drop the "smaller" duplicate before rasterizing, then backfill its attributes from the "larger" field's results afterward.
