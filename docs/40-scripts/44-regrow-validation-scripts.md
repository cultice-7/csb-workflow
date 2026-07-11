# Regrow CDL Validation Scripts

**`rasterize_regrow_to_CDL_grid.py`**: Rasterizes Regrow field polygons onto the exact CDL raster grid (after dropping near-duplicate overlapping fields and verifying every configured CDL year shares an identical grid), assigning each field a simple integer ID (`pid`) burned into the raster by pixel-center containment.

**`join_regrow_cdl.py`**: For each state/year, extracts the modal CDL crop class under each field's rasterized footprint, backfills near-duplicate fields from their paired "larger" field, and flags `cdl_valid_{year}` based on whether the modal CDL class matches any of that field's Regrow-reported crop codes for that year (excluding the `999` non-cropland sentinel from the comparison).

**`validate_regrow_crop_with_cdl.py`**: Produces the final by-state, by-year, by-crop-category accuracy summary Excel files described in [Regrow Dataset](../20-datasets/21-regrow-dataset.md) ("How is Regrow data validated?").
