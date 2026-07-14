# Script Reference

This subfolder catalogs every active script in `scripts/` (excluding `scripts/archive/`, which holds superseded versions kept for reference only). For each script: what it produces and why, the key conceptual steps, and notable assumptions. Scripts are organized by pipeline stage to match the recommended execution chains described in [How to Obtain Output Datasets](../10-output/14-how-to-obtain-output-datasets.md).

> Several scripts listed here had known bugs identified during a systematic documentation review; almost all have since been corrected in the current codebase. See [Known Issues](../50-miscellaneous/52-known-issues.md) for the one item that remains open. The descriptions below reflect current, corrected behavior.

## Contents

| File | Pipeline stage | Scripts covered |
|---|---|---|
| [41 — Download Scripts](41-download-scripts.md) | Stage A (chain A) | `download_census_tract.py`, `download_state_bound.py`, `download_county_bound.py`, `download_roads.py`, `download_watershed.py`, `download_weather.py`, `download_csb.py`, `download_cdl.py`, `download_elevation.py` |
| [42 — Raw Data Processing Scripts](42-raw-data-processing-scripts.md) | Stage A–B | `reproject_elevation.py`, `calculate_slope.py`, `clip_elevation_slope.py`, `clip_gSSURGO_mukey_rasters.py`, `clip_cdl_rasters.py`, `clean_crop_prices.py`, `clip_reproject_weather_rasters.py`, `clip_csb_shape.py`, `split_csb_shape_table.py`, `check_*.py` |
| [43 — Regrow Preparation Scripts](43-regrow-preparation-scripts.md) | Stage C | `join_regrow_2025_updates.py`, `split_save_regrow_geometry.py`, `split_regrow_monitor_by_state.py`, `clean_regrow_table.py`, `join_regrow_shape_table.py` |
| [44 — Regrow Validation Scripts](44-regrow-validation-scripts.md) | Stage G (optional) | `rasterize_regrow_to_CDL_grid.py`, `join_regrow_cdl.py`, `validate_regrow_crop_with_cdl.py` |
| [45 — DISES Preparation Scripts](45-dises-preparation-scripts.md) | Stage B | `clean_dises_table.py`, `clean_dises_shape.py`, `join_dises_shape_table.py` |
| [46 — Regrow–DISES Merge Script](46-merge-regrow-dises-scripts.md) | Stage D | `join_regrow_dises.py` |
| [47 — Regrow Supplementary Data Merge Scripts](47-merge-regrow-supplements-scripts.md) | Stage D | `rasterize_regrow_to_gSSURGO_grid.py`, `join_regrow_census_tract.py`, `join_regrow_elevation_slope.py`, `join_regrow_watershed.py`, `join_regrow_nearest_roads.py`, `join_regrow_neighbor_field_mgmt.py`, `join_regrow_weather.py`, `join_regrow_crop_prices.py`, `cut_regrow_crop_prices.py`, `join_regrow_soil_composition.py`, `join_regrow_ag_census.py` |
| [48 — CSB–DISES Merge Script](48-merge-csb-dises-scripts.md) | Stage F | `join_csb_dises.py` |
| [49 — CSB Supplementary Data Merge Scripts](49-merge-csb-supplements-scripts.md) | Stage F | `rasterize_CSB_to_gSSURGO_grid.py`, `join_csb_census_tract.py`, `join_csb_elevation_slope.py`, `join_csb_watershed.py`, `join_csb_nearest_roads.py`, `join_csb_neighbor_field_mgmt.py`, `join_csb_weather.py`, `join_csb_crop_prices.py`, `cut_csb_crop_prices.py`, `join_csb_soil_composition.py`, `join_csb_ag_census.py` |
