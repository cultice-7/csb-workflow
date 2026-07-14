# How to Obtain Output Datasets

## The Snakefile is the single source of truth for execution order

All scripts are wired together as Snakemake rules in `Snakefile`. Dependencies between rules are expressed through `input:`/`output:` file paths, so Snakemake automatically determines the correct execution order — you do not need to run scripts manually in sequence.

To build **everything** (both Regrow-based and CSB-based final outputs and all supplementary tables declared in `rule all`):

```bash
conda activate csb_workflow
snakemake all --cores <N>
```

`--cores <N>` tells Snakemake the maximum number of CPU cores it's allowed to use to run independent rules in parallel — it is not a unit of memory or time, just a concurrency limit. `N` must be a positive integer no larger than the number of logical cores on your machine. To check how many you have:

- **Windows**: Task Manager → Performance tab → CPU → "Logical processors", or run `echo %NUMBER_OF_PROCESSORS%` in Command Prompt / `$env:NUMBER_OF_PROCESSORS` in PowerShell.
- **macOS/Linux**: run `nproc` (Linux) or `sysctl -n hw.ncpu` (macOS) in a terminal.

> **Memory warning:** running the full DAG with a high `--cores` count can run multiple memory-heavy geospatial steps (rasterization, weather extraction, soil aggregation) in parallel and exhaust available RAM, which can crash VS Code or the shell session — this is a memory limit, not a core-count limit, so simply having many cores available does not mean you should use all of them. If you hit this, lower `--cores` (start with `--cores 1` or `2` and increase only if your machine has RAM headroom), or run the pipeline in stages using the chains below.

## Running a specific rule or chain of rules

Snakemake will automatically build all upstream dependencies of whatever rule you target. To build just one output:

```bash
snakemake <rule_name> --cores <N>
```

For example, to merge Regrow with the census tract data to get the `{state}_regrow_census_tract_table` output file:

```bash
snakemake join_regrow_census_tract --cores 1
```

## Recommended execution chains

**A. Supplementary/reference data (run once, shared by both branches):**
```
download_census_tract → download_county_bound → download_state_bound
download_elevation → reproject_elevation → calculate_slope → clip_elevation_slope
download_watershed
download_roads
download_weather → clip_reproject_weather_rasters
clean_crop_prices
clip_gSSURGO_mukey_rasters
```

**B. DISES data:**
```
clean_dises_table
clean_dises_shape
join_dises_shape_table
```

**C. Regrow data preparation:**
```
join_regrow_2025_updates
split_save_regrow_geometry
split_regrow_monitor_by_state
clean_regrow_table
join_regrow_shape_table
check_regrow_shape_table
```
(`rasterize_regrow_to_gSSURGO_grid` also belongs in this chain if you intend to run the soil composition join.)

**D. Regrow + DISES + supplementary joins:**
```
join_regrow_dises
join_regrow_census_tract
join_regrow_elevation_slope
join_regrow_watershed
join_regrow_nearest_roads
join_regrow_neighbor_field_mgmt
join_regrow_weather
join_regrow_crop_prices → cut_regrow_crop_prices
join_regrow_soil_composition  (requires rasterize_regrow_to_gSSURGO_grid first)
join_regrow_ag_census  (requires join_regrow_census_tract first)
```

**E. CSB data preparation:**
```
download_csb
clip_csb_shape
split_csb_shape_table
check_csb_shape_table
rasterize_CSB_to_gSSURGO_grid   (only needed for the soil composition join)
```

**F. CSB + DISES + supplementary joins:**
```
join_csb_dises
join_csb_census_tract
join_csb_elevation_slope
join_csb_watershed
join_csb_nearest_roads
join_csb_neighbor_field_mgmt
join_csb_weather
join_csb_crop_prices → cut_csb_crop_prices
join_csb_soil_composition
join_csb_ag_census  (requires join_csb_census_tract first)
```

**G. Regrow validation against CDL (optional, independent of the merge chains above):**
```
download_cdl
clip_cdl_rasters
rasterize_regrow_to_CDL_grid
join_regrow_cdl
validate_regrow_crop_with_cdl
```
