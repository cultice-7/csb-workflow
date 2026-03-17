# Branch: `fix/snakefile-cleanup` - Change Summary

This document describes all changes made on the `fix/snakefile-cleanup` branch compared to `main`. The changes fall into three categories: **bug fixes**, **readability improvements**, and **cross-platform compatibility**.

**Files changed:** 19 (13 modified, 5 new, 1 renamed)

---

## Table of Contents

1. [New File Structure](#1-new-file-structure)
2. [Bug Fixes (11 items)](#2-bug-fixes)
3. [Readability Improvements](#3-readability-improvements)
4. [Cross-Platform Compatibility](#4-cross-platform-compatibility)
5. [Files Changed Summary](#5-files-changed-summary)

---

## 1. New File Structure

The single 657-line `Snakefile` has been split into 6 files organized by pipeline phase:

```
BEFORE:                          AFTER:
csb-workflow/                    csb-workflow/
├── Snakefile (657 lines)        ├── Snakefile (92 lines)
├── config.yml                   ├── rules/
├── envs/env.yml                 │   ├── download.smk  (153 lines)
└── scripts/                     │   ├── clean.smk     (127 lines)
                                 │   ├── clip.smk      (129 lines)
                                 │   ├── joins.smk     ( 57 lines)
                                 │   └── supplements.smk (345 lines)
                                 ├── config.yml
                                 ├── envs/env.yml
                                 └── scripts/
```

The main `Snakefile` now contains only configuration, constants, `rule all`, and `include:` statements. All rule definitions are in `rules/*.smk` organized by pipeline phase: download, clean, clip/reproject, spatial joins, and supplementary data joins.

---

## 2. Bug Fixes

### Bug 1: Backslash in gSSURGO Path

**File:** `scripts/clip_gSSURGO_rasters.py` (line 40) and `Snakefile`

Windows-style backslash path that fails on Linux and macOS.

```python
# BEFORE (broken on Linux/macOS):
input_file_path = Path("data/gSSURGO/FY2026_gSSURGO_mukey_grid\MURASTER_30m.tif")

# AFTER (works on all platforms):
input_file_path = Path("data/gSSURGO/FY2026_gSSURGO_mukey_grid/MURASTER_30m.tif")
```

### Bug 2: F-String Nested Double Quotes

**File:** `Snakefile` (lines 26, 38, 50, 156, 168) -> now in `rules/download.smk`

Nested double quotes inside double-quoted f-strings require Python 3.12+ (PEP 701). Changed to single outer quotes for broader compatibility.

```python
# BEFORE (requires Python 3.12+):
html = f"{config["census_tract"]["base_html"]}/cb_2023_us_tract_500k.zip"

# AFTER (works on all Python 3.6+):
html = f'{config["census_tract"]["base_html"]}/cb_2023_us_tract_500k.zip'
```

### Bug 3: Script Filename Case Mismatch

**File:** `scripts/download_3DEP.py` -> `scripts/download_3dep.py`

The Snakefile referenced `scripts/download_3dep.py` but the file was named `scripts/download_3DEP.py`. On case-sensitive filesystems (Linux), Snakemake could not find the script.

```
# BEFORE:
scripts/download_3DEP.py   (actual file)
script: "scripts/download_3dep.py"   (Snakefile reference)

# AFTER:
scripts/download_3dep.py   (renamed to match Snakefile)
```

### Bug 4: Typo `matadata`

**File:** `Snakefile` -> now in `rules/download.smk`

```python
# BEFORE:
matadata = "data/CSB/NationalCSB_2017-2024_rev23_metadata.htm"

# AFTER:
metadata = "data/CSB/NationalCSB_2017-2024_rev23_metadata.htm"
```

### Bug 5: `.geojson` vs `.parquet` Format Mismatch

**Files:** 7 scripts + Snakefile rule declarations

`join_regrow_shape_table.py` outputs `.parquet` files, but 7 downstream scripts tried to read `.geojson`. This broke the dependency chain.

| Script | Before | After |
|--------|--------|-------|
| `validate_regrow_shape.py` | `gp.read_file("...geojson")` | `gp.read_parquet("...parquet")` |
| `join_regrow_supplement_1.py` | `gpd.read_file("...geojson")` | `gpd.read_parquet("...parquet")` |
| `join_regrow_supplement_2.py` | `gpd.read_file("...geojson")` | `gpd.read_parquet("...parquet")` |
| `join_regrow_supplement_3.py` | `gpd.read_file("...geojson")` | `gpd.read_parquet("...parquet")` |
| `join_regrow_supplement_4.py` | `gpd.read_file("...geojson")` | `gpd.read_parquet("...parquet")` |
| `join_regrow_supplement_5.py` | `gpd.read_file("...geojson")` | `gpd.read_parquet("...parquet")` |
| `join_regrow_supplement_6.py` | `gpd.read_file("...geojson")` | `gpd.read_parquet("...parquet")` |

Example from `join_regrow_supplement_1.py`:

```python
# BEFORE:
input_path_Regrow = os.path.join(input_folder_Regrow, f"{state}_regrow_shape_table.geojson")
regrow_shape = gpd.read_file(input_path_Regrow)

# AFTER:
input_path_Regrow = input_folder_Regrow / f"{state}_regrow_shape_table.parquet"
regrow_shape = gpd.read_parquet(input_path_Regrow)
```

### Bug 6: Missing Weather Clipped Outputs

**File:** `Snakefile` -> now in `rules/clip.smk`

The `clip_reproject_weather_rasters` rule clipped all 9 weather variables but only declared 7 in its `output:` section. `soltotal` and `solslope` were missing.

```python
# BEFORE (only 7 of 9 variables declared):
output:
    ppt_clipped = expand(...)
    tmax_clipped = expand(...)
    tmin_clipped = expand(...)
    tmean_clipped = expand(...)
    tdmean_clipped = expand(...)
    vpdmax_clipped = expand(...)
    vpdmin_clipped = expand(...)
    # soltotal and solslope MISSING

# AFTER (all 9 variables declared):
output:
    ppt_clipped = expand(...)
    tmax_clipped = expand(...)
    tmin_clipped = expand(...)
    tmean_clipped = expand(...)
    tdmean_clipped = expand(...)
    vpdmax_clipped = expand(...)
    vpdmin_clipped = expand(...)
    soltotal_clipped = expand(...)   # ADDED
    solslope_clipped = expand(...)   # ADDED
```

### Bug 7: Weather Supplement 6 - Only 3 States Instead of 7

**File:** `Snakefile` -> now in `rules/supplements.smk`

The `join_regrow_supplement_6` and `join_csb_supplement_6` rules only declared clipped weather input rasters for 3 states (`["OH", "IN", "MI"]`), but the scripts process all 7 states. Fixed to expand for all 7 states.

```python
# BEFORE (only 3 states in input):
input:
    ppt_clipped = expand("...{state}...clipped.tif", state=["OH", "IN", "MI"])

# AFTER (all 7 states in input):
input:
    ppt_clipped = expand("...{state}...clipped.tif", state=STATES)
```

### Bug 8: Orphan `convert_csv_to_parquet` Rule

**File:** `Snakefile` -> now in `rules/supplements.smk`

This rule had no `input:` or `output:`, making it invisible to the DAG. Now documented as a manual utility rule with a clear comment explaining how to run it.

### Bug 9: `validate_regrow_shape.py` Missing Michigan (MI)

**File:** `scripts/validate_regrow_shape.py`

The script hardcoded 6 states and ignored `snakemake.params.states`. Michigan was never validated.

```python
# BEFORE (hardcoded 6 states, ignored Snakemake params):
states = ["OH", "MN", "WI", "IA", "IN", "IL"]   # MI missing!

# AFTER (uses Snakemake params - all 7 states including MI):
states = snakemake.params.states
```

Additionally converted from `os.path` to `pathlib.Path` and from `os.makedirs` to `Path.mkdir()`.

### Bug 10: Inconsistent `Grain Price` vs `Grain price` Casing

**Files:** `clean_grain_price.py`, `join_regrow_supplement_7.py`, `join_csb_supplement_7.py`

On case-sensitive filesystems (Linux), `data/edited/Grain price/` and `data/edited/Grain Price/` are different directories. The producer wrote to lowercase `price`, consumers read from both. Standardized to `Grain Price` everywhere.

```python
# BEFORE (inconsistent - broke on Linux):
elevator_location.to_file(f'data/edited/Grain price/{crop}_elevator_location.geojson')
# In other scripts:
elevator_location = gpd.read_file(f"data/edited/Grain price/{crop}_elevator_location.geojson")

# AFTER (consistent casing):
elevator_location.to_file(f'data/edited/Grain Price/{crop}_elevator_location.geojson')
elevator_location = gpd.read_file(f"data/edited/Grain Price/{crop}_elevator_location.geojson")
```

### Bug 11: `clip_csb_shape.py` Ignored Snakemake Params

**File:** `scripts/clip_csb_shape.py`

The script hardcoded FIPS codes and state names twice (once per CSB version) and ignored `snakemake.params.states`. Rewritten to use params with a FIPS lookup dictionary, eliminating code duplication.

```python
# BEFORE (hardcoded, duplicated, ignored params):
target_fips = ['17', '18', '19', '26', '27', '39', '55']
state_name =  ['IL', 'IN', 'IA', 'MI', 'MN', 'OH', 'WI']
# ... same block repeated for CSB1724 ...

# AFTER (uses snakemake.params, single FIPS mapping):
FIPS_TO_STATE = {'17': 'IL', '18': 'IN', '19': 'IA', '26': 'MI', '27': 'MN', '39': 'OH', '55': 'WI'}
states = snakemake.params.states
STATE_TO_FIPS = {v: k for k, v in FIPS_TO_STATE.items()}
target_fips = [STATE_TO_FIPS[s] for s in states]
```

---

## 3. Readability Improvements

### 3a. Constants Defined Once

The state list `["OH", "MN", "WI", "IA", "IN", "IL", "MI"]` was copy-pasted ~50 times throughout the old Snakefile with inconsistent ordering and quote styles. Now defined once in the main `Snakefile`:

```python
# BEFORE (repeated ~50 times with inconsistent ordering):
state=["OH", "MN", "WI", "IA", "IN", "IL", "MI"]    # in some rules
state=['IN', 'IL', 'MN', 'IA', 'OH', 'WI', 'MI']    # in other rules

# AFTER (defined once, referenced everywhere):
STATES = config["elevation"]["states"]
DISES_STATES = ["OH", "IN", "MI"]
YEARS = range(2014, 2025)
CSB_YEARS = ["1623", "1724"]
WEATHER_VARS = config["weather"]["weather_variables"]
MONTHS = [f"{m:02d}" for m in range(1, 13)]
```

### 3b. Complete `rule all`

```python
# BEFORE (only 7 intermediate outputs - most of the pipeline was unreachable):
rule all:
    input:
        tract = "data/Census/census_tract/cb_2023_us_tract_500k.shp",
        roads = "data/Census/road/prisecroads.shp",
        elevation = "data/Geo/elevation/elevation.tif",
        slope = "data/Geo/elevation/slope.tif",
        watershed = "data/Geo/watershed/watershed.shp",
        raw_cdl = expand("...", year=range(2014,2025)),
        validated_regrow_shape = expand("...", state=["OH", "MN", "WI", "IA", "IN", "IL"]),

# AFTER (all final pipeline outputs - triggers complete workflow):
rule all:
    input:
        # Regrow CDL validation (7 states)
        expand("data/edited/Regrow/validation/{state}_regrow_with_cdl_validation.geojson", state=STATES),
        expand("data/edited/Regrow/validation/{state}_regrow_validity_summary_by_year.csv", state=STATES),
        # Regrow-DISES spatial join (3 states)
        expand("data/edited/Regrow/{state}_regrow_dises_table.parquet", state=DISES_STATES),
        # CSB-DISES spatial join (3 states x 2 year ranges)
        expand("data/edited/CSB/{state}_CSB{years}_dises_table.parquet", state=DISES_STATES, years=CSB_YEARS),
        # Regrow supplement tables 1-7 (CSV) + 8 (Parquet)
        expand("data/edited/Regrow/{state}_regrow_supplement_{n}_table.csv", state=STATES, n=range(1, 8)),
        expand("data/edited/Regrow/{state}_regrow_supplement_8_table.parquet", state=STATES),
        # CSB supplement tables 1-7 (CSV) + 8 (Parquet)
        expand("data/edited/CSB/{state}_CSB1724_supplement_{n}_table.csv", state=STATES, n=range(1, 8)),
        expand("data/edited/CSB/{state}_CSB1724_supplement_8_table.parquet", state=STATES),
```

Running `snakemake` now triggers the complete pipeline to produce **142 final output files** (compared to only 7 intermediate files before):

| Final Output | States | Files | Format |
|-------------|--------|-------|--------|
| Regrow CDL validation (geojson) | 7 | 7 | GeoJSON |
| Regrow CDL validation summaries | 7 x 2 | 14 | CSV |
| Regrow-DISES join tables | 3 (OH, IN, MI) | 3 | Parquet |
| CSB-DISES join tables | 3 states x 2 year ranges | 6 | Parquet |
| Regrow supplement tables 1-7 | 7 states x 7 supplements | 49 | CSV |
| Regrow supplement table 8 (soil) | 7 | 7 | Parquet |
| CSB supplement tables 1-7 | 7 states x 7 supplements | 49 | CSV |
| CSB supplement table 8 (soil) | 7 | 7 | Parquet |
| **Total** | | **142** | |

Snakemake resolves the full DAG backwards from these 142 final outputs through all intermediate steps (downloading, cleaning, clipping, spatial joining) all the way to the raw data sources. Every rule in the pipeline is now reachable from `rule all`.

### 3c. Section Organization

Rules are organized by pipeline phase with consistent headers:

```python
# BEFORE (inconsistent dividers):
########### For running all rules
#---# Download datasets
########## Manage and clean downloaded data ########
#---#---#
#------------# Join supplement 1

# AFTER (consistent section headers in each .smk file):
# =============================================================================
# DOWNLOAD RULES
# Download external datasets
# =============================================================================
```

### 3d. Header Docstring

A 30-line header docstring was added to the main `Snakefile` explaining:
- What the pipeline does
- The two main datasets (Regrow + CSB) and 8 supplement layers
- Prerequisites (manual downloads required)
- How to run the pipeline

### 3e. Log Directives and Resource Hints

Every rule now has a `log:` directive and compute-heavy rules have `threads:` and `resources:`:

```python
# BEFORE:
rule clip_cdl_rasters:
    input: ...
    output: ...
    script: "scripts/clip_cdl_rasters.py"

# AFTER:
rule clip_cdl_rasters:
    input: ...
    output: ...
    log: "logs/clip_cdl_rasters.log"
    threads: 4
    resources:
        mem_mb = 8000
    script: "../scripts/clip_cdl_rasters.py"
```

### 3f. Removed Dead Code

Commented-out output declarations (e.g., `#..._shape = expand(...)`) scattered throughout the old Snakefile have been removed. They are preserved in git history on `main`.

---

## 4. Cross-Platform Compatibility

All modified scripts were converted from `os.path.join()` and `import os` to `pathlib.Path()`:

```python
# BEFORE:
import os
input_folder_Regrow = "data/edited/Regrow/"
input_path = os.path.join(input_folder_Regrow, f"{state}_regrow_shape_table.geojson")
os.makedirs("data/edited/Regrow/validation", exist_ok=True)

# AFTER:
from pathlib import Path
input_folder_Regrow = Path("data/edited/Regrow")
input_path = input_folder_Regrow / f"{state}_regrow_shape_table.parquet"
Path("data/edited/Regrow/validation").mkdir(parents=True, exist_ok=True)
```

Scripts that were not modified (because they had no bugs) were left as-is. Their existing `os.path.join()` and forward-slash string paths are already cross-platform compatible.

---

## 5. Files Changed Summary

| File | Change Type | What Changed |
|------|-------------|-------------|
| `Snakefile` | Rewritten | Split into modules, added constants, complete rule all, header docstring |
| `rules/download.smk` | New | 10 download rules (from old Snakefile) |
| `rules/clean.smk` | New | 9 clean/prepare rules (from old Snakefile) |
| `rules/clip.smk` | New | 7 clip/reproject rules + Bug 1, 6 fixes |
| `rules/joins.smk` | New | 3 validation & DISES join rules |
| `rules/supplements.smk` | New | 17 supplement rules + Bug 7, 8 fixes |
| `scripts/download_3DEP.py` | Renamed | -> `scripts/download_3dep.py` (Bug 3) |
| `scripts/clip_gSSURGO_rasters.py` | Modified | Backslash fix (Bug 1) |
| `scripts/validate_regrow_shape.py` | Modified | Added MI, snakemake.params, pathlib, .parquet (Bugs 5, 9) |
| `scripts/clip_csb_shape.py` | Modified | Integrated snakemake.params (Bug 11) |
| `scripts/clean_grain_price.py` | Modified | Grain Price casing (Bug 10) |
| `scripts/join_regrow_supplement_1.py` | Modified | .parquet + pathlib (Bug 5) |
| `scripts/join_regrow_supplement_2.py` | Modified | .parquet + pathlib (Bug 5) |
| `scripts/join_regrow_supplement_3.py` | Modified | .parquet + pathlib (Bug 5) |
| `scripts/join_regrow_supplement_4.py` | Modified | .parquet + pathlib (Bug 5) |
| `scripts/join_regrow_supplement_5.py` | Modified | .parquet + pathlib (Bug 5) |
| `scripts/join_regrow_supplement_6.py` | Modified | .parquet + pathlib (Bug 5) |
| `scripts/join_regrow_supplement_7.py` | Modified | Grain Price casing (Bug 10) |
| `scripts/join_csb_supplement_7.py` | Modified | Grain Price casing (Bug 10) |
