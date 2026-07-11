# csb-workflow

**Land-Based Carbon SAS Project — Farmer Decision Model**

This repository implements the SAS Project data pipeline: a Snakemake-orchestrated workflow that builds agricultural field-level datasets describing land management practices, farmer characteristics, and farmland conditions across a 7-state Midwest study region for a multi-year period.

> **Status of this document:** This is a draft of README document supporting the code pipepline developed for SAS project. The file is assembled from the project README (10Mar2026), the representative-field methodology memo, the script order memo, the Snakefile, and a full read-through of every active script in `scripts/`.

---

## Table of Contents

1. [Project Overview](#1-project-overview)
2. [Organization Structure of Output Datasets](#2-organization-structure-of-output-datasets)
3. [How to Obtain Output Datasets](#3-how-to-obtain-output-datasets)
4. [Dataset Documentation (FAQ)](#4-dataset-documentation-faq)
   - [4.1 Regrow Dataset](#41-regrow-dataset)
   - [4.2 USDA Crop Sequence Boundaries (CSB) Dataset](#42-usda-crop-sequence-boundaries-csb-dataset)
   - [4.3 DISES Data](#43-dises-data)
   - [4.4 Supplementary Farmland Characteristics](#44-supplementary-farmland-characteristics)
   - [4.5 Merging Regrow, DISES, and Supplementary Data](#45-merging-regrow-dises-and-supplementary-data)
   - [4.6 Merging CSB, DISES, and Supplementary Data](#46-merging-csb-dises-and-supplementary-data)
5. [Documentation](#5-documentation)
6. [Explanation of Individual Scripts](#6-explanation-of-individual-scripts)
7. [Maumee Watershed Filtering](#7-maumee-watershed-filtering)
8. [Known Issues](#8-known-issues)

---

## 1. Project Overview

### 1.1 What this project does

This project provides field-level data on land management activities carried out by individual farmers, farmer characteristics, and supplementary data that may influence farmers' decision-making regarding adoption of specific land management practices.

#### This pipeline develops two parallel datasets: Regrow-based and CSB-based datasets

There is no single dataset that offers both rich land-management detail and unrestricted public distribution. **Regrow** is a proprietary, remote-sensing-derived dataset with detailed tillage, cover crop, and main-crop information, but its distribution is restricted. **CSB** is a public USDA dataset with main-crop information only, but no licensing restrictions. Rather than committing to a single source, this project develops the pipeline as **two parallel branches**, one built around Regrow fields and one built around CSB fields. Each branch is independently merged with the same DISES farmer-survey data and the same supplementary farmland characteristics, using a mirrored (but separately implemented) set of scripts. This lets the team compare results across both field-delineation sources, and choose — or combine — whichever is more suitable for a given downstream analysis, rather than being locked into the limitations of either source alone.

- Land management activity data are dynamic from **2014–2025** (Regrow) or **2017–2024** (CSB).
- Farmer characteristics (DISES survey) are available for a **single reference year (2023 season, collected early 2024)**.
- Supplementary farmland characteristics cover 8 topics: census geography, elevation/slope, watershed geography, road proximity, neighboring-field land management, weather, crop price, and soil composition — see [4.4](#44-supplementary-farmland-characteristics) for the full breakdown.

The project produces **two parallel, independently usable output datasets**, built around two different field-delineation sources:

| | Regrow-based dataset | CSB-based dataset |
|---|---|---|
| Field unit | "Regrow field" (proprietary ML-delineated boundary) | "CSB field" (USDA algorithm-delineated boundary, from CDL) |
| Land management detail | Tillage, cover crop, main crop, full cultivation-cycle timeline | Main crop rotation only |
| License | Proprietary (Regrow Ag) | Public (USDA NASS/ERS) |
| Primary key | `field_id` | `CSBID` |
| Timeline | 2014–2025 | 2017–2024 |

Both datasets are merged with DISES farmer-survey data and the same 8 blocks of supplementary farmland characteristics, using a parallel but separately-implemented set of scripts (the Regrow-side and CSB-side scripts mirror each other in structure — see [Section 6](#6-explanation-of-individual-scripts)).

### 1.2 Technological prerequisites

- **Anaconda / Miniconda** — for environment management.
- **Python** — managed via the conda environment below;
- **Git Bash** (Windows) or any POSIX-compatible shell (macOS/Linux) — for running `git`/command-line steps below.
- **Snakemake** — workflow orchestration; every processing step in this repository is a Snakemake rule.
- **Visual Studio Code** (or any editor) — used by the project authors to write and run the code; not strictly required, any Python-capable editor works.
- Core scientific/geospatial stack: `geopandas`, `pandas`, `numpy`, `rasterio`, `rasterstats`, `shapely`, `gdal`, `requests` — all pinned in `envs/env.yml`.

### 1.3 Cloning the repository and setting up the environment

#### Option A: Git command line (Git Bash / terminal)

1. Install Git if you don't already have it. A helpful guide covering three key steps is available here: [Git and GitHub Setup Guide](https://dev.to/sharon_m/installing-git-bash-on-windows-a-beginners-guide-2k9p). The guide walks you through: (1) How to install Git Bash on your computer; (2) Connecting Git Bash to GitHub; and (3) Understanding Git version control.
2. Open Git Bash (Windows) or your terminal (macOS/Linux).
3. Navigate to the folder where you want the project to live:
   ```bash
   cd /path/to/your/projects/folder
   ```
4. Clone the repository:
   ```bash
   git clone <repository-url>
   ```
5. Move into the cloned folder:
   ```bash
   cd csb-workflow
   ```

#### Option B: Visual Studio Code

1. Open VS Code.
2. Open the Command Palette (`Ctrl+Shift+P`) and run **Git: Clone**.
3. Paste the repository URL when prompted and select a local destination folder.
4. When prompted, choose **Open** to open the cloned repository in VS Code.

#### Setting up the conda environment

From an Anaconda Prompt (or any shell with `conda` on the PATH), with your working directory set to the cloned `csb-workflow` folder:

```bash
# Navigate to the cloned repository folder
cd /path/to/your/csb-workflow

# Create the environment from the provided spec (first-time setup only)
conda env create -f envs/env.yml

# Activate the environment
conda activate csb_workflow
```

If you've already created the environment previously and `envs/env.yml` has since changed, update it instead of recreating it:

```bash
conda env update -f envs/env.yml --prune
```

The environment name is `csb_workflow` (defined in `envs/env.yml`). All Snakemake commands in this README assume this environment is active.

#### Accessing the protected project fileshare

> **[TEAMMATE TO FILL IN]**: Add step-by-step instructions here for how to request/obtain access to the project's protected fileshare (e.g., who to contact, required permissions or accounts, VPN steps, link to the fileshare itself), since the proprietary raw datasets referenced below must be downloaded from it manually.

#### Placing manually-downloaded raw data

Several raw datasets are proprietary or otherwise not downloadable by the pipeline itself, and must be obtained manually from the project's protected fileshare and placed into specific subfolders of `data/` (relative to the repository root) before running the corresponding Snakemake rules. Create these folders yourself if they don't already exist:

- **Regrow** — create `data/Regrow/2014-2024/` and `data/Regrow/2025/` folders. Place the legacy (2014–2024) field boundary GeoJSONs and Monitor_data exports in `data/Regrow/2014-2024/`, and the 2025 boundary GeoJSONs and 2025 Monitor data exports in `data/Regrow/2025/`.
- **DISES** — create `data/DISES/` folder. Place the survey export (`combined_data_clean.csv`) and the parcel shapefile (`DISES_All_Parcels_11.12.25.shp`) and all of its sidecar files (`.dbf`, `.shx`, `.prj`, etc.) directly in this folder.
- **gSSURGO** — create `data/gSSURGO/` folder. Inside it, place the national mukey grid folder (`FY2026_gSSURGO_mukey_grid/`, containing `MURASTER_30m.tif` and its sidecar files) and one `gSSURGO_{state}/` subfolder per study state, each containing that state's `gSSURGO_{state}.gdb` geodatabase.
- **Barchart elevator price data** — create `data/Grain Price/` folder. Place the raw per-state, per-crop Excel exports here (e.g. `{state}_corn_elevator_level.xlsx`, `{state}_corn_county_level.xlsx`, `{state}_soybeans_elevator_level.xlsx`, `{state}_soybeans_county_level.xlsx`, `{state}_wheat_elevator_level.xlsx` — wheat has no county-level file), plus any manually-curated elevator-location correction spreadsheet referenced by `clean_grain_price.py`.
---

## 2. Organization Structure of Output Datasets

The pipeline produces output in two parallel branches — one based on **Regrow**, one based on **CSB** — each centered on its own field-delineation unit. Within each branch, most merge operations produce their own standalone output file (rather than one monolithic table), so that users can pull in only the supplementary blocks they need.

### 2.1 Regrow-based outputs

| Script(s) that produce the file | Output file pattern | Contents |
|---|---|---|
| `split_save_regrow_geometry` | `{state}_regrow_fieldID_geometry.parquet` / `.gpkg` | `field_id` + geometry only |
| `join_regrow_shape_table` | `{state}_regrow_table.parquet` | All Regrow attribute columns (tillage, cover crop, main crop, cultivation cycle dates/confidence), no geometry |
| `join_regrow_dises` | `{state}_regrow_dises_table.parquet` | DISES survey variables (suffix `_dises`) + match-quality + representative-field indicators, joined to each Regrow field |
| `join_regrow_supplement_1` … `join_regrow_supplement_8` (+ `cut_regrow_supplement_7`) | `{state}_regrow_supplement_{1-8}_table.parquet` (supplement 7 also has a `_table_reduced` variant) | One file per supplementary data block (see [2.3](#23-supplementary-data-codes)). The `_table_reduced` variant of supplement 7 exists because the full file contains nearest-single-elevator price data that is sensitive and not cleared for distribution — `cut_regrow_supplement_7.py` drops those columns before the reduced file is shared (see the corresponding script for exactly which columns are removed) |
| `join_regrow_supplement_9` | `{state}_regrow_supplement_9_table.parquet` | **County-level**, not field-level — see [2.3](#23-supplementary-data-codes) and [4.4](#44-supplementary-farmland-characteristics) |
| `join_regrow_cdl` / `validate_regrow_crop_with_cdl` | `{state}_regrow_cdl_validation.parquet`, `{state}_regrow_cdl_summary_by_crop_category.xlsx` | Internal validation of Regrow main-crop calls against USDA CDL — not a "supplementary data" file per se, but a QA output |

States covered: `IA, IL, IN, MI, MN, OH, WI`. DISES join is only produced for `IN, MI, OH` (the states where DISES survey data exists).

### 2.2 CSB-based outputs

Built around the **CSB field** (`CSBID`) instead of the Regrow field, with an added `{years}` token in every filename (currently `1724`, i.e. the 2017–2024 CSB release window):

| Script(s) | Output file pattern | Contents |
|---|---|---|
| `split_csb_shape_table` | `{state}_CSB{years}_CSBID_geometry.parquet` / `.gpkg` | `CSBID` + geometry only |
| `split_csb_shape_table` | `{state}_CSB{years}_table.parquet` | All CSB attribute columns (`CSBACRES` + one `CDL{year}` main-crop column per year in the release window, `CDL2017`…`CDL2024`), no geometry |
| `join_csb_dises` | `{state}_CSB{years}_dises_table.parquet` | DISES survey variables (suffix `_dises`) + match-quality + representative-field indicators, joined to each CSB field |
| `join_csb_supplement_1` … `join_csb_supplement_8` (+ `cut_csb_supplement_7`) | `{state}_CSB{years}_supplement_{1-8}_table.parquet` | One file per supplementary data block (see [2.3](#23-supplementary-data-codes)). As with the Regrow branch, the `_table_reduced` variant of supplement 7 drops nearest-single-elevator price data that is not cleared for distribution |
| `join_csb_supplement_9` | `{state}_CSB{years}_supplement_9_table.parquet` | **County-level**, not field-level — see [2.3](#23-supplementary-data-codes) and [4.4](#44-supplementary-farmland-characteristics) |

### 2.3 Supplementary data codes

| Code | Topic |
|---|---|
| 1 | Census tract / county / state boundaries |
| 2 | Elevation & slope |
| 3 | Watershed (HUC8/10/12) |
| 4 | Distance to primary/secondary roads |
| 5 | Neighboring-field land management |
| 6 | Weather (PRISM) |
| 7 | Crop price (elevator & county level) |
| 8 | Soil composition |
| 9 | USDA Census of Agriculture (county-level) |

### 2.4 Dataset schema reference links

> **`[TODO — paste links here]`**
>
> **Regrow Data**
> - Regrow Output Dataset Schema:
> - Regrow Dataset Overview (slide deck):
> - Regrow_DISES_Supplementary_1-9 Codebook:
> - Regrow Raw Data Explainer: 
>
> **CSB Data**
> - CSB Output Dataset Schema:
> - CSB Dataset Overview (slide deck):
> - CSB1724_DISES_Supplementary_1-9 Codebook:
>
> **DISES Data**
> - DISES Codebook:
> - DISES Survey Data Overview (slide deck):
>
> **Soil Supplement (Supplement 8)**
> - Supplementary Information on Soil Dataset Design:

---

## 3. How to Obtain Output Datasets

### 3.1 The Snakefile is the single source of truth for execution order

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

### 3.2 Running a specific rule or chain of rules

Snakemake will automatically build all upstream dependencies of whatever rule you target. To build just one output:

```bash
snakemake <rule_name> --cores <N>
```

For example, to merge Regrow with supplement 1 to get the `{state}_regrow_supplement_1_table` output file:

```bash
snakemake join_regrow_supplement_1 --cores 1
```

### 3.3 Recommended execution chains

**A. Supplementary/reference data (run once, shared by both branches):**
```
download_census_tract → download_county_bound → download_state_bound
download_elevation → reproject_elevation → calculate_slope → clip_elevation_slope
download_watershed
download_roads
download_weather → clip_reproject_weather_rasters
clean_grain_price
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
(`rasterize_regrow_to_gSSURGO_grid` also belongs in this chain if you intend to run supplement 8.)

**D. Regrow + DISES + supplementary joins:**
```
join_regrow_dises
join_regrow_supplement_1 … join_regrow_supplement_6
join_regrow_supplement_7 → cut_regrow_supplement_7
join_regrow_supplement_8  (requires rasterize_regrow_to_gSSURGO_grid first)
join_regrow_supplement_9  (requires join_regrow_supplement_1 first)
```

**E. CSB data preparation:**
```
download_csb
clip_csb_shape
split_csb_shape_table
check_csb_shape_table
rasterize_CSB_to_gSSURGO_grid   (only needed for supplement 8)
```

**F. CSB + DISES + supplementary joins:**
```
join_csb_dises
join_csb_supplement_1 … join_csb_supplement_6
join_csb_supplement_7 → cut_csb_supplement_7
join_csb_supplement_8
join_csb_supplement_9  (requires join_csb_supplement_1 first)
```

**G. Regrow validation against CDL (optional, independent of the merge chains above):**
```
download_cdl
clip_cdl_rasters
rasterize_regrow_to_CDL_grid
join_regrow_cdl
validate_regrow_crop_with_cdl
```

---

## 4. Dataset Documentation (FAQ)

### 4.1 Regrow Dataset

<br>

**What is the Regrow dataset?**

Regrow is a proprietary dataset built from satellite remote sensing (Landsat 5–9 and Sentinel-1/2, with coverage dating back to 1985, though this project uses 2014–2025 land management data). Regrow's **Parcel ID** model delineates agricultural field boundaries directly from satellite imagery — a "Regrow field" is the output of that delineation algorithm, not a property/tax boundary. For each delineated field, Regrow's **MonitorML** practice-detection models report a timeline of land preparation, tillage intensity, cover crop presence, and main commodity crop, for each state in the study region, from 2014–2025.

<br>

**What is a "Regrow field"?**

A Regrow field is a parcel-shaped output of Regrow's in-house field-delineation model (Parcel ID), constructed by analyzing satellite imagery to find areas with homogeneous cropping sequences. It is **not** the same kind of object as a "DISES parcel" (see [4.3](#43-dises-data)), which is derived from tax/ownership records — the two can and do diverge, sometimes substantially, in shape and extent for the same physical land. Background on the delineation methodology: Regrow's [Introduction to Remote Sensing for Agriculture Monitoring](https://help.regrow.ag/introduction-to-remote-sensing-for-agriculture-monitoring) and [A Deep Dive into MonitorML](https://help.regrow.ag/a-deep-dive-into-monitorml).

<br>

**What is the geographic span of the Regrow data?**

Regrow data is processed for 7 states: Illinois, Indiana, Iowa, Michigan, Minnesota, Ohio, Wisconsin. For Michigan specifically, coverage is limited to two counties adjacent to Ohio, rather than the full state.

<br>

**What land management practices are available in the Regrow dataset?**

Four practice categories, organized around the concept of a **cultivation cycle** (see below):

- **Post-harvest tillage** — tillage intensity in the observation window following a harvest
- **Cover crop** — whether a cover crop was present between two commodity crops
- **Pre-plant tillage** — tillage intensity in the observation window preceding main crop planting
- **Main (commodity) crop** — which crop was grown, with planting/harvest dates

<br>

<br>

**How does Regrow detect and define tillage classes?**

Regrow does not detect tillage *events* directly — it estimates **crop residue cover** using indices like the Normalized Difference Tillage Index (NDTI) and Crop Residue Cover Index (CRC), derived from satellite imagery, observed over two 8-week windows per cultivation cycle: one starting at harvest ("post-harvest") and one ending at planting ("pre-plant"). Tillage class is assigned from residue cover following USDA's standard residue-to-tillage-practice thresholds, which differ for fragile vs. non-fragile crops:

| | Fragile crops | Non-fragile crops |
|---|---|---|
| No-till | >30% residue | >60% residue |
| Reduced till | 15–30% residue | 30–60% residue |
| Conventional till | 0–15% residue | 0–30% residue |
| Example crops | Soybean, cotton, wheat, canola, dry bean, sugar beet | Corn, alfalfa, sugarcane, sorghum, switchgrass |

Per Regrow's [Remote Sensing Practice Detection Technology](https://help.regrow.ag/remote-sensing-practice-detection-technology), the actual determination process is:

> "Weekly residue cover percentages are calculated over the two 8-week observation periods, pre-plant and post-harvest. Observations with the highest confidence (where the greatest area of the field was observed and where residue cover estimations within the field were consistent) are identified, and a tillage intensity classification is made based on the median residue percent of those observations. This provides two tillage intensity estimations for each crop cycle."

Additional information on how Regrow defines tillage practices and the peculiarities of tillage detection is available in Regrow's [Tillage Definitions Explained](https://help.regrow.ag/tillage-definitions-explained) guide.

<br>

**How does Regrow detect and define cover crop classes?**

Regrow reports a cover crop **presence indicator**, not a cover crop *type* — fields are evaluated for sustained "green cover" during the window between one commodity crop's harvest and the next planting, requiring at least 8 weeks of observation (to exclude temporary regrowth, weeds, or volunteer plants). Greenness is measured via Normalized Difference Vegetation Index (NDVI), sampled roughly every 5–10 days and averaged monthly, then compared against a **regional NDVI threshold calculated annually** from surrounding natural herbaceous vegetation, to account for local and year-to-year variation in weather (winter conditions can push this threshold toward a seasonal minimum).

![Regrow cover crop NDVI classification thresholds](https://help.regrow.ag/hs-fs/hubfs/Screenshot%202025-04-08%20at%2011-53-18%E2%80%AFAM-png.png?width=933&height=1143&name=Screenshot%202025-04-08%20at%2011-53-18%E2%80%AFAM-png.png)  
Source: https://help.regrow.ag/remote-sensing-practice-detection-technology

There are five cover crop categories in Regrow data: three reflect sustained greenness, ordered from lowest to highest, and two reflect a different kind of situation (a data-availability gap, or a field where the observation doesn't apply at all) rather than a point on the greenness scale:

- **Cover crop not tracked** — average NDVI over the full observation window is low, comparable to having no green cover at all
- **Potential cover crop** — some sustained greenness detected, but not strong/consistent enough to confidently call a cover crop
- **Cover crop** — NDVI sustained above the regional/seasonal threshold throughout the observation window — confident detection of a cover crop
- **No green cover data** — there was not enough high-quality remote-sensing imagery available during the observation window to make a meaningful cover-crop determination (analogous to the "no tillage data" outcome for tillage — a data-availability gap, not a low-NDVI classification)
- **Cover crop not applicable** — the field is not expected to have a cover crop at all in a given year, because a perennial crop is grown on that field that year

One source of ambiguity worth noting: a field where a cover crop was intentionally seeded but failed to establish (e.g., winterkill, or non-sustained emergence) will also show low average NDVI over the observation period, so it can land in either the "Cover crop not tracked" or "Potential cover crop" class — the algorithm cannot distinguish "no cover crop was attempted" from "a cover crop was attempted but didn't survive," since both look like low, non-sustained greenness. 

Additional information on cover crop detection is available via this link: https://help.regrow.ag/remote-sensing-practice-detection-technology

<br>

**How does Regrow detect and define commodity crops?**

Regrow's **CropID** is a machine-learning model trained to detect and identify common crops in a field. It analyzes multi-spectral Sentinel-2 satellite data throughout the growing season to determine crop type, based on historical data and crop-specific spectral signatures (learned reflectance values for each crop). By learning these reflectance values, CropID can predict crop type just one to two months after harvest, and in some cases even during the growing season.

The model also uses satellite imagery and vegetation indices — NDVI and EVI, which measure the chlorophyll content in crops — together with region-specific crop calendars, which account for local farming practices, to identify crop emergence and senescence. The reported "planting date" is the **detected emergence date**, not a seeding date: the actual planting date may differ from the detected emergence, and harvest timing may vary depending on how much green vegetation remains.

Additional information on commodity crop detection is available via this link: https://help.regrow.ag/remote-sensing-practice-detection-technology

<br>

**Transformation of the Regrow commodity crop variable**

`clean_regrow_table.py` recodes Regrow's string crop names (e.g. `"corn"`, `"soybean"`, `"wheat_winter"`) into integer codes via a hardcoded crosswalk (`main_crop_map`). This crosswalk is deliberately aligned with the **USDA CDL/CSB numeric crop code scheme** (e.g. corn=1, soybean=5, winter wheat=24) so that Regrow crop calls can be directly compared against CDL pixel values without further remapping, and so that Regrow and CSB main-crop codes are directly comparable to each other. One sentinel value, `999`, is used for `non_cropland`/`berry` — `999` is **not** a real CDL code; it's a pipeline-internal "doesn't map cleanly" marker.

Performing this transformation also serves a data-quality purpose beyond enabling CDL comparison: building the crosswalk requires explicitly enumerating every distinct crop-name string Regrow uses, which surfaces cases where the same underlying crop is recorded under two different raw names (e.g. `"flax"` and `"flaxseed"`, or `"evergreen"` and `"forest_evergreen"`) so they can be deliberately mapped to the same numeric code rather than being silently treated as different crops downstream.

<br>

**How are raw and processed Regrow datasets structured? What is a cultivation cycle?**

Raw Regrow data arrives as a **long table, one row per cultivation cycle per field**. `clean_regrow_table.py` reorganizes this into the project's final format: a **wide table, one row per field**, with practice/crop values spread across columns indexed by calendar year and cultivation cycle number.

A **cultivation cycle** is the unit that raw Regrow data is built around: it groups the full sequence of land management activities — post-harvest tillage → cover crop → pre-plant tillage → main crop — that contribute to growing **one** commodity crop, following a **harvest-to-harvest model** (not a calendar-year model). A cycle starts the day after the previous crop's harvest and ends at the next commodity's harvest. Because a cultivation cycle is defined around a single commodity crop by construction, **a field cannot have more than one main crop within a single cultivation cycle** — that scenario is excluded by the definition of a cultivation cycle. What *can* happen is multiple complete cultivation cycles occurring within the same calendar year (e.g. two crops grown in succession in one growing season); the pipeline handles this by adding a cultivation-cycle-number suffix to the variable name, rather than by allowing multiple crops inside one cycle. Fallow periods are themselves reported as a "fallow" cultivation cycle, so the data structure stays uniform even when a field sits idle. Source: [Monitor Cultivation Cycles](https://help.regrow.ag/monitor-cultivation-cycles).

Reorganizing the raw long table into the final wide table involves three steps:

1. **Assign each cultivation cycle a calendar year.** The pipeline uses `end_year` — the calendar year in which that cycle's main crop was harvested — as the cycle's year label. `end_year` is used (rather than, say, the cycle's start date) because the harvest date is the one fixed, unambiguous point in a harvest-to-harvest cycle: a cycle's *start* is just whenever the *previous* crop happened to be harvested, which could fall in the prior calendar year, but its *end* is always the harvest of the crop the cycle is actually about.
2. **Number the cultivation cycles within each field-year.** Within a given field and `end_year`, cycles are numbered in chronological order as `cycle_count` (1st, 2nd, or 3rd cycle ending in that calendar year — capped at 3). This `cycle_count` is the mechanism that captures multiple cultivation cycles landing in the same calendar year, since (per above) it's never the case that one cycle itself contains more than one main crop.
3. **Pivot from long to wide table.** The long table is pivoted so that each `(activity, end_year, cycle_count)` combination becomes its own column, producing the variable-naming convention used throughout this project:

   ```
   {land management activity}_{year}_{cultivation cycle number}
   ```
   where `land management activity ∈ {PHtill, cover, PPtill, crop}`, `year` is the 2-digit `end_year`, and `cultivation cycle number ∈ {1, 2, 3}`.

**Handling multiple records within the same cycle.** Although a cycle can only have one main crop, it can have more than one *recorded tillage event* of the same type (e.g. a field tilled twice within one cycle's post-harvest window). Before pivoting, such repeated records are aggregated using these rules: if the repeated records share identical dates and identical values, the value is kept as a single observation; if dates are identical but the recorded values differ, the result is set to missing (the conflicting records can't both be right); if the dates genuinely differ (i.e. these are two distinct tillage events, not duplicate records of the same event), the two values are joined into one **compound code** representing the pair — for example, a field tilled twice within one cycle, once no-till and once reduced-till, is recorded as a single combined tillage code (5) representing that specific pair, rather than as two separate tillage columns. The full set of these paired codes is in the table below.

Because more than one tillage *event* can occur within a single cycle, the tillage code is not a strict 1/2/3 ordinal scale — it also encodes same-or-mixed-intensity *pairs*:

| Code | Meaning |
|---|---|
| 1 | No-till (single event, or two no-till events) |
| 2 | Reduced till (single event, or no-till + reduced-till pair) |
| 3 | Conventional till (single event, or no-till + conventional-till pair) |
| 4 | Reduced-till + reduced-till pair |
| 5 | Reduced-till + conventional-till pair |
| 6 | Conventional-till + conventional-till pair |

**Cover crop codes**: `1` = no cover crop, `2` = potential cover crop, `3` = detected cover crop, `999` = cover crop not applicable.

<br>

**What is the Regrow confidence score?**

Confidence quantifies how reliable a given practice classification is, and is affected by remote-sensing data availability/quality (cloud cover degrades it). For **main crop**, confidence is on a 0–100 scale and is explicitly probabilistic: "a confidence score of 80 means roughly 8 out of 10 times the model gets this crop type right" — Regrow recommends treating ≥75 as reliable. For **tillage** and **cover crop**, confidence is on a coarser 1–3 scale; Regrow recommends relying on observations with confidence = 3. Source: Monitor Data Explainer files, [Remote Sensing Practice Detection Technology](https://help.regrow.ag/remote-sensing-practice-detection-technology).

<br>

**How is confidence score determined for tillage, cover crop, and main crop?**

- **Tillage**: confidence is rated 1–3 for each weekly residue observation in the 8-week window. A given week's confidence reflects how much of the field was actually observed (cloud cover reduces this) and how spatially consistent the residue-cover estimates were across the observed area (low variability across the field = higher confidence); confidence is also reduced when green crop canopy covers the soil, since residue is harder to detect underneath it. Per Regrow's documentation, the observations with the highest confidence in the 8-week window are identified, and the tillage classification is based on the median residue percent of those highest-confidence observations — producing two tillage estimations per cultivation cycle (one pre-plant, one post-harvest).
- **Cover crop**: confidence (0–3) is the average percent of the observation window actually observed via satellite (more cloud-free, high-quality observations → higher confidence).
- **Main crop**: confidence (0–100) reflects the statistical reliability of the crop-type call itself, from the CropID model.

<br>

**How does Regrow address uncertainty from remote sensing data and algorithm errors?**

Cloud cover and sensor limitations are the primary sources of observational error; Regrow mitigates this by using multi-week observation windows rather than single snapshots and explicit confidence scoring at every level (tillage, cover crop, crop). It also creates **intermediate classes** for cases of genuine ambiguity rather than forcing a binary call — e.g. "potential cover crop" sits between "not tracked" and "cover crop detected" specifically to flag low-confidence green-cover signals.

<br>

**How is Regrow data validated?**

Regrow's own validation (per their public documentation) compares crop-type predictions against ground-truth survey/government data using various quality metrics on a representative benchmark subset, with reported accuracy ranging ~70–95% depending on crop and region; tillage detection has been benchmarked at ~71% accuracy distinguishing conventional vs. conservation tillage over 22,000+ observations, and cover crop detection at ~85% accuracy over 25,000+ field observations (continental US).

**In addition, this project runs its own independent validation** of Regrow's main-crop calls against USDA's Cropland Data Layer (CDL), implemented in `rasterize_regrow_to_CDL_grid.py` → `join_regrow_cdl.py` → `validate_regrow_crop_with_cdl.py`:

1. Regrow field polygons are rasterized onto the exact CDL pixel grid (pixel-center rule), after dropping near-duplicate overlapping polygons and resolving the grid alignment across all CDL years.
2. For each state/year, the **modal** CDL crop class observed under each field's footprint is computed.
3. A field is flagged `cdl_valid_{year} = 1` if that modal CDL class matches *any* of the field's Regrow-reported main crop codes for that year (across all recorded cultivation cycles), else `0`.
4. `validate_regrow_crop_with_cdl.py` summarizes percent-valid fields by year and by crop category (all crops; corn only; soybean only; winter wheat only; corn+soy+all-wheat combined), saved as one Excel file per state.

<br>

**What are the key attributes of the Regrow dataset?**

Regrow's key attributes are: field size (area), the land management practice variables themselves — tillage intensity, cover crop presence, and main crop — and the attributes associated with each of those land management values: their start/end dates and confidence scores. Alongside these are the field's unique identifier (`field_id`) and geographic location (geometry).

<br>

**What naming conventions are used for Regrow variables?**

`{activity}_{year}_{cultivation cycle number}` for practice/crop values (e.g. `crop_23_1`), with `_start`, `_end`, `_conf` suffixes for the corresponding date/confidence companions of each (e.g. `crop_start_23_1`, `crop_conf_23_1`).

<br>

**Where can I find more information about Regrow dataset attributes?**

`[TODO — paste links]`
- Regrow Codebook:
- Regrow Dataset Overview (slide deck):

---

### 4.2 USDA Crop Sequence Boundaries (CSB) Dataset

<br>

**What is the CSB dataset?**

CSB is a **public** USDA dataset (developed with USDA's Economic Research Service) that produces field boundary, crop acreage, and crop-rotation estimates across the contiguous United States. CSB boundaries are derived algorithmically from the **Cropland Data Layer (CDL)** — an open-source method (available on GitHub) that analyzes multi-year crop-rotation consistency in CDL pixels to infer field boundaries. CSB releases are published annually, each covering a rolling 8-year window; this project uses the 2017–2024 release. Because CSB boundaries are derived purely from satellite-based crop classification, they contain **no personally identifying information** and represent "crops grown only," not ownership or tax-parcel boundaries.

<br>

**The definition of a "CSB field"**

A CSB field (`CSBID`) is a polygon produced by the CSB algorithm's crop-rotation-consistency analysis over CDL pixels for a given multi-year window — conceptually similar to a Regrow field (both are remote-sensing-derived, not ownership-derived) but built with a different methodology and, in practice, different boundaries on the ground for the same physical land.

<br>

**What is the geographic span of CSB data?**

CSB itself covers the contiguous United States; within this project, CSB data is processed for the same 7 states as Regrow: Illinois, Indiana, Iowa, Michigan, Minnesota, Ohio, Wisconsin, with data available for all Michigan counties.

<br>

**What data is available within the CSB dataset?**

Each CSB field has one main (commodity) crop classification per calendar year, for every year in the release window (`CDL2017`…`CDL2024`), plus the field's acreage (`CSBACRES`) and geographic location.

<br>

**Does the CSB dataset have data on tillage intensity and cover crops?**

No. CSB contains **main crop rotation only** — no tillage intensity, no cover crop indicator, and no within-year practice timeline. This is the single largest informational gap relative to Regrow.

<br>

**What land management practices are available in the CSB dataset?**

Only one: **main (commodity) crop**, classified once per calendar year per field (see "How does CSB classify main crop?" below). There is no post-harvest tillage, cover crop, or pre-plant tillage equivalent.

<br>

**How are raw and processed CSB datasets structured?**

Raw CSB data is distributed as a single national geodatabase per release window (e.g. `CSB1724.gdb`), containing one polygon per field with one `CDL{year}` column per year in the window. USDA publishes CSB in rolling 8-year windows (e.g. 2016–2023, 2017–2024, 2018–2025), each released as a new geodatabase; this pipeline is currently configured for the 2017–2024 window (`CSB1724`, set via `CSB_years` in `config.yml`). `clip_csb_shape.py` extracts each state's portion via an attribute filter on `STATEFIPS`, and `split_csb_shape_table.py` splits each state's file into a geometry-only file (`CSBID` + polygon), an attributes-only file (`CSBID`, `CSBACRES`, and the `CDL{year}` columns), and a combined shape+table file. Unlike Regrow, there is no long-to-wide reshaping step, since CSB's raw structure (one row per field, one column per year) is already in the project's target wide format.

<br>

**How does CSB classify main crop?**

CSB assigns one annual crop code per field per year (`CDL2017`…`CDL2024`), taken as the majority CDL pixel class within that field's polygon footprint for that year — i.e., one classification per calendar year, not per cultivation cycle. The crop-code scheme is the standard USDA CDL numeric classification (0=background, 1–60 major row/specialty crops, 61–65 non-crop/fallow/pasture/forest/shrub, 66–80 tree crops, 81–109 other land cover, 110–195 NLCD-derived classes, 204–254 specialty crops) — the same scheme Regrow's crosswalk targets, which is what makes Regrow and CSB main-crop codes directly comparable.

<br>

**How is CSB data validated by USDA?**

`[TODO: ADD THIS INFORMATION]`

<br>

**What are the key differences between Regrow and CSB?**

| | Regrow | CSB |
|---|---|---|
| License | Proprietary | Public |
| Tillage | Yes (3 classes, two windows/year) | No |
| Cover crop | Yes (presence indicator, 3 classes) | No |
| Main crop | Yes, per cultivation cycle (can be multiple/year) | Yes, one classification per calendar year |
| Timeline detail | Planting/harvest dates per cycle | Annual only |
| Year coverage (this project) | 2014–2025 | 2017–2024 |
| Michigan coverage | Limited to 2 counties adjacent to Ohio | All counties |

Both datasets are conceptually similar (remote-sensing-derived field boundaries with crop classification) but use different methodologies, so the same physical field will generally have different boundaries and possibly different crop calls in the two datasets.

<br>

**What are the key attributes of the CSB dataset?**

Field size (`CSBACRES`) and the main crop classification for each year in the release window (`CDL{year}`), alongside the field's unique identifier (`CSBID`) and geographic location (geometry).

<br>

**What naming conventions are used for CSB variables?**

The outcome variable is `CDL{year}` for `year ∈ {2017...2024}`.

<br>

**Where can I find more information about CSB dataset attributes?**

`[TODO — paste links]`
- Crop Sequence Boundaries Metadata:
- CSB1724_DISES_Supplementary_1-8 Codebook:

---

### 4.3 DISES Data

<br>

**What is the DISES dataset?**

DISES (Dynamics of Integrated Socio-Environmental Systems) aggregates farmer survey responses collected in early 2024 about the **2023 growing season**. Respondents were onboarded via stratified random sampling of corn-soy owner-operators in the western Lake Erie Basin (Maumee watershed), across different farm sizes, restricted to farmers with at least one owned field (or contiguous owned fields) so their land could be spatially identified. **645 respondents** participated.

DISES is currently the project's only source of direct farmer characteristics. At the initial stage of our project, the team is using a smaller subset of Regrow/CSB fields with farmer-characteristic data to evaluate how predictive those characteristics are of land management decisions; results will inform whether/how to expand future survey rounds.

<br>

**Which area is covered by DISES?**

DISES survey responses come specifically from the Maumee watershed / western Lake Erie Basin portion of Indiana, Michigan, and Ohio (`states_DISES` in `config.yml`) — not the entirety of those three states. This is a sub-region of the broader 7-state Regrow/CSB study area.

<br>

**What is a "DISES parcel"?**

A DISES parcel is a tax-record-derived land unit. During dataset processing, all individual tax parcels belonging to the same survey respondent (`comp_id`) are **dissolved into one combined multi-part geometry** representing that farmer's entire landholding footprint — so the unit that actually gets spatially joined against Regrow/CSB fields is "one farmer's combined holdings," not "one tax parcel."

<br>

**What is the difference between a "Regrow field" and a "DISES parcel"?**

A Regrow field is the output of a remote-sensing field-delineation *algorithm* based on observed crop-choice sequences; a DISES parcel is based on **tax/ownership records**, collected for a separate project with no direct association to this one. For the same physical land, the two may disagree in: **shape** (the polygon outline itself), **extent** (how much area is covered — a DISES holding may be larger or smaller than the corresponding Regrow fields), and **count** (one farmer's DISES holding may correspond to one Regrow field, several Regrow fields, or only partially overlap any Regrow field at all). This is exactly why the merging procedure in [4.5](#45-merging-regrow-dises-and-supplementary-data) needs explicit match-quality scoring rather than assuming a clean 1:1 correspondence. Because DISES was built independently, users should expect — and the match-quality variables are designed to surface — inconsistencies between the two datasets.

<br>

**What are the key attributes of the DISES dataset?**

`[TODO — ADD INFORMATION]`

<br>

**Where can I find more information about DISES dataset attributes?**

`[TODO — paste links]`
- DISES Codebook:
- DISES Survey Data Overview (slide deck):

---

### 4.4 Supplementary Farmland Characteristics

<br>

**Why do we need additional farmland characteristics?**

Land management decisions are plausibly influenced by physical, economic, and locational context beyond what Regrow/CSB/DISES capture directly — soil quality, terrain, weather, market access, nearby farming behavior, and crop prices. These 9 supplementary blocks attach that context to every Regrow/CSB field so it can be used as covariates in downstream analysis.

<br>

**What supplementary farmland characteristics are available?**

| # | Topic | Source | Static or dynamic | Time coverage | Resolution |
|---|---|---|---|---|---|
| 1 | Census tract/county/state boundaries | TIGER 2023 | Static | 2023 release | Tract-level polygons |
| 2 | Elevation & slope | USGS 3DEP | Static | Current 3DEP survey | 1 arc-second (~30m) |
| 3 | Watershed (HUC8/10/12) | USGS 3DHP / NHD | Static | Current NHD/3DHP delineation | Polygon (subbasin/watershed/subwatershed) |
| 4 | Roads (primary/secondary) | TIGER PRISECROADS 2023 | Static | 2023 release | Vector lines |
| 5 | Neighboring-field land management | Derived from Regrow/CSB itself | Dynamic | Same years as host dataset (Regrow: 2014–2025, CSB: 2017–2024) | Field-level, derived from the host dataset itself |
| 6 | Weather (PRISM) | PRISM AN Monthly | Dynamic | Monthly, 2014–2025 | 800m |
| 7 | Crop price (elevator & county) | Barchart (proprietary) | Dynamic | Monthly; price series with fewer than a configured minimum number of observations between 2015 and 2025 are dropped | Point-level (elevator) / county-level (index), spatially interpolated to fields |
| 8 | Soil composition | gSSURGO | Static | Current gSSURGO survey release | 30m (map-unit raster) |
| 9 | USDA Census of Agriculture (county-level) | USDA NASS QuickStats | Static | 2017 & 2022 Census of Agriculture | County-level (not field-level) |

"Static" means the dataset reflects conditions at a single point in time (e.g. one published survey/release), rather than tracking change over time. "2023 release" means the specific TIGER vintage published by Census in 2023 — i.e., the dataset as it existed at that one publication date, not a 2023 measurement of a changing quantity.

**What datasets are used to generate each supplementary block? — detail**

- **Supplement 1 (Census tract)**: TIGER 2023 Cartographic Boundary Files. Fields are matched to their containing census tract by a centroid-in-polygon spatial join (point-in-polygon, faster than full overlay); fields whose centroid falls outside all tract polygons due to boundary precision issues fall back to a nearest-tract join. Adds state/county/tract identifiers and tract land/water area.
- **Supplement 2 (Elevation & slope)**: USGS 3DEP 1-arc-second elevation, reprojected to the project's equal-area CRS (EPSG:5070) and used to derive slope (degrees) via GDAL's DEM processing. Field-level values are the **zonal mean** elevation/slope under each field's polygon footprint.
- **Supplement 3 (Watershed)**: USGS 3DHP/NHD hydrologic unit boundaries at three nested levels — subbasin (HUC8), watershed (HUC10), subwatershed (HUC12) — matched to fields via the same centroid-in-polygon approach as supplement 1.
- **Supplement 4 (Roads)**: TIGER 2023 PRISECROADS (primary + secondary roads). Computes nearest-road distance per field via nearest-neighbor spatial join, with primary roads (MTFCC class S1100) preferred over secondary (S1200) when distances tie.
- **Supplement 5 (Neighboring-field land management)**: Derived entirely from the host dataset (Regrow or CSB) itself — no external source. Fields sharing a boundary or overlapping with a given field are identified via spatial self-join; tillage/cover-crop values (Regrow only) are averaged across neighbors, and main-crop composition is summarized as an area-weighted share of neighboring acreage across crop groups (corn / soybean / wheat / other). The CSB-side version computes only the crop-composition share, since CSB has no tillage/cover-crop fields to average.
- **Supplement 6 (Weather)**: PRISM AN Monthly Time Series, 6 variables (precipitation, min/max temperature, mean dew point, min/max vapor pressure deficit), 2014–2025, 800m resolution. Field-level values are extracted as a **point sample at the field centroid** (not a full zonal average) for performance — column names use `_centroid_` to reflect this.
- **Supplement 7 (Crop price)**: Barchart elevator-level monthly cash prices for corn/soybeans/wheat, and county-level monthly cash prices for corn/soybeans (no county series for wheat). Field-level values are computed via a K-D tree nearest-neighbor search against geocoded elevator/county locations, both as a "nearest with gap-filling from the next-nearest" series and as a distance-weighted average across the `N` nearest locations (own-county is always prioritized first in the county-level search, when available). The single-nearest-elevator series is dropped in a post-processing step (`cut_*_supplement_7`), keeping only the N-nearest weighted average and the county-level series in the final `_table_reduced` file.
- **Supplement 8 (Soil/gSSURGO)**: Gridded SSURGO database (USDA NRCS), 30m map-unit (mukey) raster, joined to detailed component/horizon tabular soil data. For each field, the dominant soil component's categorical attributes (`drainagecl`, `cropprodindex`, `resdept_r`) are taken from the single most areally-significant soil component per map unit; continuous physical/chemical properties (`sandtotal_r`, `claytotal_r`, `ph1to1h2o_r`) are depth-weighted to a configurable target depth (`soil_depth_cm`, currently 30cm) and then averaged across all map units intersecting the field, weighted by pixel coverage (i.e., area-weighted); `slopegraddcp` is taken directly as a pre-aggregated mukey-level attribute.
- **Supplement 9 (USDA Census of Agriculture)**: USDA NASS QuickStats county-level data from the 2017 and 2022 Census of Agriculture — farmland/cropland/pastureland/woodland acreage and ownership tenure, harvested acreage for corn/soybeans/winter wheat, conservation tillage and cover crop acreage, and producer demographic/decision-making characteristics. **Unlike supplements 1–8, this is a county-level dataset, not a field-level one** — `join_regrow_supplement_9.py` / `join_csb_supplement_9.py` do not attach these values to individual `field_id`/`CSBID` records. Instead, each script reads the corresponding `supplement_1` table to get the `county_name`/`state_name` already assigned to each Regrow field or CSB field, builds the distinct, alphabetically-sorted list of `county_state_name` values (`county_name + "_" + state_name`) actually present in that state, and joins the wide Census of Agriculture table onto that county list. `county_state_name` is therefore the matching key to use when joining supplement 9 to any other data block that has a county/state identifier.

<br>

**Where can I find additional information about supplementary data attributes?**

`[TODO — paste Supplementary Information on Soil Dataset Design link, and any other per-supplement schema documentation]`

---

### 4.5 Merging Regrow, DISES, and Supplementary Data

#### 4.5.1 How are Regrow and DISES spatially merged?

Implemented in `join_regrow_dises.py`. The merge is built **around the Regrow field** — Regrow field is the primary unit of the output dataset, and DISES attributes are attached to it (not the other way around):

1. Each Regrow field polygon is slightly buffered (a small negative margin, `buffer_margin` in config) and intersected against all DISES parcel-holding polygons.
2. Where a Regrow field overlaps multiple DISES holdings, only the **largest-overlap** DISES match is kept per Regrow field.
3. Overlap area (in acres) and overlap share (as a percent of the Regrow field's own area) are computed from the *original* (unbuffered) geometries.
4. A `field_assigned_dises` flag (`Y`/`N`) marks whether any DISES match was found at all.

#### 4.5.2 What are the outcomes of the matching?

`[TODO — fill in the counts below]`

| State | Regrow fields with a DISES match (`field_assigned_dises == "Y"`) | Regrow fields with a DISES survey (`survey_responded_dises == "Y"`) | Total Regrow fields
|---|---|---|---|
| IN | | | |
| MI | | | |
| OH | | | |
| **Total** | | | |

Unique DISES holdings with at least one overlapping Regrow field: `[ ]` out of `[ ]` total DISES holdings (`[ ]`%).

#### 4.5.3 How is matching quality evaluated?

Once a Regrow field has a DISES match, `join_regrow_dises.py` scores match quality using two independently-available attributes: **field area** and **2023 main crop**.

- `area_match_dises` = 1 if the Regrow field's area falls within `[lower, upper] × DISES-reported field size` (currently 0.75×–1.25×, `area_match_coefs` in config), else 0.
- `crop_match_dises` = 1 if the Regrow field's 2023 main crop (either of its first two cultivation cycles) matches the DISES-reported 2023 crop, else 0.
- These combine into `match_quality_dises`:

| `match_quality_dises` | Meaning |
|---|---|
| `A` | Both area and crop match |
| `B_crop` | Crop matches, area does not (or is unavailable) |
| `B_area` | Area matches, crop does not (or is unavailable) |
| `F` | Neither matches |
| *(missing)* | Neither comparison was possible (both inputs missing) |

This is also where the **representative field** indicator is created (see next question).

#### 4.5.4 What is a "representative field" attribute? Why do we need it?

The DISES survey asks farmers a block of questions about their "representative field," without specifying *which* field they meant. Since DISES holdings are dissolved per-farmer (see [4.3](#43-dises-data)) and can correspond to multiple candidate Regrow fields, this attribute identifies the most plausible candidate by comparing each candidate Regrow field's area and 2023 main crop against the farmer's representative-field survey answers.

The procedure (from `Representative Field Attribute. Supporting Documentation.md`, implemented in `assign_representative_field()` function inside `join_regrow_dises.py`):

1. **Eligibility filter**: only Regrow fields that are DISES-assigned, belong to a farmer who completed the survey, and have `overlap_area_share_dises > 50%` are considered candidates at all — this conservative threshold excludes ambiguous cases like a DISES holding entirely containing a much larger Regrow field, or a Regrow field with only a sliver of overlap with a DISES holding.
2. **Tiered classification** into 4 levels of likelihood (combining `match_quality_dises` with the Regrow 2023 crop-confidence score, thresholded at 75%):

   - **Level 1** (highest likelihood): `match_quality_dises = "A"` and crop confidence > 75%.
   - **Level 2** (medium): `"A"` with confidence ≤ 75%; or `"B_area"` with the DISES crop answer missing; or `"B_crop"` with the DISES size answer missing and confidence > 75%.
   - **Level 3** (low-medium): `"B_crop"` with size missing and confidence ≤ 75%; or `"B_area"` with a crop mismatch and confidence ≤ 75%.
   - **Level 4** (low): `"F"` with size missing, a crop mismatch, and confidence ≤ 75%; or both size and crop answers missing but the farmer has only a single DISES holding (`n_parcels_dises = 1`), making that the only possible candidate by elimination.
3. **One representative field per farmer**: among all of a farmer's eligible candidates, the highest-likelihood-level field is kept; ties are broken by closest area match, then by largest overlap share.

To get the full sample of representative fields, combine fields across all 4 levels (`RF_level_1_dises` … `RF_level_4_dises`).

> **Note on the CSB-side equivalent (`join_csb_dises.py`)**: the same tiered logic applies, but **without a confidence-score dimension**, since CDL/CSB has no per-field confidence analog to Regrow's ML crop confidence. CSB's `RF_level_1` is therefore simply `match_quality_dises == "A"` with no confidence split, and `RF_level_2`/`RF_level_3` collapse the Regrow version's confidence-based sub-branches into their non-confidence conditions only. This is an intentional simplification driven by data availability, not an oversight.

#### 4.5.5 What are the outcomes of the matching quality procedure?

`[TODO — fill in the counts below]`

| `match_quality_dises` | Count |
|---|---|
| `A` | |
| `B_crop` | |
| `B_area` | |
| `F` | |
| **Total matched** | |

| Representative field level | Count |
|---|---|
| `RF_level_1_dises` | |
| `RF_level_2_dises` | |
| `RF_level_3_dises` | |
| `RF_level_4_dises` | |
| **Total representative fields** | |

#### 4.5.6 How to generate Regrow sub-samples with available DISES data

| Sub-sample | Filter |
|---|---|
| Regrow fields with an assigned DISES field (survey or not) | `field_assigned_dises == "Y"` |
| Regrow fields with an assigned DISES field **and** survey data | `field_assigned_dises == "Y"` and `survey_responded_dises == "Y"` |
| Regrow fields matching the assigned DISES field by size or crop | the above, **and** `match_quality_dises` is one of `"A"`, `"B_area"`, `"B_crop"` |

#### 4.5.7 How to extract DISES-related variables

Every column originating from or derived from the DISES dataset carries the suffix `_dises` (e.g. `field_crop_23_dises`, `match_quality_dises`, `RF_level_1_dises`).

#### 4.5.8 Merging Regrow with each supplementary dataset

Each supplementary join is implemented as its own script and produces its own standalone output table (joined to the host dataset only by `field_id`, not chained to each other) — see [4.4](#44-supplementary-farmland-characteristics) for the per-supplement methodology, and [Section 6](#6-explanation-of-individual-scripts) for full script-level detail. The one exception is **supplement 9**, which is keyed on `county_state_name` rather than `field_id` — see [4.4](#44-supplementary-farmland-characteristics) for how it's built and matched.

#### 4.5.9 Output dataset structure and integration

Output files generally follow the pattern `{state}_regrow_<block>_table.parquet` — one file per state per data block, e.g. `{state}_regrow_table.parquet` for the base table, `{state}_regrow_dises_table.parquet` for the DISES join, and `{state}_regrow_supplement_{n}_table.parquet` for supplement `n`.

To combine multiple blocks into one field-level dataset, merge them on `field_id` — every block's table has one row per Regrow field (fields with no DISES or supplement match still have a row, just with missing values), so `field_id` is always the key to use when joining any two of these tables together.

---

### 4.6 Merging CSB, DISES, and Supplementary Data

#### 4.6.1 How are CSB and DISES spatially merged?

Identical methodology to [4.5.1](#451-how-are-regrow-and-dises-spatially-merged), implemented in `join_csb_dises.py`, with `CSBID` as the join key and the merge built around the **CSB field** as the primary output unit. One implementation difference: CSB processing loops over `CSB_years` (currently just `"1724"`) as an outer dimension, since Regrow has no year-window concept.

#### 4.6.2 What are the outcomes of the matching?

`[TODO — fill in the counts below]`

| State | CSB fields with a DISES match (`field_assigned_dises == "Y"`) | CSB fields with a DISES survey (`survey_responded_dises == "Y"`) | Total CSB fields |
|---|---|---|---|
| IN | | | |
| MI | | | |
| OH | | | |
| **Total** | | | |

Unique DISES holdings with at least one overlapping CSB field: `[ ]` out of `[ ]` total DISES holdings (`[ ]`%).

#### 4.6.3 How is matching quality evaluated?

Same `area_match_dises` / `crop_match_dises` → `match_quality_dises` (`A`/`B_crop`/`B_area`/`F`) logic as the Regrow side, with one necessary difference: CSB's crop match compares the DISES 2023 crop answer against the single `CDL2023` column (CSB has one annual crop call, not multiple within-year cultivation cycles like Regrow).

#### 4.6.4 What is a representative-field attribute on the CSB side?

Same purpose and tiering as [4.5.4](#454-what-is-a-representative-field-attribute-why-do-we-need-it), implemented in `join_csb_dises.py`'s `assign_representative_field()`. As noted above, the CSB version has no crop-confidence dimension (CDL has no per-field ML confidence score), so its tiers are coarser:

- **Level 1**: `match_quality_dises == "A"`.
- **Level 2**: `"B_area"` with DISES crop answer missing, or `"B_crop"` with DISES size answer missing.
- **Level 3**: `"B_area"` with a crop mismatch.
- **Level 4**: `"F"` with size missing and a crop mismatch; or both size and crop missing with only a single DISES holding.

#### 4.6.5 What are the outcomes of the matching quality procedure?

`[TODO — fill in the counts below]`

| `match_quality_dises` | Count |
|---|---|
| `A` | |
| `B_crop` | |
| `B_area` | |
| `F` | |
| **Total matched** | |

| Representative field level | Count |
|---|---|
| `RF_level_1_dises` | |
| `RF_level_2_dises` | |
| `RF_level_3_dises` | |
| `RF_level_4_dises` | |
| **Total representative fields** | |

#### 4.6.6 How to generate CSB sub-samples with available DISES data

Same filter logic as [4.5.6](#456-how-to-generate-regrow-sub-samples-with-available-dises-data), substituting `CSBID` for `field_id`.

#### 4.6.7 How to extract DISES-related variables

Same `_dises` suffix convention as the Regrow branch.

#### 4.6.8 Merging CSB with each supplementary dataset

Structurally near-identical to the Regrow-side scripts (same spatial-join methodology per supplement), with `CSBID` substituted for `field_id` and an added outer loop over `CSB_years`. Note that **Supplement 9** is structurally different from the rest: it's keyed on `county_state_name` rather than `CSBID` — see [4.4](#44-supplementary-farmland-characteristics).

#### 4.6.9 Output dataset structure and integration

Output files generally follow the pattern `{state}_CSB{years}_<block>_table.parquet` — e.g. `{state}_CSB{years}_table.parquet` for the base table, `{state}_CSB{years}_dises_table.parquet` for the DISES join, and `{state}_CSB{years}_supplement_{n}_table.parquet` for supplement `n`.

To combine multiple blocks into one field-level dataset, merge them on `CSBID` — every block's table has one row per CSB field (fields with no DISES or supplement match still have a row, just with missing values), so `CSBID` is always the key to use when joining any two of these tables together.

---

## 5. Documentation

`[TODO — you indicated you'll provide instructions for this section later.]`

---

## 6. Explanation of Individual Scripts

This section catalogs every active script in `scripts/` (excluding `scripts/archive/`, which holds superseded versions kept for reference only). For each: what it produces and why, the key conceptual steps, and notable assumptions. Organized by pipeline stage to match [Section 3](#3-how-to-obtain-output-datasets).

> Several scripts listed here had known bugs identified during this documentation pass; almost all have since been corrected in the current codebase (see [Known Issues](#8-known-issues) for the one that is still open). The descriptions below reflect current, corrected behavior.

### 6.1 Download scripts

**`download_census_tract.py` / `download_state_bound.py` / `download_county_bound.py`**
- *Overview*: Download single Census TIGER 2023 boundary zip files (tract, state, county respectively) and extract them.
- *Key steps*: build the download URL from config, stream-download with a file-size tolerance check, extract the zip.
- *Assumptions*: SSL verification is deliberately disabled for census.gov downloads only. Header access for `Content-Length` now safely defaults to `0` if absent (see [Known Issues](#8-known-issues) for the remaining division-by-zero edge case on the size-tolerance check in these three scripts specifically).

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

### 6.2 Geo-processing: clip / reproject / clean

**`reproject_elevation.py`**: Reprojects the merged national elevation mosaic to the project's equal-area CRS (EPSG:5070) via GDAL Warp.

**`calculate_slope.py`**: Computes slope in degrees from the reprojected elevation raster via GDAL DEM processing. Slope is deliberately computed *after* reprojection to an equal-area CRS, so the gradient calculation uses planar (not geographic lat/lon) units, avoiding the distortion that would result from differentiating directly over a geographic grid.

**`clip_elevation_slope.py`**: Clips the national reprojected elevation and slope rasters to each state boundary, with an outward buffer (`state_buffer_margin` in config) to retain context just past the state line.

**`clip_gSSURGO_mukey_rasters.py`**: Clips the national gSSURGO 30m mukey raster to each state boundary (same buffer parameter) to produce the per-state soil map-unit grids used later in supplement 8.

**`clip_cdl_rasters.py`**: Clips national CDL rasters per state per year (same buffer parameter), producing the `{state}_{year}_30m_cdls_clipped.tif` files used for Regrow's CDL validation. Missing-file years are skipped silently (by design, since not every year may be downloaded yet); any other processing error is re-raised rather than silently skipped.

**`clean_grain_price.py`**: Cleans raw Barchart Excel elevator/county price spreadsheets into monthly average price series and geocodes elevator/county-index locations (via OpenStreetMap Nominatim, with a manually-curated correction file for known bad addresses) for use in supplement 7's nearest-neighbor price join. Drops price series (elevators or county indices) with fewer than a configured minimum number of non-null monthly observations.

**`clip_reproject_weather_rasters.py`**: Reprojects and clips monthly PRISM rasters to each state boundary, producing per-state/variable/month clipped GeoTIFFs used by supplement 6.

**`clip_csb_shape.py`**: Extracts each state's portion of the national CSB geodatabase via an attribute filter on `STATEFIPS`, saving one parquet per state per CSB year-window.

**`split_csb_shape_table.py`**: Splits each state's clipped CSB parquet into geometry-only, attributes-only, and combined shape+table outputs, reprojected to the target CRS.

**`check_csb_shape_table.py` / `check_regrow_shape_table.py`**: Detects pairs of heavily-overlapping (near-duplicate) field polygons within the same dataset via a self-overlay, flags the smaller of each overlapping pair, and saves a small lookup table. This matters specifically for rasterization: if two polygons cover nearly the same area, naive pixel-by-pixel rasterization would arbitrarily assign each pixel to whichever polygon happens to be processed last, producing inconsistent field assignment. The downstream rasterization scripts use this lookup to drop the "smaller" duplicate before rasterizing, then backfill its attributes from the "larger" field's results afterward.

### 6.3 Regrow data preparation chain

**`join_regrow_2025_updates.py`**: Merges legacy (2014–2024) and newly-released (2025) raw Regrow exports before any other processing — combining geometry files (preferring the legacy boundary where both exist) and concatenating the tabular Monitor data.

**`split_save_regrow_geometry.py`**: Produces the canonical `{state}_regrow_fieldID_geometry.parquet`/`.gpkg` files: reprojects to the target CRS, repairs any invalid geometries, and renames the raw `boundary_id` key to the project's canonical `field_id`.

**`split_regrow_monitor_by_state.py`**: Splits the legacy era's multi-state-concatenated Monitor tables back out into true per-state files, by matching each state's own set of field IDs.

**`clean_regrow_table.py`**: The core Regrow cleaning/reshaping script — see the full walkthrough in [4.1](#41-regrow-dataset) ("How are raw and processed Regrow datasets structured?" and "Transformation of the Regrow commodity crop variable"). In brief: aggregates multiple land-management records within a cultivation cycle, assigns each cycle a calendar year and cycle count (capped at 3/year), pivots from long to the final wide `{activity}_{year}_{cycle}` format, and recodes crop/tillage/cover-crop strings to the project's numeric coding schemes. Also unifies the column schema across all states at the end, so every state's output file has an identical (alphabetically sorted) set of columns even if, e.g., one state never recorded a 3rd cultivation cycle in any year.

**`join_regrow_shape_table.py`**: Joins the canonical field geometry with the cleaned wide attribute table on `field_id`, and computes `area_acre` from the polygon area.

### 6.4 CDL validation chain

**`rasterize_regrow_to_CDL_grid.py`**: Rasterizes Regrow field polygons onto the exact CDL raster grid (after dropping near-duplicate overlapping fields and verifying every configured CDL year shares an identical grid), assigning each field a simple integer ID (`pid`) burned into the raster by pixel-center containment.

**`join_regrow_cdl.py`**: For each state/year, extracts the modal CDL crop class under each field's rasterized footprint, backfills near-duplicate fields from their paired "larger" field, and flags `cdl_valid_{year}` based on whether the modal CDL class matches any of that field's Regrow-reported crop codes for that year (excluding the `999` non-cropland sentinel from the comparison).

**`validate_regrow_crop_with_cdl.py`**: Produces the final by-state, by-year, by-crop-category accuracy summary Excel files described in [4.1](#41-regrow-dataset).

### 6.5 DISES preparation chain

**`clean_dises_table.py`**: Cleans the raw combined DISES survey CSV, renames key survey fields, and computes three farmer-typology indices (productivism, conservationism, civic engagement) as row-means of specific Likert-scale survey items.

**`clean_dises_shape.py`**: Cleans the raw DISES parcel shapefile, counts the number of distinct tax parcels per farmer (`n_parcels`), then **dissolves** all of a farmer's parcels into one combined multi-part geometry — meaning all downstream Regrow/CSB-to-DISES spatial joins match against a farmer's entire landholding footprint, not individual tax parcels.

**`join_dises_shape_table.py`**: Joins the dissolved DISES geometry with the cleaned survey table on farmer ID, adding a `survey_responded` flag.

### 6.6 Regrow ↔ DISES join

**`join_regrow_dises.py`** / **`join_csb_dises.py`**: See the full methodology in [4.5](#45-merging-regrow-dises-and-supplementary-data) / [4.6](#46-merging-csb-dises-and-supplementary-data) — spatial overlay with buffer margin, largest-overlap selection, area/crop match-quality scoring, and representative-field tiering.

### 6.7 Supplementary data joins (Regrow side)

**`rasterize_regrow_to_gSSURGO_grid.py`**: Same rasterization methodology as the CDL grid script, but burning field IDs onto the gSSURGO mukey raster's grid instead — a prerequisite for supplement 8's soil aggregation.

**`join_regrow_supplement_1.py`** (census tract), **`join_regrow_supplement_2.py`** (elevation/slope), **`join_regrow_supplement_3.py`** (watershed), **`join_regrow_supplement_4.py`** (road distance), **`join_regrow_supplement_5.py`** (neighboring-field land management), **`join_regrow_supplement_6.py`** (weather), **`join_regrow_supplement_7.py`** + **`cut_regrow_supplement_7.py`** (crop price), **`join_regrow_supplement_8.py`** (soil, requires `rasterize_regrow_to_gSSURGO_grid.py` above) — see the per-supplement methodology table in [4.4](#44-supplementary-farmland-characteristics).

**`join_regrow_supplement_9.py`** (USDA Census of Agriculture, county-level — keyed on `county_state_name`, not `field_id`; uses `{state}_regrow_supplement_1_table.parquet` to get each county/state pair) — see [4.4](#44-supplementary-farmland-characteristics).

### 6.8 Supplementary data joins (CSB side)

**`join_csb_supplement_1.py`** through **`join_csb_supplement_8.py`**, **`cut_csb_supplement_7.py`**, **`rasterize_CSB_to_gSSURGO_grid.py`** — structurally mirror their Regrow-side counterparts with `CSBID`/`CSBACRES` substituted for `field_id`/`area_acre` and an added outer loop over `CSB_years`. The one methodological difference (supplement 5 lacking a tillage/cover-crop component, since CSB has none) is documented in [4.6.8](#468-merging-csb-with-each-supplementary-dataset).

**`join_csb_supplement_9.py`**: CSB-side mirror of `join_regrow_supplement_9.py` — same county-level methodology, keyed on `county_state_name`, using `{state}_CSB{year}_supplement_1_table.parquet` to get each county/state pair.

---

## 7. Maumee Watershed Filtering

The Maumee River watershed is the largest in the Great Lakes region, consisting of **seven subbasins** that drain into the Western Basin of Lake Erie. To filter fields to the Maumee watershed, use the `subbasin_id` variable from the supplement 3 data file together with the HUC-8 codes for those subbasins.

<br>

**Table 1. HUC-8 codes for the Maumee subbasins**

| HUC-8 code | Watershed name |
|---|---|
| 04100003 | St. Joseph River |
| 04100004 | St. Marys River |
| 04100005 | Upper Maumee River |
| 04100006 | Tiffin River |
| 04100007 | Auglaize River |
| 04100008 | Blanchard River |
| 04100009 | Lower Maumee River |

*Source: Ohio EPA Technical Report AMS/2020-MWN-5.*

<br>

**Figure 1. Maumee River watershed subbasins and Lake Erie drainage area**

![Maumee River watershed map](figures/maumee_watershed_2.png)

*Source: Ohio EPA Technical Report AMS/2020-MWN-5.*

---

## 8. Known Issues

These are code-level findings from a systematic review of the pipeline. As of this update, one item remains open.

### Open

1. **Three download scripts still have a latent division-by-zero / `KeyError` risk.** `download_state_bound.py`, `download_county_bound.py`, and `download_census_tract.py` safely default `Content-Length` to `0` if missing, but then divide by it unconditionally in the file-size-tolerance check, and still index `response.headers["Connection"]` directly. Fix: guard the size check with `if response_size and ...`, and use `response.headers.get("Connection", "")`.