# csb-workflow

**Land-Based Carbon SAS Project — Data Pipeline**

## Table of Contents

- [Overview](#overview)
- [Detailed Documentation](#detailed-documentation)
- [Support and Feedback](#support-and-feedback)
- [Authors](#authors)
- [License](#license)
- [Current Project Status](#current-project-status)

## Overview

This repository implements the SAS Project data pipeline: a Snakemake-orchestrated workflow that builds agricultural field-level datasets describing land management practices, farmer characteristics, and farmland conditions across a 7-state Midwest study region for a multi-year period.

This project provides field-level data on land management activities carried out by individual farmers, farmer characteristics, and supplementary data that may influence farmers' decision-making regarding adoption of specific land management practices.

### Two parallel datasets: Regrow-based and CSB-based

There is no single dataset that offers both rich land-management detail and unrestricted public distribution. **Regrow** is a proprietary, remote-sensing-derived dataset with detailed tillage, cover crop, and main-crop information, but its distribution is restricted. **CSB** is a public USDA dataset with main-crop information only, but no licensing restrictions. Rather than committing to a single source, this project develops the pipeline as **two parallel branches**, one built around Regrow fields and one built around CSB fields. Each branch is independently merged with the same DISES farmer-survey data and the same supplementary farmland characteristics, using a mirrored (but separately implemented) set of scripts. This lets the team compare results across both field-delineation sources, and choose — or combine — whichever is more suitable for a given downstream analysis, rather than being locked into the limitations of either source alone.

- Land management activity data are dynamic from **2014–2025** (Regrow) or **2017–2024** (CSB).
- Farmer characteristics (DISES survey) are available for a **single reference year (2023 season, collected early 2024)**.
- Supplementary farmland characteristics cover 9 topics: census geography, elevation/slope, watershed geography, road proximity, neighboring-field land management, weather, crop price, soil composition, and USDA Census of Agriculture.

The project produces **two parallel, independently usable output datasets**, built around two different field-delineation sources:

| | Regrow-based dataset | CSB-based dataset |
|---|---|---|
| Field unit | "Regrow field" (proprietary ML-delineated boundary) | "CSB field" (USDA algorithm-delineated boundary, from CDL) |
| Land management detail | Tillage, cover crop, main crop, full cultivation-cycle timeline | Main crop rotation only |
| License | Proprietary (Regrow Ag) | Public (USDA NASS/ERS) |
| Primary key | `field_id` | `CSBID` |
| Timeline | 2014–2025 | 2017–2024 |

Both datasets are merged with DISES farmer-survey data and the same 9 blocks of supplementary farmland characteristics, using a parallel but separately-implemented set of scripts.

---

## Detailed Documentation

<details>
<summary><strong>Setup &amp; Getting Started</strong> — Software prerequisites, cloning the repository, setting up the conda environment, and obtaining proprietary data</summary>

This section covers everything needed to get up and running with the pipeline for the first time. See the [Setup overview](00-setup/00-overview.md) for a full table of contents.

| File | What it answers |
|---|---|
| [Technological Prerequisites](00-setup/01-technological-prerequisites.md) | What software do I need to install before I can use this pipeline? |
| [How to Clone the Repository](00-setup/02-how-to-clone-repository.md) | How do I get a local copy of this repository using Git or VS Code? |
| [Repository Organization](00-setup/03-repository-organization.md) | How is this repository organized? |
| [How to Get Access to the Fileshare](00-setup/04-how-to-get-access-to-fileshare.md) | How do I get access to the protected project fileshare that holds the proprietary raw data? |
| [Which Proprietary Data to Download](00-setup/05-which-proprietary-data-to-download.md) | Which datasets must I download manually from the fileshare, and where do I place them in the project folder? |
| [How to Set Up the Conda Environment](00-setup/06-how-to-set-environment.md) | How do I create and activate the `csb_workflow` conda environment? |

</details>

<details>
<summary><strong>Output Datasets</strong> — What output files the pipeline produces, what each file contains, where schemas and codebooks live, and how to run the pipeline to generate outputs</summary>

This section describes the structure of all output files and how to generate them using Snakemake. See the [Output overview](10-output/10-overview.md) for a full table of contents.

| File | What it answers |
|---|---|
| [Regrow-Based Output](10-output/11-regrow-based-output.md) | What output files does the Regrow branch produce, and what is in each? |
| [CSB-Based Output](10-output/12-csb-based-output.md) | What output files does the CSB branch produce, and what supplementary data codes are used? |
| [Data Schema Links](10-output/13-data-schema-links.md) | Where are the codebooks, schema references, and dataset overview slides? |
| [How to Obtain Output Datasets](10-output/14-how-to-obtain-output-datasets.md) | How do I run Snakemake to build all or selected outputs? What are the recommended execution chains? |

</details>

<details>
<summary><strong>Dataset Documentation</strong> — Detailed FAQ on the Regrow, CSB, DISES, and supplementary datasets: what they are, how they are structured, and what variables they contain</summary>

This section is the primary reference for understanding the datasets that feed into and come out of this pipeline. See the [Datasets overview](20-datasets/20-overview.md) for a full table of contents.

| File | What it answers |
|---|---|
| [Regrow Dataset](20-datasets/21-regrow-dataset.md) | What is Regrow? How does it detect tillage, cover crops, and commodity crops? How is confidence defined? How is the raw data restructured into the final wide-table format? |
| [CSB Dataset](20-datasets/22-csb-dataset.md) | What is CSB? How does it classify crops? How does its structure and coverage compare to Regrow? |
| [DISES Dataset](20-datasets/23-dises-dataset.md) | What is DISES? Who are the respondents, and what geographic area do they cover? What is a DISES parcel, and how does it differ from a Regrow field? |
| [Supplementary Farmland Characteristics](20-datasets/24-supplementary-data.md) | What are the 9 supplementary data blocks (census geography, elevation, watershed, roads, neighboring fields, weather, crop prices, soil, Census of Agriculture), their sources, and how each is computed? |

</details>

<details>
<summary><strong>Workflow: Merging Datasets</strong> — How Regrow and CSB fields are spatially merged with DISES farmer data and supplementary characteristics, including match-quality scoring and representative-field assignment</summary>

This section explains the spatial merge methodology, how match quality is evaluated, what the representative-field attribute means, and how to combine output blocks into a single analysis dataset. See the [Workflow overview](30-workflow/30-overview.md) for a full table of contents.

| File | What it answers |
|---|---|
| [Merging Regrow with DISES and Supplements](30-workflow/31-merging-regrow-dises-supplements.md) | How are Regrow field polygons spatially matched to DISES farmer parcels, and how is the best match selected when a field overlaps multiple holdings? How is match quality scored? What is the representative-field attribute and how is it tiered? How are the 9 supplementary farmland characteristic blocks joined to Regrow fields, and how do I combine those standalone tables into a single analysis dataset? |
| [Merging CSB with DISES and Supplements](30-workflow/32-merging-csb-dises-supplements.md) | How are CSB field polygons spatially matched to DISES farmer parcels, and how is the best match selected when a field overlaps multiple holdings? How is match quality scored? What is the representative-field attribute and how is it adapted for CSB fields? How are the 9 supplementary farmland characteristic blocks joined to CSB fields, and how do I combine those standalone tables into a single analysis dataset? |

</details>

<details>
<summary><strong>Script Reference</strong> — Explanation of every active script in <code>scripts/</code>, organized by pipeline stage: what each script does, key steps, and assumptions</summary>

This section catalogs every active script in `scripts/` (excluding `scripts/archive/`, which holds superseded versions kept for reference only). See the [Scripts overview](40-scripts/40-overview.md) for a full table of contents.

| File | Scripts covered |
|---|---|
| [Download Scripts](40-scripts/41-download-scripts.md) | `download_census_tract.py`, `download_state_bound.py`, `download_county_bound.py`, `download_roads.py`, `download_watershed.py`, `download_weather.py`, `download_csb.py`, `download_cdl.py`, `download_elevation.py` |
| [Raw Data Processing Scripts](40-scripts/42-raw-data-processing-scripts.md) | `reproject_elevation.py`, `calculate_slope.py`, `clip_elevation_slope.py`, `clip_gSSURGO_mukey_rasters.py`, `clip_cdl_rasters.py`, `clean_crop_prices.py`, `clip_reproject_weather_rasters.py`, `clip_csb_shape.py`, `split_csb_shape_table.py`, `check_csb_shape_table.py`, `check_regrow_shape_table.py` |
| [Regrow Preparation Scripts](40-scripts/43-regrow-preparation-scripts.md) | `join_regrow_2025_updates.py`, `split_save_regrow_geometry.py`, `split_regrow_monitor_by_state.py`, `clean_regrow_table.py`, `join_regrow_shape_table.py` |
| [Regrow Validation Scripts](40-scripts/44-regrow-validation-scripts.md) | `rasterize_regrow_to_CDL_grid.py`, `join_regrow_cdl.py`, `validate_regrow_crop_with_cdl.py` |
| [DISES Preparation Scripts](40-scripts/45-dises-preparation-scripts.md) | `clean_dises_table.py`, `clean_dises_shape.py`, `join_dises_shape_table.py` |
| [Regrow–DISES Merge Script](40-scripts/46-merge-regrow-dises-scripts.md) | `join_regrow_dises.py` |
| [Regrow Supplementary Data Merge Scripts](40-scripts/47-merge-regrow-supplements-scripts.md) | `rasterize_regrow_to_gSSURGO_grid.py`, `join_regrow_census_tract.py`, `join_regrow_elevation_slope.py`, `join_regrow_watershed.py`, `join_regrow_nearest_roads.py`, `join_regrow_neighbor_field_mgmt.py`, `join_regrow_weather.py`, `join_regrow_crop_prices.py`, `cut_regrow_crop_prices.py`, `join_regrow_soil_composition.py`, `join_regrow_ag_census.py` |
| [CSB–DISES Merge Script](40-scripts/48-merge-csb-dises-scripts.md) | `join_csb_dises.py` |
| [CSB Supplementary Data Merge Scripts](40-scripts/49-merge-csb-supplements-scripts.md) | `rasterize_CSB_to_gSSURGO_grid.py`, `join_csb_census_tract.py`, `join_csb_elevation_slope.py`, `join_csb_watershed.py`, `join_csb_nearest_roads.py`, `join_csb_neighbor_field_mgmt.py`, `join_csb_weather.py`, `join_csb_crop_prices.py`, `cut_csb_crop_prices.py`, `join_csb_soil_composition.py`, `join_csb_ag_census.py` |

</details>

<details>
<summary><strong>Miscellaneous</strong> — Maumee watershed filtering guide and known pipeline issues</summary>

See the [Miscellaneous overview](50-miscellaneous/50-overview.md) for a full table of contents.

| File | What it answers |
|---|---|
| [How to Filter to the Maumee Watershed](50-miscellaneous/51-how-to-filter-Maumee-watershed.md) | Which HUC-8 codes define the Maumee River watershed study area? How do I subset my data to this watershed? |
| [Known Issues](50-miscellaneous/52-known-issues.md) | What open bugs or edge-case limitations are currently known in the pipeline code? |

</details>

---

## Support and Feedback

> **[BRIAN TO FILL IN]**: instructions on how to post an issue in the project's GitHub repository. Plus, point of contact for exceptional circumstances
---

## Authors

This data pipeline was designed and created by the SAS Project team at the Ohio State University. We would like to express our gratitude to all contributors participating in its development, maintenance, and continuous improvement.

Repository maintainers:  
Brian Cultice, Ohio State University: <cultice.7@osu.edu>  
Grigoriy Baranov, Ohio State University: <baranov.4@osu.edu>

---

## License

---

## Current Project Status

> **Status of this document:** This is the first release of the entirely updated README supporting the code pipeline developed for the SAS project. The documentation is assembled from the intital project README (10Mar2026), the representative-field methodology memo, the script order memo, the Snakefile file, and the key methodological steps implemeted in active scripts in `scripts/`.
