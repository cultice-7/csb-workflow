# USDA Crop Sequence Boundaries (CSB) Dataset

## What is the CSB dataset?

CSB is a **public** USDA dataset (developed with USDA's Economic Research Service) that produces field boundary, crop acreage, and crop-rotation estimates across the contiguous United States. CSB boundaries are derived algorithmically from the **Cropland Data Layer (CDL)** — an open-source method (available on GitHub) that analyzes multi-year crop-rotation consistency in CDL pixels to infer field boundaries. CSB releases are published annually, each covering a rolling 8-year window; this project uses the 2017–2024 release. Because CSB boundaries are derived purely from satellite-based crop classification, they contain **no personally identifying information** and represent "crops grown only," not ownership or tax-parcel boundaries.

## The definition of a "CSB field"

A CSB field (`CSBID`) is a polygon produced by the CSB algorithm's crop-rotation-consistency analysis over CDL pixels for a given multi-year window — conceptually similar to a Regrow field (both are remote-sensing-derived, not ownership-derived) but built with a different methodology and, in practice, different boundaries on the ground for the same physical land.

## What is the geographic span of CSB data?

CSB itself covers the contiguous United States; within this project, CSB data is processed for the same 7 states as Regrow: Illinois, Indiana, Iowa, Michigan, Minnesota, Ohio, Wisconsin, with data available for all Michigan counties.

## What data is available within the CSB dataset?

Each CSB field has one main (commodity) crop classification per calendar year, for every year in the release window (`CDL2017`…`CDL2024`), plus the field's acreage (`CSBACRES`) and geographic location.

## Does the CSB dataset have data on tillage intensity and cover crops?

No. CSB contains **main crop rotation only** — no tillage intensity, no cover crop indicator, and no within-year practice timeline. This is the single largest informational gap relative to Regrow.

## What land management practices are available in the CSB dataset?

Only one: **main (commodity) crop**, classified once per calendar year per field (see "How does CSB classify main crop?" below). There is no post-harvest tillage, cover crop, or pre-plant tillage equivalent.

## How are raw and processed CSB datasets structured?

Raw CSB data is distributed as a single national geodatabase per release window (e.g. `CSB1724.gdb`), containing one polygon per field with one `CDL{year}` column per year in the window. USDA publishes CSB in rolling 8-year windows (e.g. 2016–2023, 2017–2024, 2018–2025), each released as a new geodatabase; this pipeline is currently configured for the 2017–2024 window (`CSB1724`, set via `CSB_years` in `config.yml`). `clip_csb_shape.py` extracts each state's portion via an attribute filter on `STATEFIPS`, and `split_csb_shape_table.py` splits each state's file into a geometry-only file (`CSBID` + polygon), an attributes-only file (`CSBID`, `CSBACRES`, and the `CDL{year}` columns), and a combined shape+table file. Unlike Regrow, there is no long-to-wide reshaping step, since CSB's raw structure (one row per field, one column per year) is already in the project's target wide format.

## How does CSB classify main crop?

CSB assigns one annual crop code per field per year (`CDL2017`…`CDL2024`), taken as the majority CDL pixel class within that field's polygon footprint for that year — i.e., one classification per calendar year, not per cultivation cycle. The crop-code scheme is the standard USDA CDL numeric classification (0=background, 1–60 major row/specialty crops, 61–65 non-crop/fallow/pasture/forest/shrub, 66–80 tree crops, 81–109 other land cover, 110–195 NLCD-derived classes, 204–254 specialty crops) — the same scheme Regrow's crosswalk targets, which is what makes Regrow and CSB main-crop codes directly comparable.

## How is CSB data validated by USDA?

`[TODO: ADD THIS INFORMATION]`

## What are the key differences between Regrow and CSB?

| | Regrow | CSB |
| --- | --- | --- |
| License | Proprietary | Public |
| Tillage | Yes (3 classes, two windows/year) | No |
| Cover crop | Yes (presence indicator, 3 classes) | No |
| Main crop | Yes, per cultivation cycle (can be multiple/year) | Yes, one classification per calendar year |
| Timeline detail | Planting/harvest dates per cycle | Annual only |
| Year coverage (this project) | 2014–2025 | 2017–2024 |
| Michigan coverage | Limited to 2 counties adjacent to Ohio | All counties |

Both datasets are conceptually similar (remote-sensing-derived field boundaries with crop classification) but use different methodologies, so the same physical field will generally have different boundaries and possibly different crop calls in the two datasets.

## What are the key attributes of the CSB dataset?

Field size (`CSBACRES`) and the main crop classification for each year in the release window (`CDL{year}`), alongside the field's unique identifier (`CSBID`) and geographic location (geometry).

## What naming conventions are used for CSB variables?

The outcome variable is `CDL{year}` for `year ∈ {2017...2024}`.

## Where can I find more information about CSB dataset attributes?

`[TODO — paste links]`

- Crop Sequence Boundaries Metadata:
- CSB1724_DISES_Supplementary_Data Codebook:
