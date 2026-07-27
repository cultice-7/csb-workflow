# Cleaning and Processing Regrow Data

## How is 2025 Regrow data combined with the 2014-2024 data?

Implemented in `join_regrow_2025_updates.py`. Regrow originally delivered two separate raw exports — a legacy batch covering 2014–2024, and a separate 2025 update — which this script combines into one continuous file per state before any other processing happens:

- **Field geometry**: for each state, the 2014–2024 boundary GeoJSON and the 2025 boundary GeoJSON are outer-merged on `boundary_id`. Where a boundary exists in both exports, the 2014–2024 geometry is kept; only when a boundary is missing from the 2014–2024 export does the script fall back to the 2025 geometry (`combine_first`). The result is saved as `{state}_field_boundaries_2014-2025.parquet`.
- **Monitor practice data**: rather than one CSV per state, Regrow's 2014–2024 raw monitor exports bundle some smaller states together (`states_monitor` in config, e.g. `MN_WI_IA_IN_IL` as one combined file, with `OH` and `MI` each on their own). For each such bundle, the script reads the 2014–2024 CSV once, then loops over the individual states within that bundle, reading each one's separate 2025 CSV and appending it as additional rows (`pd.concat`) — since the 2025 update is simply more history for the same fields/practices, not a new set of columns. The result is saved as `Monitor_data_{state_monitor}_2014-2025.parquet`.

This combined geometry and combined monitor table are the inputs to `split_save_regrow_geometry.py` and `clean_regrow_table.py`, respectively.

## How is Regrow field geometry extracted and cleaned?

Implemented in `split_save_regrow_geometry.py`. For each state:

1. The combined 2014–2025 geometry is reprojected to the project's equal-area CRS (`target_CRS`).
2. `geometry.buffer(0)` is applied — a standard trick to repair invalid geometries (e.g. self-intersections) without altering geometries that are already valid.
3. `boundary_id` is renamed to `field_id`.
4. Only `field_id` and `geometry` are kept, saved as `{state}_regrow_fieldID_geometry.parquet` (and a `.gpkg` copy). This geometry-only file, keyed on `field_id`, is the canonical Regrow field boundary used by every later spatial join in the pipeline.

## Transformation of the Regrow land management variables

`clean_regrow_table.py` recodes each of Regrow's four land management variables — main crop, post-harvest tillage, pre-plant tillage, and cover crop — from Regrow's raw string labels into small integer codes, using three separate hardcoded crosswalks.

**Main crop** (`main_crop_map`): Regrow's string crop names (e.g. `"corn"`, `"soybean"`, `"wheat_winter"`) are recoded into integer codes deliberately aligned with the **USDA CDL/CSB numeric crop code scheme** (e.g. corn=1, soybean=5, winter wheat=24), so that Regrow crop calls can be directly compared against CDL pixel values without further remapping, and so that Regrow and CSB main-crop codes are directly comparable to each other. One sentinel value, `999`, is used for `non_cropland`/`berry` — `999` is **not** a real CDL code; it's a pipeline-internal "doesn't map cleanly" marker. Building this crosswalk also serves a data-quality purpose beyond enabling CDL comparison: it requires explicitly enumerating every distinct crop-name string Regrow uses, which surfaces cases where the same underlying crop is recorded under two different raw names (e.g. `"flax"` and `"flaxseed"`, or `"evergreen"` and `"forest_evergreen"`) so they can be deliberately mapped to the same numeric code rather than being silently treated as different crops downstream.

**Tillage** (`tillage_map`, applied identically to both post-harvest `PHtill` and pre-plant `PPtill`): Regrow's raw tillage-intensity labels — `TILLAGE_INTENSITY_NO_TILLAGE`, `TILLAGE_INTENSITY_REDUCED_TILLAGE`, `TILLAGE_INTENSITY_CONVENTIONAL_TILLAGE` — are recoded to `1`/`2`/`3` respectively, and `TILLAGE_INTENSITY_NO_TILLAGE_DATA` (no usable observation) is recoded to missing. Because a single cultivation cycle can have more than one recorded tillage event of the same type (see below), the same map also covers same-or-mixed-intensity **pairs** of raw labels, recoded to codes `4`–`6` (see the pair-code table in the next section).

**Cover crop** (`cover_crop_map`): Regrow's raw cover-crop class labels are recoded as `1` = cover crop not tracked, `2` = potential cover crop, `3` = detected cover crop, `999` = not applicable (a perennial crop is grown on the field that year, so a cover crop doesn't apply), with the raw "no data" label recoded to missing.

## How are raw and processed Regrow datasets structured? What is a cultivation cycle?

Raw Regrow data arrives as a **long table, one row per cultivation cycle per field**. `clean_regrow_table.py` reorganizes this into the project's final format: a **wide table, one row per field**, with practice/crop values spread across columns indexed by calendar year and cultivation cycle number.

A **cultivation cycle** is the unit that raw Regrow data is built around: it groups the full sequence of land management activities — post-harvest tillage → cover crop → pre-plant tillage → main crop — that contribute to growing **one** commodity crop, following a **harvest-to-harvest model** (not a calendar-year model). A cycle starts the day after the previous crop's harvest and ends at the next commodity's harvest. Because a cultivation cycle is defined around a single commodity crop by construction, **a field cannot have more than one main crop within a single cultivation cycle** — that scenario is excluded by the definition of a cultivation cycle. What *can* happen is multiple complete cultivation cycles occurring within the same calendar year (e.g. two crops grown in succession in one growing season); the pipeline handles this by adding a cultivation-cycle-number suffix to the variable name, rather than by allowing multiple crops inside one cycle. Fallow periods are themselves reported as a "fallow" cultivation cycle, so the data structure stays uniform even when a field sits idle. Source: [Monitor Cultivation Cycles](https://help.regrow.ag/monitor-cultivation-cycles).

Reorganizing the raw long table into the final wide table involves three steps:

1. **Assign each cultivation cycle a calendar year.** The pipeline uses `end_year` — the calendar year in which that cycle's main crop was harvested — as the cycle's year label. `end_year` is used (rather than, say, the cycle's start date) because the harvest date is the one fixed, unambiguous point in a harvest-to-harvest cycle: a cycle's *start* is just whenever the *previous* crop happened to be harvested, which could fall in the prior calendar year, but its *end* is always the harvest of the crop the cycle is actually about.
2. **Number the cultivation cycles within each field-year.** Within a given field and `end_year`, cycles are numbered in chronological order as `cycle_count` (1st, 2nd, or 3rd cycle ending in that calendar year — capped at 3). This `cycle_count` is the mechanism that captures multiple cultivation cycles landing in the same calendar year, since (per above) it's never the case that one cycle itself contains more than one main crop.
3. **Pivot from long to wide table.** The long table is pivoted so that each `(activity, end_year, cycle_count)` combination becomes its own column, producing the variable-naming convention used throughout this project:

   ```text
   {land management activity}_{year}_{cultivation cycle number}
   ```

   where `land management activity ∈ {PHtill, cover, PPtill, crop}`, `year` is the 2-digit `end_year`, and `cultivation cycle number ∈ {1, 2, 3}`.

**Handling multiple records within the same cycle.** Although a cycle can only have one main crop, it can have more than one *recorded tillage event* of the same type (e.g. a field tilled twice within one cycle's post-harvest window). Before pivoting, such repeated records are aggregated using these rules: if the repeated records share identical dates and identical values, the value is kept as a single observation; if dates are identical but the recorded values differ, the result is set to missing (the conflicting records can't both be right); if the dates genuinely differ (i.e. these are two distinct tillage events, not duplicate records of the same event), the two values are joined into one **compound code** representing the pair — for example, a field tilled twice within one cycle, once no-till and once reduced-till, is recorded as a single combined tillage code (5) representing that specific pair, rather than as two separate tillage columns. The full set of these paired codes is in the table below.

Because more than one tillage *event* can occur within a single cycle, the tillage code is not a strict 1/2/3 ordinal scale — it also encodes same-or-mixed-intensity *pairs*:

| Code | Meaning |
| --- | --- |
| 1 | No-till (single event, or two no-till events) |
| 2 | Reduced till (single event, or no-till + reduced-till pair) |
| 3 | Conventional till (single event, or no-till + conventional-till pair) |
| 4 | Reduced-till + reduced-till pair |
| 5 | Reduced-till + conventional-till pair |
| 6 | Conventional-till + conventional-till pair |

**Cover crop codes**: `1` = no cover crop, `2` = potential cover crop, `3` = detected cover crop, `999` = cover crop not applicable.

Missing tillage/cover/crop code columns that don't exist for every state (e.g. a rare second cultivation cycle in a given year) are backfilled as nullable `Int16` so that column dtypes stay consistent across all state files, rather than accidentally upcasting to `float64`; their `_start`/`_end` companions are backfilled as `datetime64`, and `_conf` companions as `float64`.

## How are the Regrow geometry and attribute tables joined?

Implemented in `join_regrow_shape_table.py`. For each state, the cleaned geometry-only file (`{state}_regrow_fieldID_geometry.parquet`) and the cleaned wide attribute table (`{state}_regrow_monitor_wide_coded.parquet`, produced by `clean_regrow_table.py`) are merged on `field_id`. Field area in acres (`area_acre`) is computed from the geometry at this point. The combined result is saved both as a full shape+table file (`{state}_regrow_shape_table.parquet`) and, with geometry dropped, as the attribute-only `{state}_regrow_table.parquet` used throughout the rest of the pipeline (DISES join, supplementary joins, CDL validation).
