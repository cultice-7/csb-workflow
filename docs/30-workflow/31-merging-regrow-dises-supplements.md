# Merging Regrow, DISES, and Supplementary Data

## How are Regrow and DISES spatially merged?

Implemented in `join_regrow_dises.py`. The merge is built **around the Regrow field** — Regrow field is the primary unit of the output dataset, and DISES attributes are attached to it (not the other way around):

1. Each Regrow field polygon is slightly buffered (a small negative margin, `buffer_margin` in config) and intersected against all DISES parcel-holding polygons.
2. Where a Regrow field overlaps multiple DISES holdings, only the **largest-overlap** DISES match is kept per Regrow field.
3. Overlap area (in acres) and overlap share (as a percent of the Regrow field's own area) are computed from the *original* (unbuffered) geometries.
4. A `field_assigned_dises` flag (`Y`/`N`) marks whether any DISES match was found at all.

## What are the outcomes of the matching?

`[TODO — fill in the counts below]`

| State | Regrow fields with a DISES match (`field_assigned_dises == "Y"`) | Regrow fields with a DISES survey (`survey_responded_dises == "Y"`) | Total Regrow fields
|---|---|---|---|
| IN | | | |
| MI | | | |
| OH | | | |
| **Total** | | | |

Unique DISES holdings with at least one overlapping Regrow field: `[ ]` out of `[ ]` total DISES holdings (`[ ]`%).

## How is matching quality evaluated?

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

## What is a "representative field" attribute? Why do we need it?

The DISES survey asks farmers a block of questions about their "representative field," without specifying *which* field they meant. Since DISES holdings are dissolved per-farmer (see [DISES Dataset](../20-datasets/23-dises-dataset.md)) and can correspond to multiple candidate Regrow fields, this attribute identifies the most plausible candidate by comparing each candidate Regrow field's area and 2023 main crop against the farmer's representative-field survey answers.

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

## What are the outcomes of the matching quality procedure?

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

## How to generate Regrow sub-samples with available DISES data

| Sub-sample | Filter |
|---|---|
| Regrow fields with an assigned DISES field (survey or not) | `field_assigned_dises == "Y"` |
| Regrow fields with an assigned DISES field **and** survey data | `field_assigned_dises == "Y"` and `survey_responded_dises == "Y"` |
| Regrow fields matching the assigned DISES field by size or crop | the above, **and** `match_quality_dises` is one of `"A"`, `"B_area"`, `"B_crop"` |

## How to extract DISES-related variables

Every column originating from or derived from the DISES dataset carries the suffix `_dises` (e.g. `field_crop_23_dises`, `match_quality_dises`, `RF_level_1_dises`).

## Merging Regrow with each supplementary dataset

Each supplementary join is implemented as its own script and produces its own standalone output table (joined to the host dataset only by `field_id`, not chained to each other) — see [Supplementary Farmland Characteristics](../20-datasets/24-supplementary-data.md) for the per-supplement methodology, and [Script Reference](../40-scripts/01-overview.md) for full script-level detail. The one exception is **supplement 9**, which is keyed on `county_state_name` rather than `field_id` — see [Supplementary Farmland Characteristics](../20-datasets/24-supplementary-data.md) for how it's built and matched.

## Output dataset structure and integration

Output files generally follow the pattern `{state}_regrow_<block>_table.parquet` — one file per state per data block, e.g. `{state}_regrow_table.parquet` for the base table, `{state}_regrow_dises_table.parquet` for the DISES join, and `{state}_regrow_supplement_{n}_table.parquet` for supplement `n`.

To combine multiple blocks into one field-level dataset, merge them on `field_id` — every block's table has one row per Regrow field (fields with no DISES or supplement match still have a row, just with missing values), so `field_id` is always the key to use when joining any two of these tables together.
