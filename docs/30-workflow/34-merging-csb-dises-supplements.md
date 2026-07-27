# Merging CSB, DISES, and Supplementary Data

## How are CSB and DISES spatially merged?

Implemented in `join_csb_dises.py`. The merge is built **around the CSB field** — CSB field is the primary unit of the output dataset, and DISES attributes are attached to it (not the other way around). The whole procedure is looped once per entry in `CSB_years` (currently just `"1724"`), since CSB (unlike Regrow) is versioned by year-window:

1. Before any spatial join happens, all DISES tax parcels belonging to the same farmer are consolidated into a single record: `clean_dises_shape.py` groups parcels by `comp_id` (the unique farmer identifier) and dissolves each farmer's parcels into one combined **multipolygon** representing their entire landholding footprint (see [DISES Dataset](../20-datasets/23-dises-dataset.md)). This means the unit that CSB fields are actually matched against is "one farmer's combined holdings," not individual tax parcels — so when a CSB field overlaps land belonging to one farmer, that overlap is measured (and the largest-overlap comparison in step 3 below is made) against the farmer's whole multipolygon, not against a single polygon.
2. Each CSB field polygon is slightly buffered (a small negative margin, `buffer_margin` in config) and intersected against all DISES holding multipolygons.
3. Where a CSB field overlaps multiple DISES holdings, only the **largest-overlap** DISES match is kept per CSB field (`CSBID`).
4. Overlap area (in acres) and overlap share (as a percent of the CSB field's own area) are computed from the *original* (unbuffered) geometries.
5. A `parcel_assigned_dises` flag (`Y`/`N`) marks whether any DISES match was found at all.

## What are the outcomes of the matching?

`[TODO — fill in the counts below]`

| State | CSB fields with a DISES match (`parcel_assigned_dises == "Y"`) | CSB fields with a DISES survey (`survey_responded_dises == "Y"`) | Total CSB fields |
|---|---|---|---|
| IN | | | |
| MI | | | |
| OH | | | |
| **Total** | | | |

Unique DISES holdings with at least one overlapping CSB field: `[ ]` out of `[ ]` total DISES holdings (`[ ]`%).

## How is matching quality evaluated?

Same `area_match_dises` / `crop_match_dises` → `match_quality_dises` (`A`/`B_crop`/`B_area`/`F`) logic as the Regrow side, with one necessary difference: CSB's crop match compares the DISES 2023 crop answer against the single `CDL2023` column (CSB has one annual crop call, not multiple within-year cultivation cycles like Regrow).

## What is a representative-field attribute on the CSB side?

Same purpose and tiering as [the Regrow side](33-merging-regrow-dises-supplements.md#what-is-a-representative-field-attribute-why-do-we-need-it), implemented in `join_csb_dises.py`'s `assign_representative_field()`. As noted above, the CSB version has no crop-confidence dimension (CDL has no per-field ML confidence score), so its tiers are coarser:

- **Level 1**: `match_quality_dises == "A"`.
- **Level 2**: `"B_area"` with DISES crop answer missing, or `"B_crop"` with DISES size answer missing.
- **Level 3**: `"B_area"` with a crop mismatch.
- **Level 4**: `"F"` with size missing and a crop mismatch; or both size and crop missing with only a single DISES holding.

## What are the outcomes of the matching quality procedure?

`[TODO — fill in the counts below]`

| `match_quality_dises` | Count |
|---|---|
| `A` | |
| `B_crop` | |
| `B_area` | |
| `F` | |
| **Total matched** | |

| `RF_assignment_dises` | Count |
|---|---|
| `Level 1` | |
| `Level 2` | |
| `Level 3` | |
| `Level 4` | |
| **Total representative fields** | |

## How to generate CSB sub-samples with available DISES data

Same filter logic as [the Regrow side](33-merging-regrow-dises-supplements.md#how-to-generate-regrow-sub-samples-with-available-dises-data), substituting `CSBID` for `field_id`.

## How to extract DISES-related variables

Same `_dises` suffix convention as the Regrow branch.

## Merging CSB with each supplementary dataset

Structurally near-identical to the Regrow-side scripts (same spatial-join methodology per supplement), with `CSBID` substituted for `field_id` and an added outer loop over `CSB_years`. Note that **the USDA Census of Agriculture block** (`join_csb_ag_census`) is structurally different from the rest: it's keyed on `county_state_name` rather than `CSBID` — see [Supplementary Farmland Characteristics](../20-datasets/24-supplementary-data.md).

## Output dataset structure and integration

Output files generally follow the pattern `{state}_CSB{years}_<block>_table.parquet` — e.g. `{state}_CSB{years}_table.parquet` for the base table, `{state}_CSB{years}_dises_table.parquet` for the DISES join, and `{state}_CSB{years}_census_tract_table.parquet` / `{state}_CSB{years}_ag_census_table.parquet` for the supplementary blocks.

To combine multiple blocks into one field-level dataset, merge them on `CSBID` — every block's table has one row per CSB field (fields with no DISES or supplement match still have a row, just with missing values), so `CSBID` is always the key to use when joining any two of these tables together.
