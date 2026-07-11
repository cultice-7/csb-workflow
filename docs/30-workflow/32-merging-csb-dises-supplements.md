# Merging CSB, DISES, and Supplementary Data

## How are CSB and DISES spatially merged?

Identical methodology to [Merging Regrow with DISES and Supplements](31-merging-regrow-dises-supplements.md), implemented in `join_csb_dises.py`, with `CSBID` as the join key and the merge built around the **CSB field** as the primary output unit. One implementation difference: CSB processing loops over `CSB_years` (currently just `"1724"`) as an outer dimension, since Regrow has no year-window concept.

## What are the outcomes of the matching?

`[TODO — fill in the counts below]`

| State | CSB fields with a DISES match (`field_assigned_dises == "Y"`) | CSB fields with a DISES survey (`survey_responded_dises == "Y"`) | Total CSB fields |
|---|---|---|---|
| IN | | | |
| MI | | | |
| OH | | | |
| **Total** | | | |

Unique DISES holdings with at least one overlapping CSB field: `[ ]` out of `[ ]` total DISES holdings (`[ ]`%).

## How is matching quality evaluated?

Same `area_match_dises` / `crop_match_dises` → `match_quality_dises` (`A`/`B_crop`/`B_area`/`F`) logic as the Regrow side, with one necessary difference: CSB's crop match compares the DISES 2023 crop answer against the single `CDL2023` column (CSB has one annual crop call, not multiple within-year cultivation cycles like Regrow).

## What is a representative-field attribute on the CSB side?

Same purpose and tiering as [the Regrow side](31-merging-regrow-dises-supplements.md#what-is-a-representative-field-attribute-why-do-we-need-it), implemented in `join_csb_dises.py`'s `assign_representative_field()`. As noted above, the CSB version has no crop-confidence dimension (CDL has no per-field ML confidence score), so its tiers are coarser:

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

| Representative field level | Count |
|---|---|
| `RF_level_1_dises` | |
| `RF_level_2_dises` | |
| `RF_level_3_dises` | |
| `RF_level_4_dises` | |
| **Total representative fields** | |

## How to generate CSB sub-samples with available DISES data

Same filter logic as [the Regrow side](31-merging-regrow-dises-supplements.md#how-to-generate-regrow-sub-samples-with-available-dises-data), substituting `CSBID` for `field_id`.

## How to extract DISES-related variables

Same `_dises` suffix convention as the Regrow branch.

## Merging CSB with each supplementary dataset

Structurally near-identical to the Regrow-side scripts (same spatial-join methodology per supplement), with `CSBID` substituted for `field_id` and an added outer loop over `CSB_years`. Note that **Supplement 9** is structurally different from the rest: it's keyed on `county_state_name` rather than `CSBID` — see [Supplementary Farmland Characteristics](../20-datasets/24-supplementary-data.md).

## Output dataset structure and integration

Output files generally follow the pattern `{state}_CSB{years}_<block>_table.parquet` — e.g. `{state}_CSB{years}_table.parquet` for the base table, `{state}_CSB{years}_dises_table.parquet` for the DISES join, and `{state}_CSB{years}_supplement_{n}_table.parquet` for supplement `n`.

To combine multiple blocks into one field-level dataset, merge them on `CSBID` — every block's table has one row per CSB field (fields with no DISES or supplement match still have a row, just with missing values), so `CSBID` is always the key to use when joining any two of these tables together.
