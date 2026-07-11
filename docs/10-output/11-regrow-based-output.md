# Regrow-Based Output Files

The pipeline produces output in two parallel branches — one based on **Regrow**, one based on **CSB** — each centered on its own field-delineation unit. Within each branch, most merge operations produce their own standalone output file (rather than one monolithic table), so that users can pull in only the supplementary blocks they need.

## Output files

| Script(s) that produce the file | Output file pattern | Contents |
|---|---|---|
| `split_save_regrow_geometry` | `{state}_regrow_fieldID_geometry.parquet` / `.gpkg` | `field_id` + geometry only |
| `join_regrow_shape_table` | `{state}_regrow_table.parquet` | All Regrow attribute columns (tillage, cover crop, main crop, cultivation cycle dates/confidence), no geometry |
| `join_regrow_dises` | `{state}_regrow_dises_table.parquet` | DISES survey variables (suffix `_dises`) + match-quality + representative-field indicators, joined to each Regrow field |
| `join_regrow_supplement_1` … `join_regrow_supplement_8` (+ `cut_regrow_supplement_7`) | `{state}_regrow_supplement_{1-8}_table.parquet` (supplement 7 also has a `_table_reduced` variant) | One file per supplementary data block (see supplementary data codes below). The `_table_reduced` variant of supplement 7 exists because the full file contains nearest-single-elevator price data that is sensitive and not cleared for distribution — `cut_regrow_supplement_7.py` drops those columns before the reduced file is shared |
| `join_regrow_supplement_9` | `{state}_regrow_supplement_9_table.parquet` | **County-level**, not field-level — see supplementary data codes below and [Supplementary Farmland Characteristics](../20-datasets/24-supplementary-data.md) |
| `join_regrow_cdl` / `validate_regrow_crop_with_cdl` | `{state}_regrow_cdl_validation.parquet`, `{state}_regrow_cdl_summary_by_crop_category.xlsx` | Internal validation of Regrow main-crop calls against USDA CDL — not a "supplementary data" file per se, but a QA output |

States covered: `IA, IL, IN, MI, MN, OH, WI`. DISES join is only produced for `IN, MI, OH` (the states where DISES survey data exists).

## Supplementary data codes

Supplementary blocks are referenced by numeric code throughout both the Regrow-based and CSB-based output file names.

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
