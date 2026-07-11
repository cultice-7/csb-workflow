# CSB-Based Output Files

## Output files

Built around the **CSB field** (`CSBID`) instead of the Regrow field, with an added `{years}` token in every filename (currently `1724`, i.e. the 2017–2024 CSB release window):

| Script(s) | Output file pattern | Contents |
|---|---|---|
| `split_csb_shape_table` | `{state}_CSB{years}_CSBID_geometry.parquet` / `.gpkg` | `CSBID` + geometry only |
| `split_csb_shape_table` | `{state}_CSB{years}_table.parquet` | All CSB attribute columns (`CSBACRES` + one `CDL{year}` main-crop column per year in the release window, `CDL2017`…`CDL2024`), no geometry |
| `join_csb_dises` | `{state}_CSB{years}_dises_table.parquet` | DISES survey variables (suffix `_dises`) + match-quality + representative-field indicators, joined to each CSB field |
| `join_csb_supplement_1` … `join_csb_supplement_8` (+ `cut_csb_supplement_7`) | `{state}_CSB{years}_supplement_{1-8}_table.parquet` | One file per supplementary data block (see supplementary data codes below). As with the Regrow branch, the `_table_reduced` variant of supplement 7 drops nearest-single-elevator price data that is not cleared for distribution |
| `join_csb_supplement_9` | `{state}_CSB{years}_supplement_9_table.parquet` | **County-level**, not field-level — see supplementary data codes below and [Supplementary Farmland Characteristics](../20-datasets/24-supplementary-data.md) |

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
