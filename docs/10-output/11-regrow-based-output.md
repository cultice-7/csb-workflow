# Regrow-Based Output Files

The pipeline produces output in two parallel branches — one based on **Regrow**, one based on **CSB** — each centered on its own field-delineation unit. Within each branch, most merge operations produce their own standalone output file (rather than one monolithic table), so that users can pull in only the supplementary blocks they need.

## Output files

| Script(s) that produce the file | Output file pattern | Contents |
|---|---|---|
| `split_save_regrow_geometry` | `{state}_regrow_fieldID_geometry.parquet` / `.gpkg` | `field_id` + geometry only |
| `join_regrow_shape_table` | `{state}_regrow_table.parquet` | All Regrow attribute columns (tillage, cover crop, main crop, cultivation cycle dates/confidence), no geometry |
| `join_regrow_dises` | `{state}_regrow_dises_table.parquet` | DISES survey variables (suffix `_dises`) + match-quality + representative-field indicators, joined to each Regrow field |
| `join_regrow_census_tract` | `{state}_regrow_census_tract_table.parquet` | Census tract/county/state boundaries |
| `join_regrow_elevation_slope` | `{state}_regrow_elevation_slope_table.parquet` | Elevation & slope |
| `join_regrow_watershed` | `{state}_regrow_watershed_table.parquet` | Watershed (HUC8/10/12) |
| `join_regrow_nearest_roads` | `{state}_regrow_nearest_roads_table.parquet` | Distance to primary/secondary roads |
| `join_regrow_neighbor_field_mgmt` | `{state}_regrow_neighbor_field_mgmt_table.parquet` | Neighboring-field land management |
| `join_regrow_weather` | `{state}_regrow_weather_table.parquet` | Weather (PRISM) |
| `join_regrow_crop_prices` (+ `cut_regrow_crop_prices`) | `{state}_regrow_crop_prices_table.parquet` (+ `_table_reduced` variant) | Crop price (elevator & county level). The `_table_reduced` variant exists because the full file contains nearest-single-elevator price data that is sensitive and not cleared for distribution — `cut_regrow_crop_prices.py` drops those columns before the reduced file is shared |
| `join_regrow_soil_composition` | `{state}_regrow_soil_composition_table.parquet` | Soil composition |
| `join_regrow_ag_census` | `{state}_regrow_ag_census_table.parquet` | USDA Census of Agriculture — **county-level**, not field-level — see [Supplementary Farmland Characteristics](../20-datasets/24-supplementary-data.md) |
| `join_regrow_cdl` / `validate_regrow_crop_with_cdl` | `{state}_regrow_cdl_validation_table.parquet`, `{state}_regrow_cdl_summary_by_crop_category.xlsx` | Internal validation of Regrow main-crop calls against USDA CDL — not a "supplementary data" file per se, but a QA output |

States covered: `IA, IL, IN, MI, MN, OH, WI`. DISES join is only produced for `IN, MI, OH` (the states where DISES survey data exists).
