# CSB-Based Output Files

## Output files

Built around the **CSB field** (`CSBID`) instead of the Regrow field, with an added `{years}` token in every filename (currently `1724`, i.e. the 2017–2024 CSB release window):

| Script(s) | Output file pattern | Contents |
|---|---|---|
| `split_csb_shape_table` | `{state}_CSB{years}_CSBID_geometry.parquet` / `.gpkg` | `CSBID` + geometry only |
| `split_csb_shape_table` | `{state}_CSB{years}_table.parquet` | All CSB attribute columns (`CSBACRES` + one `CDL{year}` main-crop column per year in the release window, `CDL2017`…`CDL2024`), no geometry |
| `join_csb_dises` | `{state}_CSB{years}_dises_table.parquet` | DISES survey variables (suffix `_dises`) + match-quality + representative-field indicators, joined to each CSB field |
| `join_csb_census_tract` | `{state}_CSB{years}_census_tract_table.parquet` | Census tract/county/state boundaries |
| `join_csb_elevation_slope` | `{state}_CSB{years}_elevation_slope_table.parquet` | Elevation & slope |
| `join_csb_watershed` | `{state}_CSB{years}_watershed_table.parquet` | Watershed (HUC8/10/12) |
| `join_csb_nearest_roads` | `{state}_CSB{years}_nearest_roads_table.parquet` | Distance to primary/secondary roads |
| `join_csb_neighbor_field_mgmt` | `{state}_CSB{years}_neighbor_field_mgmt_table.parquet` | Neighboring-field land management |
| `join_csb_weather` | `{state}_CSB{years}_weather_table.parquet` | Weather (PRISM) |
| `join_csb_crop_prices` (+ `cut_csb_crop_prices`) | `{state}_CSB{years}_crop_prices_table.parquet` (+ `_table_reduced` variant) | Crop price (elevator & county level). The `_table_reduced` variant drops nearest-single-elevator price data that is not cleared for distribution |
| `join_csb_soil_composition` | `{state}_CSB{years}_soil_composition_table.parquet` | Soil composition |
| `join_csb_ag_census` | `{state}_CSB{years}_ag_census_table.parquet` | USDA Census of Agriculture — **county-level**, not field-level — see [Supplementary Farmland Characteristics](../20-datasets/24-supplementary-data.md) |
