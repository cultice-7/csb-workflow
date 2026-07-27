# Cleaning DISES Data

## How is the raw DISES survey table cleaned?

Implemented in `clean_dises_table.py`. Starting from the raw combined survey export (`combined_data_clean.csv`):

1. Rows that are entirely empty, and rows with a missing `Comprehensive_ID`, are dropped.
2. Key columns are renamed for consistency with the rest of the pipeline: `Comprehensive_ID` → `comp_id`, `field_crop` → `field_crop_23`, `field_tillage` → `field_till_23`, `field_CC` → `field_cover_23`.
3. **Education is unified onto a single scale.** The primary survey distributed during pilot and main waves and shortened follow-up survey for non-respondents asked about education using different scales — a 7-point `Education` scale in the main survey and a 4-point `NR_education` scale in the shortened one. `Education` is first recoded onto `NR_education`'s scale (`1,2→1`; `3,4→2`; `5→3`; `6→4`; `7→5`), then the two are combined into one `Education` column: `Education`'s own (recoded) value is used whenever present, and `NR_education` only fills in where `Education` is missing. The original `NR_education` column is dropped afterward to avoid keeping two names for the same concept.
4. Three farmer-identity indices are computed as row-means of specific Likert-scale Farmer Identity (`FI_*`) survey items: `productivism_index` (`FI_1, FI_4, FI_7, FI_9, FI_10`), `conservationism_index` (`FI_2, FI_3, FI_5, FI_6, FI_8`), `civic_index` (`FI_11`–`FI_16`).

See [DISES Dataset](../20-datasets/23-dises-dataset.md#what-are-the-key-attributes-of-the-dises-dataset) for which of these (and other survey) variables are marked as key farmer characteristics.

## How is the raw DISES parcel shapefile cleaned and consolidated?

Implemented in `clean_dises_shape.py`, starting from the raw parcel shapefile (`DISES_All_Parcels_11.12.25.shp`):

1. Parcels with a placeholder comprehensive-ID value of `0` are dropped — these don't correspond to a real survey respondent.
2. The comprehensive-ID column is renamed to `comp_id`, the unique farmer identifier.
3. Before dissolving, the number of distinct tax parcels per farmer is counted and kept as `n_parcels`.
4. All of a farmer's individual tax parcels are **dissolved** into one combined multi-part geometry per `comp_id` — this is the step that turns "one row per tax parcel" into "one row per farmer's entire landholding footprint" (see [DISES Dataset](../20-datasets/23-dises-dataset.md#what-is-a-dises-parcel)), which is what every downstream Regrow/CSB-to-DISES spatial join actually matches against.
5. The dissolved geometry is reprojected to the project's equal-area CRS, and field area (`area_acre`) is computed from it.
6. `n_parcels` is merged back on, and the result — one row per farmer, with `comp_id`, dissolved `geometry`, `area_acre`, and `n_parcels` — is saved as `DISES_shape_cleaned.parquet`, the input to `join_dises_shape_table.py`.
