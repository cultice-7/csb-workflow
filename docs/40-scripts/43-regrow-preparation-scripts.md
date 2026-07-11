# Regrow Data Preparation Scripts

**`join_regrow_2025_updates.py`**: Merges legacy (2014–2024) and newly-released (2025) raw Regrow exports before any other processing — combining geometry files (preferring the legacy boundary where both exist) and concatenating the tabular Monitor data.

**`split_save_regrow_geometry.py`**: Produces the canonical `{state}_regrow_fieldID_geometry.parquet`/`.gpkg` files: reprojects to the target CRS, repairs any invalid geometries, and renames the raw `boundary_id` key to the project's canonical `field_id`.

**`split_regrow_monitor_by_state.py`**: Splits the legacy era's multi-state-concatenated Monitor tables back out into true per-state files, by matching each state's own set of field IDs.

**`clean_regrow_table.py`**: The core Regrow cleaning/reshaping script — see [Regrow Dataset](../20-datasets/21-regrow-dataset.md) for the full walkthrough ("How are raw and processed Regrow datasets structured?" and "Transformation of the Regrow commodity crop variable"). In brief: aggregates multiple land-management records within a cultivation cycle, assigns each cycle a calendar year and cycle count (capped at 3/year), pivots from long to the final wide `{activity}_{year}_{cycle}` format, and recodes crop/tillage/cover-crop strings to the project's numeric coding schemes. Also unifies the column schema across all states at the end, so every state's output file has an identical (alphabetically sorted) set of columns even if, e.g., one state never recorded a 3rd cultivation cycle in any year.

**`join_regrow_shape_table.py`**: Joins the canonical field geometry with the cleaned wide attribute table on `field_id`, and computes `area_acre` from the polygon area.
