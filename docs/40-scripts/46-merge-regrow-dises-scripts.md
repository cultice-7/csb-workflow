# Regrow–DISES Merge Script

**`join_regrow_dises.py`**: Implements the spatial overlay joining Regrow fields to DISES farmer holdings, and scores the quality of each match. Runs once per state in `states`.

## Spatial merge

Both datasets are reprojected to the target equal-area CRS, and the original (unbuffered) geometries are set aside before matching. A copy of each Regrow field is buffered by `buffer_margin` (a small negative margin, to avoid spurious matches from shared edges) and intersected against every DISES farmer multipolygon. Where a Regrow field's buffered geometry overlaps more than one DISES holding, only the **largest-overlap** holding is kept per Regrow field (`field_id`). The result is left-joined back onto the full Regrow table, so unmatched fields are retained with nulls, and `parcel_assigned_dises` (`Y`/`N`) records whether any match was found. `overlap_area_dises` (acres) and `overlap_area_share_dises` (share of the Regrow field's own area) are then recomputed from the *original*, unbuffered geometries — the buffered copies were only a matching device.

## Match-quality scoring

Once a Regrow field has a DISES match, two independent comparisons are made:

- **Crop match**: DISES's 2023 crop answer is remapped onto Regrow's crop coding, then checked against either of the Regrow field's first two 2023 cultivation-cycle crops (`crop_23_1`, `crop_23_2`). `crop_match_dises` is `1` if either matches, `0` otherwise, missing if either side is unavailable.
- **Area match**: `area_match_dises` is `1` if the Regrow field's area falls within `area_match_coefs` × the DISES-reported field size (currently 0.75×–1.25×), `0` otherwise, missing if either side is unavailable.

These combine into `match_quality_dises`: `"A"` (both match), `"B_crop"` (crop only), `"B_area"` (area only), `"F"` (neither), or missing (neither comparison was possible).

## Representative-field assignment (`assign_representative_field()`)

Run once on all states' results combined, after the per-state loop finishes. Candidates are first restricted to Regrow fields that are DISES-assigned, belong to a farmer who completed the survey, and have `overlap_area_share_dises > 50%`. Each candidate is then classified into one of 4 likelihood tiers by combining `match_quality_dises` with the Regrow field's 2023 crop-confidence score (`crop_conf_col`, the minimum of its two cultivation-cycle confidence scores, thresholded at 75%):

- `RF_level_1_dises`: `match_quality_dises == "A"` and `crop_conf_col > 75`.
- `RF_level_2_dises`: `"A"` with `crop_conf_col <= 75`, or `"B_area"` with the DISES crop answer missing, or `"B_crop"` with the DISES size answer missing and `crop_conf_col > 75`.
- `RF_level_3_dises`: `"B_crop"` with size missing and `crop_conf_col <= 75`, or `"B_area"` with a crop mismatch and `crop_conf_col <= 75`.
- `RF_level_4_dises`: `"F"` with size missing, a crop mismatch, and `crop_conf_col <= 75`; or both size and crop missing with only a single DISES holding (`n_parcels_dises == 1`, an elimination case).

Within each farmer (`comp_id_dises`), only the single best candidate is kept — highest tier first, ties broken by closest area match, then by largest overlap share — and merged back onto each state's output table as `RF_level_1_dises`…`RF_level_4_dises`.

See [Merging Regrow with DISES and Supplements](../30-workflow/31-merging-regrow-dises-supplements.md) for the underlying methodology and rationale.
