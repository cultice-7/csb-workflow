# CSB–DISES Merge Script

**`join_csb_dises.py`**: Implements the spatial overlay joining CSB fields to DISES farmer holdings, and scores the quality of each match. Same methodology as `join_regrow_dises.py`, but keyed on `CSBID`, and with an added outer loop over `CSB_years` (currently just `"1724"`), since CSB is versioned by year-window and Regrow is not.

## Spatial merge

Both datasets are reprojected to the target equal-area CRS, and the original (unbuffered) geometries are set aside before matching. A copy of each CSB field is buffered by `buffer_margin` (a small negative margin, to avoid spurious matches from shared edges) and intersected against every DISES farmer multipolygon. Where a CSB field's buffered geometry overlaps more than one DISES holding, only the **largest-overlap** holding is kept per CSB field (`CSBID`). The result is left-joined back onto the full CSB table, so unmatched fields are retained with nulls, and `parcel_assigned_dises` (`Y`/`N`) records whether any match was found. `overlap_area_dises` (acres) and `overlap_area_share_dises` (share of the CSB field's own area) are then recomputed from the *original*, unbuffered geometries — the buffered copies were only a matching device.

## Match-quality scoring

Once a CSB field has a DISES match, two independent comparisons are made:

- **Crop match**: DISES's 2023 crop answer is remapped onto CSB's crop coding, then checked against the field's single `CDL2023` column — CSB has one annual crop call rather than multiple within-year cultivation cycles like Regrow. `crop_match_dises` is `1` if it matches, `0` otherwise, missing if either side is unavailable.
- **Area match**: `area_match_dises` is `1` if the CSB field's `CSBACRES` falls within `area_match_coefs` × the DISES-reported field size (currently 0.75×–1.25×), `0` otherwise, missing if either side is unavailable.

These combine into `match_quality_dises`: `"A"` (both match), `"B_crop"` (crop only), `"B_area"` (area only), `"F"` (neither), or missing (neither comparison was possible).

## Representative-field assignment (`assign_representative_field()`)

Run once per year, on all states' results combined for that year, after the per-state loop finishes. Candidates are first restricted to CSB fields that are DISES-assigned, belong to a farmer who completed the survey, and have `overlap_area_share_dises > 50%`. Each candidate is then classified into one of 4 likelihood tiers based on `match_quality_dises` alone — CSB has **no crop-confidence dimension**, since CDL has no per-field ML confidence analog to Regrow's, so the tiers are coarser than the Regrow side's:

- `"Level 1"`: `match_quality_dises == "A"`.
- `"Level 2"`: `"B_area"` with the DISES crop answer missing, or `"B_crop"` with the DISES size answer missing.
- `"Level 3"`: `"B_area"` with a crop mismatch.
- `"Level 4"`: `"F"` with size missing and a crop mismatch, or both size and crop missing with only a single DISES holding (`n_parcels_dises == 1`, an elimination case).

The 4 tiers are collapsed into a single `RF_assignment_dises` column via `np.select` (first matching condition wins, so `"Level 1"` takes priority if a row somehow qualifies for more than one tier), valued `"Level 1"`…`"Level 4"` or missing. Within each farmer (`comp_id_dises`), only the single best candidate is kept — highest tier first, ties broken by closest area match, then by largest overlap share — and merged back onto each state's output table as `RF_assignment_dises`.

See [Merging CSB with DISES and Supplements](../30-workflow/32-merging-csb-dises-supplements.md) for the underlying methodology and rationale.
