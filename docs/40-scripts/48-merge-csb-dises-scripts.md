# CSB–DISES Merge Script

**`join_csb_dises.py`**: Implements the spatial overlay joining CSB fields to DISES farmer holdings — same methodology as `join_regrow_dises.py`, but with `CSBID` as the join key, an added outer loop over `CSB_years`, and no crop-confidence dimension in the representative-field tiering (CDL/CSB has no per-field ML confidence score analog to Regrow's). See [Merging CSB with DISES and Supplements](../30-workflow/32-merging-csb-dises-supplements.md) for the full methodology and CSB-specific differences.
