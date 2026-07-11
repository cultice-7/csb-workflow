# Known Issues

These are code-level findings from a systematic review of the pipeline. As of this update, one item remains open.

## Open

1. **Three download scripts still have a latent division-by-zero / `KeyError` risk.** `download_state_bound.py`, `download_county_bound.py`, and `download_census_tract.py` safely default `Content-Length` to `0` if missing, but then divide by it unconditionally in the file-size-tolerance check, and still index `response.headers["Connection"]` directly. Fix: guard the size check with `if response_size and ...`, and use `response.headers.get("Connection", "")`.
