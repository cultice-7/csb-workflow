# DISES Preparation Scripts

**`clean_dises_table.py`**: Cleans the raw combined DISES survey CSV, renames key survey fields, and computes three farmer-typology indices (productivism, conservationism, civic engagement) as row-means of specific Likert-scale survey items.

**`clean_dises_shape.py`**: Cleans the raw DISES parcel shapefile, counts the number of distinct tax parcels per farmer (`n_parcels`), then **dissolves** all of a farmer's parcels into one combined multi-part geometry — meaning all downstream Regrow/CSB-to-DISES spatial joins match against a farmer's entire landholding footprint, not individual tax parcels.

**`join_dises_shape_table.py`**: Joins the dissolved DISES geometry with the cleaned survey table on farmer ID, adding a `survey_responded` flag (`"Y"` for farmers present in the survey table, `"N"` for farmers who only have a landholding/parcel record but no survey response).
