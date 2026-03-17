# =============================================================================
# CLEAN & PREPARE RULES
# Clean and transform raw datasets into analysis-ready formats
# =============================================================================


rule clean_grain_price:
    input:
        corn_price    = expand("data/Grain Price/{state}_corn_{level}_level.xlsx",
                               state=STATES, level=["elevator", "county"]),
        soybean_price = expand("data/Grain Price/{state}_soybeans_{level}_level.xlsx",
                               state=STATES, level=["elevator", "county"]),
        wheat_price   = expand("data/Grain Price/{state}_wheat_elevator_level.xlsx",
                               state=STATES),
    output:
        elevator_location_geojson      = expand("data/edited/Grain Price/{crop}_elevator_location.geojson",
                                                crop=["corn", "soybeans", "wheat"]),
        elevator_location_xlsx         = expand("data/edited/Grain Price/{crop}_elevator_location.xlsx",
                                                crop=["corn", "soybeans", "wheat"]),
        index_county_location_geojson  = expand("data/edited/Grain Price/{crop}_index_county_location.geojson",
                                                crop=["corn", "soybeans"]),
        index_county_location_xlsx     = expand("data/edited/Grain Price/{crop}_index_county_location.xlsx",
                                                crop=["corn", "soybeans"]),
    params:
        states = config["price"]["states"],
        crops = config["price"]["crops"]
    log: "logs/clean_grain_price.log"
    script:
        "../scripts/clean_grain_price.py"


rule clean_dises_table:
    input:
        dises_table = "data/DISES/combined_data_clean.csv"
    output:
        dises_table_short = "data/edited/DISES/combined_data_clean_all_columns.csv"
    log: "logs/clean_dises_table.log"
    script:
        "../scripts/clean_dises_table.py"


rule clean_dises_shape:
    input:
        dises_shape = "data/DISES/DISES_All_Parcels_11.12.25.shp",
        dises_table_short = "data/edited/DISES/combined_data_clean_all_columns.csv"
    output:
        dises_shape_consolidated = "data/edited/DISES/dises_consolidated.gpkg"
    log: "logs/clean_dises_shape.log"
    threads: 2
    resources:
        mem_mb = 8000
    script:
        "../scripts/clean_dises_shape.py"


rule clean_regrow_main_crop:
    input:
        regrow_main_crop_OH = "data/Regrow/OH_main_crop_june24.csv"
    output:
        regrow_main_crop_wide_OH = "data/edited/Regrow/OH_main_crop_wide_coded.csv"
    log: "logs/clean_regrow_main_crop.log"
    script:
        "../scripts/clean_regrow_main_crop.py"


rule clean_regrow_tillage:
    input:
        regrow_tillage_OH = "data/Regrow/OH_tillage_june24.csv"
    output:
        regrow_tillage_wide_OH = "data/edited/Regrow/OH_tillage_wide_coded.csv"
    log: "logs/clean_regrow_tillage.log"
    script:
        "../scripts/clean_regrow_tillage.py"


rule clean_regrow_cover_crop:
    input:
        regrow_cover_crop_OH = "data/Regrow/OH_green_cover_crop_july7.csv"
    output:
        regrow_cover_crop_wide_OH = "data/edited/Regrow/OH_cover_crop_wide_coded.csv"
    log: "logs/clean_regrow_cover_crop.log"
    script:
        "../scripts/clean_regrow_cover_crop.py"


rule clean_regrow_shape:
    input:
        regrow_boundary = "data/Regrow/OH_field_boundaries.geojson",
        regrow_main_crop_table_wide = "data/edited/Regrow/OH_main_crop_wide_coded.csv",
        regrow_tillage_table_wide = "data/edited/Regrow/OH_tillage_wide_coded.csv",
        regrow_cover_crop_table_wide = "data/edited/Regrow/OH_cover_crop_wide_coded.csv"
    output:
        regrow_shape_clean = "data/edited/Regrow/OH_regrow_clean.geojson"
    log: "logs/clean_regrow_shape.log"
    script:
        "../scripts/clean_regrow_shape.py"


rule clean_regrow_table:
    input:
        regrow_OH = "data/Regrow/Monitor_data_OH.csv",
        regrow_MN_WI_IA_IN_IL = "data/Regrow/Monitor_data_MN_WI_IA_IN_IL.csv",
        regrow_MI = "data/Regrow/Monitor_data_MI.csv"
    output:
        regrow_wide_OH_MN_WI_IA_IN_IL_MI = "data/edited/Regrow/OH_MN_WI_IA_IN_IL_MI_regrow_wide_coded.parquet"
    log: "logs/clean_regrow_table.log"
    script:
        "../scripts/clean_regrow_table.py"


rule join_regrow_shape_table:
    input:
        regrow_boundary = expand("data/Regrow/{state}_field_boundaries.geojson", state=STATES),
        regrow_table_wide = "data/edited/Regrow/OH_MN_WI_IA_IN_IL_MI_regrow_wide_coded.parquet"
    output:
        regrow_shape_table           = expand("data/edited/Regrow/{state}_regrow_shape_table.parquet",         state=STATES),
        regrow_fieldID_geometry_parq = expand("data/edited/Regrow/{state}_regrow_fieldID_geometry.parquet",    state=STATES),
        regrow_fieldID_geometry_gpkg = expand("data/edited/Regrow/{state}_regrow_fieldID_geometry.gpkg",       state=STATES),
        regrow_table                 = expand("data/edited/Regrow/{state}_regrow_table.parquet",               state=STATES),
    params:
        states = STATES
    log: "logs/join_regrow_shape_table.log"
    threads: 2
    resources:
        mem_mb = 16000
    script:
        "../scripts/join_regrow_shape_table.py"
