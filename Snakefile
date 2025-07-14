configfile: "config.yml"

# Turning config file and wildcardsinto parameters
year7 = {wildcards.year} - 7

rule download_csb:
    input: 
    output:
        raw_csb = "data/NationalCSB_{year7}-{year}.gpkg"
    params:
        html = f"{config["csb"]["base_html"]}/NationalCSB_{year7}-{year}_rev23.zip"
    script:
        "scripts/download_csb.py"

rule clean_csb:
    input:
        raw_csb = "data/NationalCSB_{year7}_{year}.gpkg"
    output:
        csb = "data/SelectedCSB_{year7}_{year}.gpkg"
    params:
        states = config["csb"]["states"]
    script:
        "scripts/clean_csb.py"

rule clean_regrow:
    input:
        raw_regrow = "inputs/"
    output:
        regrow = "data/regrow_{year}.csv"
    params:
        expand({})

rule join_regrow_csb:
    input:
        csb = "data/SelectedCSB_{year7}_{year}.gpkg"
        regrow = "data/regrow_{year}.csv"
    output:
        csb_regrow_walk = "data/crosswalks/csb_regrow_assign.parquet"
    script:
        "scripts/join_csb_regrow.py"

# Karn's snake rules: need to recheck with Brian later
rule clean_dises_table:
    input:
        dises_table = "data/DISES/combined_data_clean.csv"
    output:
        dises_table_short = "data/edited/DISES/combined_data_clean_short.csv"
    script:
        "scripts/clean_dises_table.py"

rule clean_dises_shape:
    input:
        dises_shape = "data/DISES/DISES_All_Parcels_05.15.25.shp",
        dises_table_short = "data/edited/DISES/combined_data_clean_short.csv"
    output:
        dises_shape_consolidated = "data/edited/DISES/dises_consolidated.gpkg"
    script:
        "scripts/clean_dises_shape.py"

rule clean_regrow_table:
    input:
        regrow_main_crop_table = "data/Regrow/OH_main_crop_june24.csv"
    output:
        regrow_main_crop_table_wide = "data/edited/Regrow/OH_main_crop_wide_coded.csv"
    script:
        "scripts/clean_regrow_table.py"

rule clean_regrow_shape:
    input:
        regrow_boundary = "data/Regrow/OSU_field_boudaries.geojson",
        regrow_main_crop_table_wide = "data/edited/Regrow/OH_main_crop_wide_coded.csv"
    output:
        regrow_shape_clean = "data/edited/Regrow/regrow_clean.geojson"
    script:
        "scripts/clean_regrow_shape.py"

rule clip_cdl_rasters:
    input:
        cdl = "data/CDL/{year}_30m_cdls.tif",
        states = "data/Census/stat_bound/cb_2018_us_state_500k.shp"
    output:
        clipped_cdl = "data/edited/CDL/{year}_30m_cdls_clipped.tif"
    script:
        "scripts/clip_cdl_rasters.py"

rule validate_regrow_shape:
    input: 
        regrow_shape_clean = "data/edited/Regrow/regrow_clean.geojson",
        clipped_cdl_rasters = expand("data/edited/CDL/{year}_30m_cdls_clipped.tif", year=range(2014, 2025))
    output:
        validated_regrow_shape = "data/edited/Regrow/regrow_with_cdl_validation.geojson",
        summary_regrow_validation = "data/edited/Regrow/regrow_validity_summary_by_year.csv"
        summary_regrow_validation_cdl_1_5 = "data/edited/Regrow/regrow_validity_summary_cdl_1_5.csv"
    script:
        "scripts/validate_regrow_shape.py"

rule join_regrow_dises:
    input:
        regrow_shape = "data/edited/Regrow/regrow_clean.geojson",
        dises_shape = "data/edited/DISES/dises_consolidated.gpkg"
    output:
        joined_geojson = "data/edited/Regrow/regrow_dises_spatialjoin.geojson",
        joined_table = "data/edited/Regrow/regrow_dises_spatialjoin_table.csv"
    script:
        "scripts/join_regrow_dises.py"
        