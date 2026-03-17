# =============================================================================
# SPATIAL JOIN RULES
# Validation and DISES spatial joins for Regrow and CSB datasets
# =============================================================================


rule validate_regrow_shape:
    input:
        regrow_shape_clean = expand("data/edited/Regrow/{state}_regrow_shape_table.parquet", state=STATES),
        clipped_cdl_rasters = expand("data/edited/CDL/{year}_30m_cdls_clipped.tif", year=YEARS),
    output:
        validated_regrow_shape     = expand("data/edited/Regrow/validation/{state}_regrow_with_cdl_validation.geojson",                state=STATES),
        summary_regrow_validation  = expand("data/edited/Regrow/validation/{state}_regrow_validity_summary_by_year.csv",               state=STATES),
        summary_regrow_cdl_1_5     = expand("data/edited/Regrow/validation/{state}_regrow_validity_summary_by_year_cdl_1_5.csv",       state=STATES),
    params:
        states = STATES
    log: "logs/validate_regrow_shape.log"
    threads: 4
    resources:
        mem_mb = 16000
    script:
        "../scripts/validate_regrow_shape.py"


rule join_regrow_dises:
    input:
        regrow_shape = expand("data/edited/Regrow/{state}_regrow_shape_table.parquet", state=STATES),
        dises_shape = "data/edited/DISES/dises_consolidated.gpkg"
    output:
        regrow_joined_table = expand("data/edited/Regrow/{state}_regrow_dises_table.parquet", state=DISES_STATES)
    params:
        states = DISES_STATES
    log: "logs/join_regrow_dises.log"
    threads: 2
    resources:
        mem_mb = 16000
    script:
        "../scripts/join_regrow_dises.py"


rule join_csb_dises:
    input:
        csb_clipped = expand("data/edited/CSB/{state}_CSB{years}_shape_table.parquet",
                             state=STATES, years=CSB_YEARS),
        dises_shape = "data/edited/DISES/dises_consolidated.gpkg"
    output:
        csb_joined_table = expand("data/edited/CSB/{state}_CSB{years}_dises_table.parquet",
                                  state=DISES_STATES, years=CSB_YEARS)
    params:
        states = DISES_STATES,
        years = CSB_YEARS
    log: "logs/join_csb_dises.log"
    threads: 2
    resources:
        mem_mb = 16000
    script:
        "../scripts/join_csb_dises.py"
