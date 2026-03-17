# =============================================================================
# CLIP & REPROJECT RULES
# Clip rasters/vectors to state boundaries and reproject as needed
# =============================================================================


rule calculate_slope:
    input:
        elevation = "data/Geo/elevation/elevation.tif"
    output:
        slope = "data/Geo/elevation/slope.tif"
    log: "logs/calculate_slope.log"
    script:
        "../scripts/calculate_slope.py"


rule reproject_slope_elevation:
    input:
        elevation = "data/Geo/elevation/elevation.tif",
        slope = "data/Geo/elevation/slope.tif"
    output:
        elevation_proj = "data/Geo/elevation/elevation_reproj.tif",
        slope_proj = "data/Geo/elevation/slope_reproj.tif"
    log: "logs/reproject_slope_elevation.log"
    threads: 2
    resources:
        mem_mb = 8000
    script:
        "../scripts/reproject_slope_elevation.py"


rule clip_cdl_rasters:
    input:
        cdl = expand("data/CDL/{year}_30m_cdls/{year}_30m_cdls.tif", year=YEARS),
        states = "data/Census/state_bound/cb_2023_us_state_500k.shp"
    output:
        clipped_cdl = expand("data/edited/CDL/{year}_30m_cdls_clipped.tif", year=YEARS)
    log: "logs/clip_cdl_rasters.log"
    threads: 4
    resources:
        mem_mb = 8000
    script:
        "../scripts/clip_cdl_rasters.py"


rule clip_csb_shape:
    input:
        csb1623 = "data/CSB/CSB1623.gdb",
        csb1724 = "data/CSB/CSB1724.gdb"
    output:
        csb1623_clipped = expand("data/edited/CSB/{state}_CSB1623_clipped.gpkg", state=STATES),
        csb1724_clipped = expand("data/edited/CSB/{state}_CSB1724_clipped.gpkg", state=STATES),
    params:
        states = STATES
    log: "logs/clip_csb_shape.log"
    threads: 2
    resources:
        mem_mb = 32000
    script:
        "../scripts/clip_csb_shape.py"


rule split_csb_shape:
    input:
        csb1623_clipped = expand("data/edited/CSB/{state}_CSB1623_clipped.gpkg", state=STATES),
        csb1724_clipped = expand("data/edited/CSB/{state}_CSB1724_clipped.gpkg", state=STATES),
    output:
        csb_shape_table       = expand("data/edited/CSB/{state}_CSB{years}_shape_table.parquet",       state=STATES, years=CSB_YEARS),
        csb_CSBID_geom_parq   = expand("data/edited/CSB/{state}_CSB{years}_CSBID_geometry.parquet",    state=STATES, years=CSB_YEARS),
        csb_CSBID_geom_gpkg   = expand("data/edited/CSB/{state}_CSB{years}_CSBID_geometry.gpkg",       state=STATES, years=CSB_YEARS),
        csb_table             = expand("data/edited/CSB/{state}_CSB{years}_table.parquet",              state=STATES, years=CSB_YEARS),
    params:
        states = STATES,
        years = CSB_YEARS
    log: "logs/split_csb_shape.log"
    threads: 2
    resources:
        mem_mb = 16000
    script:
        "../scripts/split_csb_shape.py"


rule clip_reproject_weather_rasters:
    input:
        ppt      = expand("data/Weather/ppt/prism_ppt_us_30s_{year}{month}.tif",         year=YEARS, month=MONTHS),
        tmax     = expand("data/Weather/tmax/prism_tmax_us_30s_{year}{month}.tif",        year=YEARS, month=MONTHS),
        tmin     = expand("data/Weather/tmin/prism_tmin_us_30s_{year}{month}.tif",        year=YEARS, month=MONTHS),
        tmean    = expand("data/Weather/tmean/prism_tmean_us_30s_{year}{month}.tif",      year=YEARS, month=MONTHS),
        tdmean   = expand("data/Weather/tdmean/prism_tdmean_us_30s_{year}{month}.tif",    year=YEARS, month=MONTHS),
        vpdmax   = expand("data/Weather/vpdmax/prism_vpdmax_us_30s_{year}{month}.tif",    year=YEARS, month=MONTHS),
        vpdmin   = expand("data/Weather/vpdmin/prism_vpdmin_us_30s_{year}{month}.tif",    year=YEARS, month=MONTHS),
        soltotal = expand("data/Weather/soltotal/prism_soltotal_us_30s_{year}{month}.tif", year=YEARS, month=MONTHS),
        solslope = expand("data/Weather/solslope/prism_solslope_us_30s_{year}{month}.tif", year=YEARS, month=MONTHS),
        states   = "data/Census/state_bound/cb_2023_us_state_500k.shp",
    output:
        ppt_clipped      = expand("data/edited/Weather/ppt/{state}_prism_ppt_us_30s_{year}{month}_clipped.tif",         year=YEARS, month=MONTHS, state=STATES),
        tmax_clipped     = expand("data/edited/Weather/tmax/{state}_prism_tmax_us_30s_{year}{month}_clipped.tif",        year=YEARS, month=MONTHS, state=STATES),
        tmin_clipped     = expand("data/edited/Weather/tmin/{state}_prism_tmin_us_30s_{year}{month}_clipped.tif",        year=YEARS, month=MONTHS, state=STATES),
        tmean_clipped    = expand("data/edited/Weather/tmean/{state}_prism_tmean_us_30s_{year}{month}_clipped.tif",      year=YEARS, month=MONTHS, state=STATES),
        tdmean_clipped   = expand("data/edited/Weather/tdmean/{state}_prism_tdmean_us_30s_{year}{month}_clipped.tif",    year=YEARS, month=MONTHS, state=STATES),
        vpdmax_clipped   = expand("data/edited/Weather/vpdmax/{state}_prism_vpdmax_us_30s_{year}{month}_clipped.tif",    year=YEARS, month=MONTHS, state=STATES),
        vpdmin_clipped   = expand("data/edited/Weather/vpdmin/{state}_prism_vpdmin_us_30s_{year}{month}_clipped.tif",    year=YEARS, month=MONTHS, state=STATES),
        soltotal_clipped = expand("data/edited/Weather/soltotal/{state}_prism_soltotal_us_30s_{year}{month}_clipped.tif", year=YEARS, month=MONTHS, state=STATES),
        solslope_clipped = expand("data/edited/Weather/solslope/{state}_prism_solslope_us_30s_{year}{month}_clipped.tif", year=YEARS, month=MONTHS, state=STATES),
    params:
        states = STATES,
        weather_variables = WEATHER_VARS
    log: "logs/clip_reproject_weather_rasters.log"
    threads: 4
    resources:
        mem_mb = 16000
    script:
        "../scripts/clip_reproject_weather_rasters.py"


rule clip_gSSURGO_rasters:
    input:
        gSSURGO_mukey_grid = "data/gSSURGO/FY2026_gSSURGO_mukey_grid/MURASTER_30m.tif",
        states = "data/Census/state_bound/cb_2023_us_state_500k.shp"
    output:
        mukey_clipped = expand("data/edited/Soil/gSSURGO Mukey Grid/{state}_MURASTER_30m.tif", state=STATES)
    params:
        states = STATES
    log: "logs/clip_gSSURGO_rasters.log"
    threads: 2
    resources:
        mem_mb = 8000
    script:
        "../scripts/clip_gSSURGO_rasters.py"
