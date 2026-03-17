# =============================================================================
# SUPPLEMENTARY DATA RULES
# Enrich Regrow and CSB datasets with 8 supplementary data layers
# =============================================================================


# =============================================================================
# Supplement 1: Census tract-level data
# =============================================================================

rule join_regrow_supplement_1:
    input:
        regrow_joined_shape = expand("data/edited/Regrow/{state}_regrow_shape_table.parquet", state=STATES),
        tract_boundaries = "data/Census/census_tract/cb_2023_us_tract_500k.shp"
    output:
        regrow_supplement_1_shape = expand("data/edited/Regrow/{state}_regrow_supplement_1_spatial.geojson", state=STATES),
        regrow_supplement_1_table = expand("data/edited/Regrow/{state}_regrow_supplement_1_table.csv",      state=STATES),
    params:
        states = STATES
    log: "logs/join_regrow_supplement_1.log"
    threads: 2
    resources:
        mem_mb = 16000
    script:
        "../scripts/join_regrow_supplement_1.py"


rule join_csb_supplement_1:
    input:
        csb1724_shape = expand("data/edited/CSB/{state}_CSB1724_clipped.gpkg", state=STATES),
        tract_boundaries = "data/Census/census_tract/cb_2023_us_tract_500k.shp"
    output:
        csb1724_supplement_1_shape = expand("data/edited/CSB/{state}_CSB1724_supplement_1_spatial.geojson", state=STATES),
        csb1724_supplement_1_table = expand("data/edited/CSB/{state}_CSB1724_supplement_1_table.csv",      state=STATES),
    params:
        states = STATES
    log: "logs/join_csb_supplement_1.log"
    threads: 2
    resources:
        mem_mb = 16000
    script:
        "../scripts/join_csb_supplement_1.py"


# =============================================================================
# Supplement 2: Slope & elevation (from 3DEP)
# =============================================================================

rule join_regrow_supplement_2:
    input:
        regrow_joined_shape = expand("data/edited/Regrow/{state}_regrow_shape_table.parquet", state=STATES),
        elevation_proj = "data/Geo/elevation/elevation_reproj.tif",
        slope_proj = "data/Geo/elevation/slope_reproj.tif"
    output:
        regrow_supplement_2_table = expand("data/edited/Regrow/{state}_regrow_supplement_2_table.csv", state=STATES),
    params:
        states = STATES
    log: "logs/join_regrow_supplement_2.log"
    threads: 2
    resources:
        mem_mb = 16000
    script:
        "../scripts/join_regrow_supplement_2.py"


rule join_csb_supplement_2:
    input:
        csb1724_shape = expand("data/edited/CSB/{state}_CSB1724_clipped.gpkg", state=STATES),
        elevation_proj = "data/Geo/elevation/elevation_reproj.tif",
        slope_proj = "data/Geo/elevation/slope_reproj.tif"
    output:
        csb1724_supplement_2_table = expand("data/edited/CSB/{state}_CSB1724_supplement_2_table.csv", state=STATES),
    params:
        states = STATES
    log: "logs/join_csb_supplement_2.log"
    threads: 2
    resources:
        mem_mb = 16000
    script:
        "../scripts/join_csb_supplement_2.py"


# =============================================================================
# Supplement 3: Watershed data (NHD)
# =============================================================================

rule join_regrow_supplement_3:
    input:
        regrow_joined_shape = expand("data/edited/Regrow/{state}_regrow_shape_table.parquet", state=STATES),
        subbasin = "data/Geo/watershed/subbasin.shp",
        watershed = "data/Geo/watershed/watershed.shp",
        subwatershed = "data/Geo/watershed/subwatershed.shp"
    output:
        regrow_supplement_3_table = expand("data/edited/Regrow/{state}_regrow_supplement_3_table.csv", state=STATES),
    params:
        states = STATES
    log: "logs/join_regrow_supplement_3.log"
    threads: 2
    resources:
        mem_mb = 16000
    script:
        "../scripts/join_regrow_supplement_3.py"


rule join_csb_supplement_3:
    input:
        csb1724_shape = expand("data/edited/CSB/{state}_CSB1724_clipped.gpkg", state=STATES),
        subbasin = "data/Geo/watershed/subbasin.shp",
        watershed = "data/Geo/watershed/watershed.shp",
        subwatershed = "data/Geo/watershed/subwatershed.shp"
    output:
        csb1724_supplement_3_table = expand("data/edited/CSB/{state}_CSB1724_supplement_3_table.csv", state=STATES),
    params:
        states = STATES
    log: "logs/join_csb_supplement_3.log"
    threads: 2
    resources:
        mem_mb = 16000
    script:
        "../scripts/join_csb_supplement_3.py"


# =============================================================================
# Supplement 4: Nearest distance to road (TIGER)
# =============================================================================

rule join_regrow_supplement_4:
    input:
        regrow_joined_shape = expand("data/edited/Regrow/{state}_regrow_shape_table.parquet", state=STATES),
        roads = "data/Census/road/prisecroads.shp"
    output:
        regrow_supplement_4_table       = expand("data/edited/Regrow/{state}_regrow_supplement_4_table.csv",       state=STATES),
        regrow_points_on_road_shape     = expand("data/edited/road/{state}_regrow_points_on_road.geojson",         state=STATES),
        regrow_points_on_road_table     = expand("data/edited/road/{state}_regrow_points_on_road_table.csv",       state=STATES),
    params:
        states = STATES
    log: "logs/join_regrow_supplement_4.log"
    threads: 2
    resources:
        mem_mb = 16000
    script:
        "../scripts/join_regrow_supplement_4.py"


rule join_csb_supplement_4:
    input:
        csb1724_shape = expand("data/edited/CSB/{state}_CSB1724_clipped.gpkg", state=STATES),
        roads = "data/Census/road/prisecroads.shp"
    output:
        csb1724_supplement_4_table       = expand("data/edited/CSB/{state}_CSB1724_supplement_4_table.csv",       state=STATES),
        csb1724_points_on_road_shape     = expand("data/edited/road/{state}_CSB1724_points_on_road.geojson",       state=STATES),
        csb1724_points_on_road_table     = expand("data/edited/road/{state}_CSB1724_points_on_road_table.csv",     state=STATES),
    params:
        states = STATES
    log: "logs/join_csb_supplement_4.log"
    threads: 2
    resources:
        mem_mb = 16000
    script:
        "../scripts/join_csb_supplement_4.py"


# =============================================================================
# Supplement 5: Neighboring field activities
# =============================================================================

rule join_regrow_supplement_5:
    input:
        regrow_joined_shape = expand("data/edited/Regrow/{state}_regrow_shape_table.parquet", state=STATES),
    output:
        regrow_supplement_5_table = expand("data/edited/Regrow/{state}_regrow_supplement_5_table.csv", state=STATES),
    params:
        states = STATES
    log: "logs/join_regrow_supplement_5.log"
    threads: 2
    resources:
        mem_mb = 16000
    script:
        "../scripts/join_regrow_supplement_5.py"


rule join_csb_supplement_5:
    input:
        csb1724_shape = expand("data/edited/CSB/{state}_CSB1724_clipped.gpkg", state=STATES),
    output:
        csb1724_supplement_5_table = expand("data/edited/CSB/{state}_CSB1724_supplement_5_table.csv", state=STATES),
    params:
        states = STATES
    log: "logs/join_csb_supplement_5.log"
    threads: 2
    resources:
        mem_mb = 16000
    script:
        "../scripts/join_csb_supplement_5.py"


# =============================================================================
# Supplement 6: Weather data (PRISM)
# =============================================================================

rule join_regrow_supplement_6:
    input:
        regrow_joined_shape = expand("data/edited/Regrow/{state}_regrow_shape_table.parquet", state=STATES),
        ppt_clipped      = expand("data/edited/Weather/ppt/{state}_prism_ppt_us_30s_{year}{month}_clipped.tif",         year=YEARS, month=MONTHS, state=STATES),
        tmax_clipped     = expand("data/edited/Weather/tmax/{state}_prism_tmax_us_30s_{year}{month}_clipped.tif",        year=YEARS, month=MONTHS, state=STATES),
        tmin_clipped     = expand("data/edited/Weather/tmin/{state}_prism_tmin_us_30s_{year}{month}_clipped.tif",        year=YEARS, month=MONTHS, state=STATES),
        tmean_clipped    = expand("data/edited/Weather/tmean/{state}_prism_tmean_us_30s_{year}{month}_clipped.tif",      year=YEARS, month=MONTHS, state=STATES),
        tdmean_clipped   = expand("data/edited/Weather/tdmean/{state}_prism_tdmean_us_30s_{year}{month}_clipped.tif",    year=YEARS, month=MONTHS, state=STATES),
        vpdmax_clipped   = expand("data/edited/Weather/vpdmax/{state}_prism_vpdmax_us_30s_{year}{month}_clipped.tif",    year=YEARS, month=MONTHS, state=STATES),
        vpdmin_clipped   = expand("data/edited/Weather/vpdmin/{state}_prism_vpdmin_us_30s_{year}{month}_clipped.tif",    year=YEARS, month=MONTHS, state=STATES),
    output:
        regrow_supplement_6_table = expand("data/edited/Regrow/{state}_regrow_supplement_6_table.csv", state=STATES),
    params:
        states = STATES,
        weather_variables = WEATHER_VARS
    log: "logs/join_regrow_supplement_6.log"
    threads: 4
    resources:
        mem_mb = 16000
    script:
        "../scripts/join_regrow_supplement_6.py"


rule join_csb_supplement_6:
    input:
        csb1724_shape = expand("data/edited/CSB/{state}_CSB1724_clipped.gpkg", state=STATES),
        ppt_clipped      = expand("data/edited/Weather/ppt/{state}_prism_ppt_us_30s_{year}{month}_clipped.tif",         year=YEARS, month=MONTHS, state=STATES),
        tmax_clipped     = expand("data/edited/Weather/tmax/{state}_prism_tmax_us_30s_{year}{month}_clipped.tif",        year=YEARS, month=MONTHS, state=STATES),
        tmin_clipped     = expand("data/edited/Weather/tmin/{state}_prism_tmin_us_30s_{year}{month}_clipped.tif",        year=YEARS, month=MONTHS, state=STATES),
        tmean_clipped    = expand("data/edited/Weather/tmean/{state}_prism_tmean_us_30s_{year}{month}_clipped.tif",      year=YEARS, month=MONTHS, state=STATES),
        tdmean_clipped   = expand("data/edited/Weather/tdmean/{state}_prism_tdmean_us_30s_{year}{month}_clipped.tif",    year=YEARS, month=MONTHS, state=STATES),
        vpdmax_clipped   = expand("data/edited/Weather/vpdmax/{state}_prism_vpdmax_us_30s_{year}{month}_clipped.tif",    year=YEARS, month=MONTHS, state=STATES),
        vpdmin_clipped   = expand("data/edited/Weather/vpdmin/{state}_prism_vpdmin_us_30s_{year}{month}_clipped.tif",    year=YEARS, month=MONTHS, state=STATES),
    output:
        csb1724_supplement_6_table = expand("data/edited/CSB/{state}_CSB1724_supplement_6_table.csv", state=STATES),
    params:
        states = STATES,
        weather_variables = WEATHER_VARS
    log: "logs/join_csb_supplement_6.log"
    threads: 4
    resources:
        mem_mb = 16000
    script:
        "../scripts/join_csb_supplement_6.py"


# =============================================================================
# Supplement 7: Crop price data (grain elevators & county indices)
# =============================================================================

rule join_regrow_supplement_7:
    input:
        regrow_joined_shape          = expand("data/edited/Regrow/{state}_regrow_supplement_1_spatial.geojson", state=STATES),
        corn_price                   = expand("data/Grain Price/{state}_corn_{level}_level.xlsx",    state=STATES, level=["elevator", "county"]),
        soybean_price                = expand("data/Grain Price/{state}_soybeans_{level}_level.xlsx", state=STATES, level=["elevator", "county"]),
        wheat_price                  = expand("data/Grain Price/{state}_wheat_elevator_level.xlsx",   state=STATES),
        elevator_location_geojson    = expand("data/edited/Grain Price/{crop}_elevator_location.geojson",       crop=["corn", "soybeans", "wheat"]),
        index_county_location_geojson = expand("data/edited/Grain Price/{crop}_index_county_location.geojson",  crop=["corn", "soybeans"]),
    output:
        regrow_supplement_7_table = expand("data/edited/Regrow/{state}_regrow_supplement_7_table.csv", state=STATES),
    params:
        states = config["price"]["states"],
        crops = config["price"]["crops"],
        number_of_neighbors = config["price"]["N_nearest"]
    log: "logs/join_regrow_supplement_7.log"
    threads: 2
    resources:
        mem_mb = 16000
    script:
        "../scripts/join_regrow_supplement_7.py"


rule join_csb_supplement_7:
    input:
        csb1724_shape                 = expand("data/edited/CSB/{state}_CSB1724_supplement_1_spatial.geojson",    state=STATES),
        corn_price                    = expand("data/Grain Price/{state}_corn_{level}_level.xlsx",    state=STATES, level=["elevator", "county"]),
        soybean_price                 = expand("data/Grain Price/{state}_soybeans_{level}_level.xlsx", state=STATES, level=["elevator", "county"]),
        wheat_price                   = expand("data/Grain Price/{state}_wheat_elevator_level.xlsx",   state=STATES),
        elevator_location_geojson     = expand("data/edited/Grain Price/{crop}_elevator_location.geojson",       crop=["corn", "soybeans", "wheat"]),
        index_county_location_geojson = expand("data/edited/Grain Price/{crop}_index_county_location.geojson",   crop=["corn", "soybeans"]),
    output:
        csb1724_supplement_7_table = expand("data/edited/CSB/{state}_CSB1724_supplement_7_table.csv", state=STATES),
    params:
        states = config["price"]["states"],
        crops = config["price"]["crops"],
        number_of_neighbors = config["price"]["N_nearest"]
    log: "logs/join_csb_supplement_7.log"
    threads: 2
    resources:
        mem_mb = 16000
    script:
        "../scripts/join_csb_supplement_7.py"


# =============================================================================
# Supplement 8: Soil composition data (gSSURGO)
# =============================================================================

rule join_regrow_supplement_8:
    input:
        regrow_geometry = expand("data/edited/Regrow/{state}_regrow_fieldID_geometry.parquet", state=STATES),
        mukey_clipped   = expand("data/edited/Soil/gSSURGO Mukey Grid/{state}_MURASTER_30m.tif", state=STATES),
        gSSURGO_tabular = expand("data/gSSURGO/gSSURGO_{state}/gSSURGO_{state}.gdb",            state=STATES),
    output:
        regrow_supplement_8_table = expand("data/edited/Regrow/{state}_regrow_supplement_8_table.parquet", state=STATES),
    params:
        states = config["soil"]["states"],
        soil_depth_cm = config["soil"]["soil_depth_cm"]
    log: "logs/join_regrow_supplement_8.log"
    threads: 2
    resources:
        mem_mb = 16000
    script:
        "../scripts/join_regrow_supplement_8.py"


rule join_csb_supplement_8:
    input:
        csb1724_geometry = expand("data/edited/CSB/{state}_CSB1724_CSBID_geometry.parquet",      state=STATES),
        mukey_clipped    = expand("data/edited/Soil/gSSURGO Mukey Grid/{state}_MURASTER_30m.tif", state=STATES),
        gSSURGO_tabular  = expand("data/gSSURGO/gSSURGO_{state}/gSSURGO_{state}.gdb",            state=STATES),
    output:
        csb1724_supplement_8_table = expand("data/edited/CSB/{state}_CSB1724_supplement_8_table.parquet", state=STATES),
    params:
        states = config["soil"]["states"],
        soil_depth_cm = config["soil"]["soil_depth_cm"]
    log: "logs/join_csb_supplement_8.log"
    threads: 2
    resources:
        mem_mb = 16000
    script:
        "../scripts/join_csb_supplement_8.py"


# =============================================================================
# Utility: Convert CSV supplement tables to Parquet (manual rule)
# =============================================================================
# This rule is not connected to the DAG and must be run manually:
#   snakemake convert_csv_to_parquet --cores 1
# It converts CSB supplement CSV tables to Parquet format for efficiency.

rule convert_csv_to_parquet:
    log: "logs/convert_csv_to_parquet.log"
    script:
        "../scripts/convert_csv_to_parquet.py"
