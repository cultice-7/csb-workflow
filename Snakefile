configfile: "config.yml"

# =============================================================================
# Gloval variables
# =============================================================================

STATES = config["global_vars"]["states"]
STATES_CODES = config["global_vars"]["states_code"]
STATES_full_name = config["global_vars"]["states_full_name"]
STATES_DISES = config["global_vars"]["states_DISES"]
YEARS = range(config["global_vars"]["years"][0], config["global_vars"]["years"][1] + 1)
CSB_YEARS = config["global_vars"]["CSB_years"]
MONTHS = [f"{m:02d}" for m in range(1, 13)] # 01-12


# =============================================================================
# All output files
# =============================================================================
rule all:
    input:
        # Regrow field_ID and geometry
        expand("data/edited/Regrow/{state}_regrow_fieldID_geometry.parquet", state=STATES),
        # CSB1724 CSBID and geometry
        expand("data/edited/CSB/{state}_CSB{years}_CSBID_geometry.parquet", state=STATES, years=CSB_YEARS),
        # Regrow tabular data on land management practices
        expand("data/edited/Regrow/{state}_regrow_table.parquet", state=STATES),
        # CSB tabular data on land management practices
        expand("data/edited/CSB/{state}_CSB{years}_table.parquet", state=STATES, years=CSB_YEARS),
        # Regrow-DISES spatial join
        expand("data/edited/Regrow/{state}_regrow_dises_table.parquet", state=STATES_DISES),
        # CSB-DISES spatial join
        expand("data/edited/CSB/{state}_CSB{years}_dises_table.parquet", state=STATES_DISES, years=CSB_YEARS),
        # Regrow supplement tables 1-8
        expand("data/edited/Regrow/{state}_regrow_supplement_{n}_table.parquet", state=STATES, n=range(1, 9)),
        # CSB supplement tables 1-8
        expand("data/edited/CSB/{state}_CSB{years}_supplement_{n}_table.parquet", state=STATES, years=CSB_YEARS, n=range(1, 9))


# =============================================================================
# Download datasets 
# Excluding Regrow, DISES, gSSURGO and elevator price which have to be downloaded manually
# =============================================================================

# Download census tract
rule download_census_tract:
    input:
    output:
        tract = "data/Census/census_tract/cb_2023_us_tract_500k.shp"
    params:
        raw_dir = config["census_tract"]["raw_dir"],
        output_dir = config["census_tract"]["output_dir"],
        html = config["census_tract"]["base_html"]
    script:
        "scripts/download_census_tract.py"

# Download state boundary
rule download_state_bound:
    input:
    output:
        state_bound = "data/Census/state_bound/cb_2023_us_state_500k.shp"
    params:
        raw_dir = config["state_bound"]["raw_dir"],
        output_dir = config["state_bound"]["output_dir"],
        html = config["state_bound"]["base_html"]
    script:
        "scripts/download_state_bound.py"

# Download county boundary
rule download_county_bound:
    input:
    output:
        county_bound = "data/Census/county_bound/cb_2023_us_county_500k.shp"
    params:
        raw_dir = config["county_bound"]["raw_dir"],
        output_dir = config["county_bound"]["output_dir"],
        html = config["county_bound"]["base_html"]
    script:
        "scripts/download_county_bound.py"

# Download primary & secondary roads
rule download_roads:
    input:
    output:
        roads = "data/Roads/prisecroads.shp"
    params:
        state_codes = config["roads"]["state_codes"],
        raw_dir = config["roads"]["raw_dir"],
        output_dir = config["roads"]["output_dir"],
        html = config["roads"]["base_html"]
    script:
        "scripts/download_roads.py"

# Download 3DEP elevation file
rule download_elevation:
    input:
        state_bound = "data/Census/state_bound/cb_2023_us_state_500k.shp"
    output:
        elevation = "data/Geo/elevation/elevation.tif"
    params:
        states = STATES,
        y_range = config["elevation"]["y_range"],
        x_range = config["elevation"]["x_range"],
        raw_dir = config["elevation"]["raw_dir"],
        output_dir = config["elevation"]["output_dir"],
        html = config["elevation"]["base_html"]
    script:
        "scripts/download_elevation.py"

# Download watershed dataset (National Hydrology Dataset)
rule download_watershed:
    input:
    output:
        subbasin = "data/Geo/watershed/subbasin.shp",
        watershed = "data/Geo/watershed/watershed.shp",
        subwatershed = "data/Geo/watershed/subwatershed.shp"
    params:
        states = STATES_full_name,
        raw_dir = config["watershed"]["raw_dir"],
        output_dir = config["watershed"]["output_dir"],
        html = config["watershed"]["base_html"]
    script:
        "scripts/download_watershed.py"

# Download weather variables (PRISM, Oregon State University)
rule download_weather:
    input:
    output:
        weather_vars = expand("data/Weather/{weather_var}/prism_{weather_var}_us_30s_{year}{month}.tif", weather_var = config["weather"]["weather_variables"], year=YEARS, month=MONTHS)
    params:
        weather_variables = config["weather"]["weather_variables"],
        year_range = YEARS,
        raw_dir = config["weather"]["raw_dir"],
        output_dir = config["weather"]["output_dir"],
        html = config["weather"]["base_html"]
    script:
        "scripts/download_weather.py"

# Download Crop Sequence Boundary (CSB) dataset
rule download_csb:
    input: 
    output:
        csb_raw = [directory(f"data/CSB/CSB{CSB_year}.gdb") for CSB_year in CSB_YEARS],
        metadata = expand("data/CSB/NationalCSB_20{year1}-20{year2}_rev23_metadata.htm", year1=[y[:2] for y in CSB_YEARS], year2=[y[2:] for y in CSB_YEARS])
    params:
        raw_dir = config["csb"]["raw_dir"],
        output_dir = config["csb"]["output_dir"],
        base_html = config["csb"]["base_html"],
        CSB_years = CSB_YEARS
    script:
        "scripts/download_csb.py"

# Download Cropland Data Layer (CDL) dataset
rule download_cdl:
    input:
    output:
        raw_cdl = expand("data/CDL/{year}_30m_cdls/{year}_30m_cdls.tif", year=YEARS)
    params:
        years_range = YEARS,
        raw_dir = config["cdl"]["raw_dir"],
        output_dir = config["cdl"]["output_dir"],
        html = config["cdl"]["base_html"]
    script:
        "scripts/download_cdl.py"


# =============================================================================
# Data processing: cleaning, clipping, reprojecting, modificating
# =============================================================================

# Reproject elevation to equal area projection (EPSG:5070)
rule reproject_elevation:
    input:
        elevation = "data/Geo/elevation/elevation.tif"
    output:
        elevation_reproj = "data/Geo/elevation/elevation_reprojected.tif"
    params:
        resampling_method = config["elevation"]["resampling_method"],
        target_CRS = config["global_vars"]["target_CRS"]
    script:
        "scripts/reproject_elevation.py"

# Calculate slope from 3DEP elevation
rule calculate_slope:
    input: 
        elevation_reproj = "data/Geo/elevation/elevation_reprojected.tif"
    output: 
        slope_reproj = "data/Geo/slope/slope_reprojected.tif"
    script:
        "scripts/calculate_slope.py"

# Clip elevation and slope raster files
rule clip_elevation_slope:
    input:
        elevation_reproj = "data/Geo/elevation/elevation_reprojected.tif",
        slope_reproj = "data/Geo/slope/slope_reprojected.tif"
    output:
        elevation_clipped = expand("data/Geo/elevation/{state}_elevation_clipped.tif", state=STATES),
        slope_clipped = expand("data/Geo/slope/{state}_slope_clipped.tif", state=STATES)
    params:
        state_bound_dir = config["state_bound"]["output_dir"],
        elevation_input_dir = config["elevation"]["output_dir"],
        slope_input_dir = config["slope"]["output_dir"],
        elevation_output_dir = config["elevation"]["output_dir"],
        slope_output_dir = config["slope"]["output_dir"],
        states = STATES,
        target_CRS = config["global_vars"]["target_CRS"]
    script:
        "scripts/clip_elevation_slope.py"

# Clip and reproject weather raster files
rule clip_reproject_weather_rasters:
    input:
        weather_vars = expand("data/Weather/{weather_var}/prism_{weather_var}_us_30s_{year}{month}.tif", weather_var = config["weather"]["weather_variables"], year=YEARS, month=MONTHS),
        states = "data/Census/state_bound/cb_2023_us_state_500k.shp"
    output:
        weather_vars_clipped = expand("data/edited/Weather/{weather_var}/{state}_prism_{weather_var}_us_30s_{year}{month}_clipped.tif", weather_var = config["weather"]["weather_variables"], year=YEARS, month=MONTHS, state=STATES)
    params:
        state_bound_dir = config["state_bound"]["output_dir"],
        weather_input_dir = config["weather"]["output_dir"],
        weather_output_dir = config["weather"]["edited_output_dir"],
        states = STATES,
        weather_variables = config["weather"]["weather_variables"],
        target_CRS = config["global_vars"]["target_CRS"]
    script:
        "scripts/clip_reproject_weather_rasters.py"

#---# Clean price data
rule clean_grain_price:
    input:
        corn_price = expand("data/Grain Price/{state}_corn_{level}_level.xlsx", state=STATES, level=["elevator", "county"]),
        soybean_price = expand("data/Grain Price/{state}_soybeans_{level}_level.xlsx", state=STATES, level=["elevator", "county"]),
        wheat_price = expand("data/Grain Price/{state}_wheat_{level}_level.xlsx", state=STATES, level=["elevator"])
    output:
        elevator_average_price = expand("data/edited/Grain Price/{crop}_monthly_average_elevator_price.csv", crop=['corn', 'soybeans', 'wheat']),
        county_average_price = expand("data/edited/Grain Price/{crop}_monthly_average_county_price.csv", crop=['corn', 'soybeans']),
        elevator_location_geojson = expand("data/edited/Grain Price/{crop}_elevator_location.geojson", crop=['corn', 'soybeans', 'wheat']),
        index_county_location_geojson = expand("data/edited/Grain Price/{crop}_index_county_location.geojson", crop=['corn', 'soybeans'])
    params:
        states = STATES,
        crops = config["price"]["crops"],
        price_levels = config["price"]["levels"],
        drop_threshold = config["price"]["drop_threshold"],
        input_dir = config["price"]["raw_dir"],
        output_dir = config["price"]["output_dir"]
    script:
        "scripts/clean_grain_price.py"

# Clip gSSURGO mukey raster files
rule clip_gSSURGO_mukey_rasters:
    input:
        gSSURGO_mukey_grid = "data/gSSURGO/FY2026_gSSURGO_mukey_grid/MURASTER_30m.tif",
        states = "data/Census/state_bound/cb_2023_us_state_500k.shp"
    output:
        mukey_clipped = expand("data/edited/Soil/gSSURGO Mukey Grid/{state}_MURASTER_30m.tif", state=STATES)
    params:
        state_bound_dir = config["state_bound"]["output_dir"],
        soil_input_dir = config["soil"]["output_dir"],
        soil_output_dir = config["soil"]["edited_output_dir"],
        states = STATES,
        target_CRS = config["global_vars"]["target_CRS"]
    script:
        "scripts/clip_gSSURGO_mukey_rasters.py"

# Clean DISES data table    
rule clean_dises_table:
    input:
        dises_table = "data/DISES/combined_data_clean.csv"
    output:
        dises_table_cleaned = "data/edited/DISES/DISES_table_cleaned.csv"
    params:
        input_dir = config["DISES"]["raw_dir"],
        output_dir = config["DISES"]["output_dir"]
    script:
        "scripts/clean_dises_table.py"

# Clean DISES field shape file (geometry)
rule clean_dises_shape:
    input:
        dises_shape = "data/DISES/DISES_All_Parcels_11.12.25.shp",
    output:
        dises_shape_cleaned = "data/edited/DISES/DISES_shape_cleaned.parquet"
    params:
        input_dir = config["DISES"]["raw_dir"],
        output_dir = config["DISES"]["output_dir"],
        target_CRS = config["global_vars"]["target_CRS"]
    script:
        "scripts/clean_dises_shape.py"

# Join DISES shape (geometry) with tabular data
rule join_dises_shape_table:
    input:
        dises_shape_cleaned = "data/edited/DISES/DISES_shape_cleaned.parquet",
        dises_table_cleaned = "data/edited/DISES/DISES_table_cleaned.csv"
    output:
        dises_shape_table = "data/edited/DISES/DISES_shape_table.parquet"
    params:
        input_dir = config["DISES"]["output_dir"],
        output_dir = config["DISES"]["output_dir"],
        target_CRS = config["global_vars"]["target_CRS"]
    script:
        "scripts/join_dises_shape_table.py"

# Join Regrow Updates (adding 2025 data)
rule join_regrow_updates:
    input:
        regrow_raw_geometry_2014_2024 = expand("data/Regrow/2014-2024/{state}_field_boundaries.geojson", state=STATES),
        regrow_raw_geometry_2025 = expand("data/Regrow/2025/{state}_2025_boundaries.geojson", state=STATES),
        regrow_raw_table_2014_2024 = expand("data/Regrow/2014-2024/Monitor_data_{state}.csv", state=config["regrow"]["states_monitor"]),
        regrow_raw_table_2025 = expand("data/Regrow/2025/{state}_2025_Data.csv", state=STATES)
    output:
        regrow_merged_geometry = expand("data/Regrow/{state}_field_boundaries.parquet", state=STATES),
        regrow_concatenated_table = expand("data/Regrow/Monitor_data_{state}.parquet", state=config["regrow"]["states_monitor"])
    params:
        input_2014_2024_dir = config["regrow"]["raw_2014_2024_dir"],
        input_2025_dir = config["regrow"]["raw_2025_dir"],
        output_dir = config["regrow"]["raw_dir"],
        states = STATES,
        states_monitor = config["regrow"]["states_monitor"],
        target_CRS = config["global_vars"]["target_CRS"]
    script:
        "scripts/join_regrow_updates.py"

# Split Regrow Monitor data by state
rule split_regrow_monitor_by_state:
    input:
        regrow_merged_geometry = expand("data/Regrow/{state}_field_boundaries.parquet", state=STATES),
        regrow_concatenated_table = expand("data/Regrow/Monitor_data_{state}.parquet", state=config["regrow"]["states_monitor"])
    output:
        regrow_table_raw = expand("data/Regrow/{state}_Monitor_data_cleaned.parquet", state=STATES)
    params:
        input_dir = config["regrow"]["raw_dir"],
        output_dir = config["regrow"]["raw_dir"],
        states = STATES,
        states_monitor = config["regrow"]["states_monitor"]
    script:
        "scripts/split_regrow_monitor_by_state.py"

# Clean Regrow data table (Monitor data)
rule clean_regrow_table:
    input:
        regrow_table_raw = expand("data/Regrow/{state}_Monitor_data_cleaned.parquet", state=STATES)
    output:
        regrow_table_wide = expand("data/edited/Regrow/{state}_regrow_monitor_wide_coded.parquet", state=STATES)
    params:
        input_dir = config["regrow"]["raw_dir"],
        output_dir = config["regrow"]["edited_output_dir"],
        states = STATES
    script:
        "scripts/clean_regrow_table.py"

# Join Regrow geometry with attribute tables
rule join_regrow_shape_table:
    input:
        regrow_boundary = expand("data/Regrow/{state}_field_boundaries.parquet", state=STATES),
        regrow_table_wide = expand("data/edited/Regrow/{state}_regrow_monitor_wide_coded.parquet", state=STATES)
    output:
        regrow_fieldID_geometry_parquet = expand("data/edited/Regrow/{state}_regrow_fieldID_geometry.parquet", state=STATES),
        regrow_fieldID_geometry_gpkg = expand("data/edited/Regrow/{state}_regrow_fieldID_geometry.gpkg", state=STATES),
        regrow_table = expand("data/edited/Regrow/{state}_regrow_table.parquet", state=STATES)
    params:
        input_geometry_dir = config["regrow"]["raw_dir"],
        input_table_dir = config["regrow"]["edited_output_dir"],
        output_dir = config["regrow"]["edited_output_dir"],
        states = STATES,
        target_CRS = config["global_vars"]["target_CRS"]
    script:
        "scripts/join_regrow_shape_table.py"

# Clip CDL raster for validation with Regrow
rule clip_cdl_rasters:
    input:
        cdl = expand("data/CDL/{year}_30m_cdls/{year}_30m_cdls.tif", year=YEARS),
        states = "data/Census/state_bound/cb_2023_us_state_500k.shp"
    output:
        clipped_cdl = expand("data/edited/CDL/{year}_30m_cdls_clipped.tif", year=YEARS)
    params:
        state_bound_dir = config["state_bound"]["output_dir"],
        cdl_input_dir = config["cdl"]["raw_dir"],
        cdl_output_dir = config["cdl"]["edited_output_dir"],
        states_code = STATES_CODES,
        years = YEARS,
        target_CRS = config["global_vars"]["target_CRS"]
    script:
        "scripts/clip_cdl_rasters.py"

# CDL validation for Regrow
rule validate_regrow_shape:
    input: 
        regrow_shape_clean = expand("data/edited/Regrow/{state}_regrow_shape_table.parquet", state=STATES),
        clipped_cdl_rasters = expand("data/edited/CDL/{year}_30m_cdls_clipped.tif", year=YEARS)
    output:
        validated_regrow_shape = expand("data/edited/Regrow/validation/{state}_regrow_with_cdl_validation.parquet", state=STATES),
        summary_regrow_validation = expand("data/edited/Regrow/validation/{state}_regrow_validity_summary_by_year.xlsx", state=STATES),
        summary_regrow_validation_cdl_1_5 = expand("data/edited/Regrow/validation/{state}_regrow_validity_summary_by_year_cdl_1_5.xlsx", state=STATES)
    params:
        regrow_input_dir = config["regrow"]["edited_output_dir"],
        regrow_output_dir = config["regrow"]["edited_output_dir"],
        cdl_input_dir =  config["cdl"]["edited_output_dir"],
        states = STATES,
        years = YEARS
    script:
        "scripts/validate_regrow_shape.py"

# Clip CSB shape
rule clip_csb_shape:
    input:
        csb_raw = expand("data/CSB/CSB{CSB_year}.gdb", CSB_year = CSB_YEARS)
    output:
        csb_clipped = expand("data/edited/CSB/{state}_CSB{CSB_year}_clipped.parquet", state=STATES, CSB_year = CSB_YEARS)
    params:
        csb_input_dir = config["csb"]["output_dir"],
        csb_output_dir = config["csb"]["edited_output_dir"],
        states = STATES,
        states_codes = STATES_CODES,
        CSB_years = CSB_YEARS
    script:
        "scripts/clip_csb_shape.py"

# Split CSB shape
rule split_csb_shape:
    input:
        csb_clipped = expand("data/edited/CSB/{state}_CSB{CSB_year}_clipped.parquet", state=STATES, CSB_year = CSB_YEARS)
    output:
        csb_shape_table = expand("data/edited/CSB/{state}_CSB{years}_shape_table.parquet", state=STATES, years=CSB_YEARS),
        csb_CSBID_geometry_parquet = expand("data/edited/CSB/{state}_CSB{years}_CSBID_geometry.parquet", state=STATES, years=CSB_YEARS),
        csb_CSBID_geometry_gpkg = expand("data/edited/CSB/{state}_CSB{years}_CSBID_geometry.gpkg", state=STATES, years=CSB_YEARS),
        csb_table = expand("data/edited/CSB/{state}_CSB{years}_table.parquet", state=STATES, years=CSB_YEARS)
    params:
        csb_input_dir = config["csb"]["edited_output_dir"],
        csb_output_dir = config["csb"]["edited_output_dir"],
        states = STATES,
        CSB_years = CSB_YEARS,
        target_CRS = config["global_vars"]["target_CRS"]
    script:
        "scripts/split_csb_shape.py"

# Rasterize Regrow to gSSURGO mukey raster format
rule rasterize_regrow_to_gSSURGO_grid:
    input:
        mukey_clipped = expand("data/edited/Soil/gSSURGO Mukey Grid/{state}_MURASTER_30m.tif", state=STATES),
        regrow_geometry = expand("data/edited/Regrow/{state}_regrow_fieldID_geometry.parquet", state=STATES)
    output:
        regrow_raster_to_gSSURGO_grid = expand("data/edited/Regrow/{state}_regrow_raster_to_gSSURGO_grid.tif", state=STATES),
        regrow_duplicating_fields = expand("data/edited/Regrow/{state}_regrow_duplicating_fields.parquet", state=STATES)
    params:
        regrow_input_dir = config["regrow"]["edited_output_dir"],
        regrow_output_dir = config["regrow"]["edited_output_dir"],
        soil_input_dir = config["soil"]["edited_output_dir"],
        states = STATES,
        target_CRS = config["global_vars"]["target_CRS"]
    script:
        "scripts/rasterize_regrow_to_gSSURGO_grid.py"

# Rasterize CSB to gSSURGO mukey raster format
rule rasterize_CSB_to_gSSURGO_grid:
    input:
        mukey_clipped = expand("data/edited/Soil/gSSURGO Mukey Grid/{state}_MURASTER_30m.tif", state=STATES),
        csb1724_geometry = expand("data/edited/CSB/{state}_CSB{CSB_year}_CSBID_geometry.parquet", state=STATES, CSB_year = CSB_YEARS)
    output:
        CSB_raster_to_gSSURGO_grid = expand("data/edited/CSB/{state}_CSB{CSB_year}_raster_to_gSSURGO_grid.tif", state=STATES, CSB_year = CSB_YEARS),
        CSB_duplicating_fields = expand("data/edited/CSB/{state}_CSB{CSB_year}_duplicating_fields.parquet", state=STATES, CSB_year = CSB_YEARS)
    params:
        CSB_input_dir = config["csb"]["edited_output_dir"],
        CSB_output_dir = config["csb"]["edited_output_dir"],
        soil_input_dir = config["soil"]["edited_output_dir"],
        states = STATES,
        CSB_years = CSB_YEARS,
        target_CRS = config["global_vars"]["target_CRS"]
    script:
        "scripts/rasterize_CSB_to_gSSURGO_grid.py"


# =============================================================================
# Spatially joining Regrow with DISES and CSB with DISES
# =============================================================================

# Spatial loin Regrow & DISES
rule join_regrow_dises:
    input:
        regrow_geometry = expand("data/edited/Regrow/{state}_regrow_fieldID_geometry.parquet", state=STATES),
        regrow_table = expand("data/edited/Regrow/{state}_regrow_table.parquet", state=STATES),
        dises_shape_table = "data/edited/DISES/DISES_shape_table.parquet"
    output:
        #regrow_dises_joined_shape = expand("data/edited/Regrow/{state}_regrow_dises_spatial.parquet", state=STATES_DISES),
        regrow_dises_joined_table = expand("data/edited/Regrow/{state}_regrow_dises_table.parquet", state=STATES_DISES)
    params:
        regrow_input_dir = config["regrow"]["edited_output_dir"],
        DISES_input_dir = config["DISES"]["output_dir"],
        regrow_DISES_output_dir = config["regrow_DISES"]["edited_output_dir"],
        states = STATES_DISES,
        buffer_margin = config["regrow_DISES"]["buffer_margin"],
        area_match_coefs = config["regrow_DISES"]["area_match_coefs"],
        target_CRS = config["global_vars"]["target_CRS"]
    script:
        "scripts/join_regrow_dises.py"

# Spatial join CSB & DISES
rule join_csb_dises:
    input:
        csb_geometry = expand("data/edited/CSB/{state}_CSB{years}_CSBID_geometry.parquet", state=STATES, years=CSB_YEARS),
        csb_table = expand("data/edited/CSB/{state}_CSB{years}_table.parquet", state=STATES, years=CSB_YEARS),
        dises_shape_table = "data/edited/DISES/DISES_shape_table.parquet"
    output:
        #csb_dises_joined_shape = expand("data/edited/CSB/{state}_CSB{years}_dises_spatial.geojson", state=STATES_DISES, years=CSB_YEARS),
        csb_dises_joined_table = expand("data/edited/CSB/{state}_CSB{years}_dises_table.parquet", state=STATES_DISES, years=CSB_YEARS)
    params:
        CSB_input_dir = config["csb"]["edited_output_dir"],
        DISES_input_dir = config["DISES"]["output_dir"],
        CSB_DISES_output_dir = config["csb_DISES"]["edited_output_dir"],
        states = STATES_DISES,
        CSB_years = CSB_YEARS,
        buffer_margin = config["csb_DISES"]["buffer_margin"],
        area_match_coefs = config["csb_DISES"]["area_match_coefs"],
        target_CRS = config["global_vars"]["target_CRS"]
    script:
        "scripts/join_csb_dises.py"


# =============================================================================
# Adding supplementary farmland characteristics
# =============================================================================

#---# Supplementary data 1: Census tract data
rule join_regrow_supplement_1:
    input:
        regrow_geometry = expand("data/edited/Regrow/{state}_regrow_fieldID_geometry.parquet", state=STATES),
        tract_boundaries = "data/Census/census_tract/cb_2023_us_tract_500k.shp"
    output:
        #regrow_supplement_1_shape = expand("data/edited/Regrow/{state}_regrow_supplement_1_spatial.parquet", state=STATES),
        regrow_supplement_1_table = expand("data/edited/Regrow/{state}_regrow_supplement_1_table.parquet", state=STATES)
    params:
        states = STATES,
        target_CRS = config["global_vars"]["target_CRS"],
        regrow_input_dir = config["regrow"]["edited_output_dir"],
        regrow_output_dir = config["regrow_supplement"]["edited_output_dir"],
        census_tract_input_dir = config["census_tract"]["output_dir"]
    script:
        "scripts/join_regrow_supplement_1.py"

rule join_csb_supplement_1:
    input:
        csb_geometry = expand("data/edited/CSB/{state}_CSB{years}_CSBID_geometry.parquet", state=STATES, years=CSB_YEARS),
        tract_boundaries = "data/Census/census_tract/cb_2023_us_tract_500k.shp"
    output:
        #csb1724_supplement_1_shape = expand("data/edited/CSB/{state}_CSB{years}_supplement_1_spatial.parquet", state=STATES, years=CSB_YEARS),
        csb1724_supplement_1_table = expand("data/edited/CSB/{state}_CSB{years}_supplement_1_table.parquet", state=STATES, years=CSB_YEARS)
    params:
        states = STATES,
        CSB_years = CSB_YEARS,
        target_CRS = config["global_vars"]["target_CRS"],
        CSB_input_dir = config["csb"]["edited_output_dir"],
        CSB_output_dir = config["csb_supplement"]["edited_output_dir"],
        census_tract_input_dir = config["census_tract"]["output_dir"]
    script:
        "scripts/join_csb_supplement_1.py"

#---# Supplementary data 2: Slope & elevation (from 3DEP)
rule join_regrow_supplement_2:
    input:
        regrow_geometry = expand("data/edited/Regrow/{state}_regrow_fieldID_geometry.parquet", state=STATES),
        elevation_clipped = expand("data/Geo/elevation/{state}_elevation_clipped.tif", state=STATES),
        slope_clipped = expand("data/Geo/slope/{state}_slope_clipped.tif", state=STATES)
    output:
        #regrow_supplement_2_shape = expand("data/edited/Regrow/{state}_regrow_supplement_2_spatial.parquet", state=STATES),
        regrow_supplement_2_table = expand("data/edited/Regrow/{state}_regrow_supplement_2_table.parquet", state=STATES)
    params:
        states = STATES,
        regrow_input_dir = config["regrow"]["edited_output_dir"],
        regrow_output_dir = config["regrow_supplement"]["edited_output_dir"],
        elevation_input_dir = config["elevation"]["output_dir"],
        slope_input_dir = config["slope"]["output_dir"]
    script:
        "scripts/join_regrow_supplement_2.py"

rule join_csb_supplement_2:
    input:
        csb_geometry = expand("data/edited/CSB/{state}_CSB{years}_CSBID_geometry.parquet", state=STATES, years=CSB_YEARS),
        elevation_clipped = expand("data/Geo/elevation/{state}_elevation_clipped.tif", state=STATES),
        slope_clipped = expand("data/Geo/slope/{state}_slope_clipped.tif", state=STATES)
    output:
        #csb1724_supplement_2_shape = expand("data/edited/CSB/{state}_CSB{years}_supplement_2_spatial.parquet", state=STATES, years=CSB_YEARS),
        csb1724_supplement_2_table = expand("data/edited/CSB/{state}_CSB{years}_supplement_2_table.parquet", state=STATES, years=CSB_YEARS)
    params:
        states = STATES,
        CSB_years = CSB_YEARS,
        CSB_input_dir = config["csb"]["edited_output_dir"],
        CSB_output_dir = config["csb_supplement"]["edited_output_dir"],
        elevation_input_dir = config["elevation"]["output_dir"],
        slope_input_dir = config["slope"]["output_dir"]
    script:
        "scripts/join_csb_supplement_2.py"

#---# Supplementary data 3: Watershed data
rule join_regrow_supplement_3:
    input:
        regrow_geometry = expand("data/edited/Regrow/{state}_regrow_fieldID_geometry.parquet", state=STATES),
        subbasin = "data/Geo/watershed/subbasin.shp",
        watershed = "data/Geo/watershed/watershed.shp",
        subwatershed = "data/Geo/watershed/subwatershed.shp"
    output:
        #regrow_supplement_3_shape = expand("data/edited/Regrow/{state}_regrow_supplement_3_spatial.parquet", state=STATES),
        regrow_supplement_3_table = expand("data/edited/Regrow/{state}_regrow_supplement_3_table.parquet", state=STATES)
    params:
        states = STATES,
        target_CRS = config["global_vars"]["target_CRS"],
        regrow_input_dir = config["regrow"]["edited_output_dir"],
        regrow_output_dir = config["regrow_supplement"]["edited_output_dir"],
        watershed_input_dir = config["watershed"]["output_dir"]
    script:
        "scripts/join_regrow_supplement_3.py"

rule join_csb_supplement_3:
    input:
        csb_geometry = expand("data/edited/CSB/{state}_CSB{years}_CSBID_geometry.parquet", state=STATES, years=CSB_YEARS),
        subbasin = "data/Geo/watershed/subbasin.shp",
        watershed = "data/Geo/watershed/watershed.shp",
        subwatershed = "data/Geo/watershed/subwatershed.shp"
    output:
        #csb1724_supplement_3_shape = expand("data/edited/CSB/{state}_CSB{years}_supplement_3_spatial.parquet", state=STATES, years=CSB_YEARS),
        csb1724_supplement_3_table = expand("data/edited/CSB/{state}_CSB{years}_supplement_3_table.parquet", state=STATES, years=CSB_YEARS)
    params:
        states = STATES,
        CSB_years = CSB_YEARS,
        target_CRS = config["global_vars"]["target_CRS"],
        CSB_input_dir = config["csb"]["edited_output_dir"],
        CSB_output_dir = config["csb_supplement"]["edited_output_dir"],
        watershed_input_dir = config["watershed"]["output_dir"]
    script:
        "scripts/join_csb_supplement_3.py"

#---# Supplementary data 4: Nearest distance to road

# Join nearest distance to road with Regrow geospatial dataset. Generate nearest point coordinates table, and generate points on road feature for future use
rule join_regrow_supplement_4:
    input:
        regrow_geometry = expand("data/edited/Regrow/{state}_regrow_fieldID_geometry.parquet", state=STATES),
        roads = "data/Roads/prisecroads.shp"
    output:
        #regrow_supplement_4_shape = expand("data/edited/Regrow/{state}_regrow_supplement_4_spatial.parquet", state=STATES),
        regrow_supplement_4_table = expand("data/edited/Regrow/{state}_regrow_supplement_4_table.parquet", state=STATES), 
        regrow_points_on_road_shape = expand("data/edited/Roads/{state}_regrow_points_on_road.geojson", state=STATES),
        regrow_points_on_road_table = expand("data/edited/Roads/{state}_regrow_points_on_road_table.csv", state=STATES)
    params:
        states = STATES,
        target_CRS = config["global_vars"]["target_CRS"],
        regrow_input_dir = config["regrow"]["edited_output_dir"],
        regrow_output_dir = config["regrow_supplement"]["edited_output_dir"],
        roads_input_dir = config["roads"]["output_dir"],
        roads_output_dir = config["roads"]["edited_output_dir"]
    script:
        "scripts/join_regrow_supplement_4.py"

# Join nearest distance to road with CSB geospatial dataset. Generate nearest point coordinates table, and generate points on road feature for future use
rule join_csb_supplement_4:
    input:
        csb_geometry = expand("data/edited/CSB/{state}_CSB{years}_CSBID_geometry.parquet", state=STATES, years=CSB_YEARS),
        roads = "data/Roads/prisecroads.shp"
    output:
        #csb1724_supplement_4_shape = expand("data/edited/CSB/{state}_CSB{years}_supplement_4_spatial.parquet", state=STATES, years=CSB_YEARS),
        csb1724_supplement_4_table = expand("data/edited/CSB/{state}_CSB{years}_supplement_4_table.parquet", state=STATES, years=CSB_YEARS),
        csb1724_points_on_road_shape = expand("data/edited/Roads/{state}_CSB{years}_points_on_road.geojson", state=STATES, years=CSB_YEARS),
        csb1724_points_on_road_table = expand("data/edited/Roads/{state}_CSB{years}_points_on_road_table.csv", state=STATES, years=CSB_YEARS)
    params:
        states = STATES,
        CSB_years = CSB_YEARS,
        target_CRS = config["global_vars"]["target_CRS"],
        CSB_input_dir = config["csb"]["edited_output_dir"],
        CSB_output_dir = config["csb_supplement"]["edited_output_dir"],
        roads_input_dir = config["roads"]["output_dir"],
        roads_output_dir = config["roads"]["edited_output_dir"]
    script:
        "scripts/join_csb_supplement_4.py"


#---# Supplementary data 5: land management activities on neighboring fields (either intersecting or sharing a boundary)
rule join_regrow_supplement_5:
    input:
        regrow_geometry = expand("data/edited/Regrow/{state}_regrow_fieldID_geometry.parquet", state=STATES),
        regrow_table = expand("data/edited/Regrow/{state}_regrow_table.parquet", state=STATES)
    output:
        #regrow_supplement_5_shape = expand("data/edited/Regrow/{state}_regrow_supplement_5_spatial.parquet", state=STATES),
        regrow_supplement_5_table = expand("data/edited/Regrow/{state}_regrow_supplement_5_table.parquet", state=STATES)
    params:
        states = STATES,
        regrow_input_dir = config["regrow"]["edited_output_dir"],
        regrow_output_dir = config["regrow_supplement"]["edited_output_dir"]
    script:
        "scripts/join_regrow_supplement_5.py"

rule join_csb_supplement_5:
    input:
        csb_geometry = expand("data/edited/CSB/{state}_CSB{years}_CSBID_geometry.parquet", state=STATES, years=CSB_YEARS),
        csb_table = expand("data/edited/CSB/{state}_CSB{years}_table.parquet", state=STATES, years=CSB_YEARS)
    output:
        #csb1724_supplement_5_shape = expand("data/edited/CSB/{state}_CSB{years}_supplement_5_spatial.parquet", state=STATES, years=CSB_YEARS),
        csb1724_supplement_5_table = expand("data/edited/CSB/{state}_CSB{years}_supplement_5_table.parquet", state=STATES, years=CSB_YEARS)
    params:
        states = STATES,
        CSB_years = CSB_YEARS,
        CSB_input_dir = config["csb"]["edited_output_dir"],
        CSB_output_dir = config["csb_supplement"]["edited_output_dir"]
    script:
        "scripts/join_csb_supplement_5.py"

        
#---# Supplementary data 6: weather data 
rule join_regrow_supplement_6:
    input:
        regrow_geometry = expand("data/edited/Regrow/{state}_regrow_fieldID_geometry.parquet", state=STATES),
        weather_vars_clipped = expand("data/edited/Weather/{weather_var}/{state}_prism_{weather_var}_us_30s_{year}{month}_clipped.tif", weather_var = config["weather"]["weather_variables"], year=YEARS, month=MONTHS, state=STATES)
    output:
        #regrow_supplement_6_shape = expand("data/edited/Regrow/{state}_regrow_supplement_6_spatial.parquet", state=STATES),
        regrow_supplement_6_table = expand("data/edited/Regrow/{state}_regrow_supplement_6_table.parquet", state=STATES)
    params:
        states = STATES,
        weather_variables = config["weather"]["weather_variables"],
        regrow_input_dir = config["regrow"]["edited_output_dir"],
        regrow_output_dir = config["regrow_supplement"]["edited_output_dir"],
        weather_input_dir = config["weather"]["edited_output_dir"]
    script:
        "scripts/join_regrow_supplement_6.py"

rule join_csb_supplement_6:
    input:
        csb_geometry = expand("data/edited/CSB/{state}_CSB{years}_CSBID_geometry.parquet", state=STATES, years=CSB_YEARS),
        weather_vars_clipped = expand("data/edited/Weather/{weather_var}/{state}_prism_{weather_var}_us_30s_{year}{month}_clipped.tif", weather_var = config["weather"]["weather_variables"], year=YEARS, month=MONTHS, state=STATES)
    output:
        #csb1724_supplement_6_shape = expand("data/edited/CSB/{state}_CSB{years}_supplement_6_spatial.parquet", state=STATES, years=CSB_YEARS),
        csb1724_supplement_6_table = expand("data/edited/CSB/{state}_CSB{years}_supplement_6_table.parquet", state=STATES, years=CSB_YEARS)
    params:
        states = STATES,
        CSB_years = CSB_YEARS,
        weather_variables = config["weather"]["weather_variables"],
        CSB_input_dir = config["csb"]["edited_output_dir"],
        CSB_output_dir = config["csb_supplement"]["edited_output_dir"],
        weather_input_dir = config["weather"]["edited_output_dir"]
    script:
        "scripts/join_csb_supplement_6.py"

#---# Supplementary data 7: crop price data
rule join_regrow_supplement_7:
    input:
        regrow_geometry = expand("data/edited/Regrow/{state}_regrow_fieldID_geometry.parquet", state=STATES),
        regrow_supplement_1_table = expand("data/edited/Regrow/{state}_regrow_supplement_1_table.parquet", state=STATES),
        elevator_average_price = expand("data/edited/Grain Price/{crop}_monthly_average_elevator_price.csv", crop=['corn', 'soybeans', 'wheat']),
        county_average_price = expand("data/edited/Grain Price/{crop}_monthly_average_county_price.csv", crop=['corn', 'soybeans']),
        elevator_location_parquet = expand("data/edited/Grain Price/{crop}_elevator_location.geojson", crop=['corn', 'soybeans', 'wheat']),
        index_county_location_geojson = expand("data/edited/Grain Price/{crop}_index_county_location.geojson", crop=['corn', 'soybeans'])
    output:
        #regrow_supplement_7_shape = expand("data/edited/Regrow/{state}_regrow_supplement_7_spatial.geojson", state=STATES),
        regrow_supplement_7_table = expand("data/edited/Regrow/{state}_regrow_supplement_7_table.parquet", state=STATES)
    params:
        states = STATES,
        crops = config["price"]["crops"],
        price_levels = config["price"]["levels"],
        number_of_neighbors = config["price"]["N_nearest"],
        K = config["price"]["K"],
        target_CRS = config["global_vars"]["target_CRS"],
        regrow_input_dir = config["regrow"]["edited_output_dir"],
        regrow_output_dir = config["regrow_supplement"]["edited_output_dir"],
        crop_price_input_dir = config["price"]["output_dir"]
    script:
        "scripts/join_regrow_supplement_7.py"

rule join_csb_supplement_7:
    input:
        csb_geometry = expand("data/edited/CSB/{state}_CSB{years}_CSBID_geometry.parquet", state=STATES, years=CSB_YEARS),
        csb_supplement_1_table = expand("data/edited/CSB/{state}_CSB{years}_supplement_1_table.parquet", state=STATES, years=CSB_YEARS),
        elevator_average_price = expand("data/edited/Grain Price/{crop}_monthly_average_elevator_price.csv", crop=['corn', 'soybeans', 'wheat']),
        county_average_price = expand("data/edited/Grain Price/{crop}_monthly_average_county_price.csv", crop=['corn', 'soybeans']),
        elevator_location_geojson = expand("data/edited/Grain Price/{crop}_elevator_location.geojson", crop=['corn', 'soybeans', 'wheat']),
        index_county_location_geojson = expand("data/edited/Grain Price/{crop}_index_county_location.geojson", crop=['corn', 'soybeans'])
    output:
        #csb1724_supplement_7_shape = expand("data/edited/CSB/{state}_CSB{years}_supplement_7_spatial.parquet", state=STATES, years=CSB_YEARS),
        csb1724_supplement_7_table = expand("data/edited/CSB/{state}_CSB{years}_supplement_7_table.parquet", state=STATES, years=CSB_YEARS)
    params:
        states = STATES,
        CSB_years = CSB_YEARS,
        crops = config["price"]["crops"],
        price_levels = config["price"]["levels"],
        number_of_neighbors = config["price"]["N_nearest"],
        K = config["price"]["K"],
        target_CRS = config["global_vars"]["target_CRS"],
        CSB_input_dir = config["csb"]["edited_output_dir"],
        CSB_output_dir = config["csb_supplement"]["edited_output_dir"],
        crop_price_input_dir = config["price"]["output_dir"]
    script:
        "scripts/join_csb_supplement_7.py"

#---# Supplementary data 8: soil composition data
rule join_regrow_supplement_8:
    input:
        regrow_geometry = expand("data/edited/Regrow/{state}_regrow_fieldID_geometry.parquet", state=STATES),
        regrow_raster_to_gSSURGO_grid = expand("data/edited/Regrow/{state}_regrow_raster_to_gSSURGO_grid.tif", state=STATES),
        regrow_duplicating_fields = expand("data/edited/Regrow/{state}_regrow_duplicating_fields.parquet", state=STATES),
        gSSURGO_tabular = expand("data/gSSURGO/gSSURGO_{state}/gSSURGO_{state}.gdb", state=STATES)
    output:
        #regrow_supplement_8_shape = expand("data/edited/Regrow/{state}_regrow_supplement_8_spatial.parquet", state=STATES),
        regrow_supplement_8_table = expand("data/edited/Regrow/{state}_regrow_supplement_8_table.parquet", state=STATES)
    params:
        states = STATES,
        soil_depth_cm = config["soil"]["soil_depth_cm"],
        regrow_input_dir = config["regrow"]["edited_output_dir"],
        regrow_output_dir = config["regrow_supplement"]["edited_output_dir"],
        soil_input_dir = config["soil"]["output_dir"],
        soil_output_dir = config["soil"]["edited_output_dir"],
        mukey_input_dir = config["soil"]["edited_output_dir"]
    script:
        "scripts/join_regrow_supplement_8.py"

rule join_csb_supplement_8:
    input:
        csb_geometry = expand("data/edited/CSB/{state}_CSB{years}_CSBID_geometry.parquet", state=STATES, years=CSB_YEARS),
        CSB_raster_to_gSSURGO_grid = expand("data/edited/CSB/{state}_CSB1724_raster_to_gSSURGO_grid.tif", state=STATES),
        CSB_duplicating_fields = expand("data/edited/CSB/{state}_CSB1724_duplicating_fields.parquet", state=STATES),
        gSSURGO_tabular = expand("data/gSSURGO/gSSURGO_{state}/gSSURGO_{state}.gdb", state=STATES)
    output:
        #csb1724_supplement_8_shape = expand("data/edited/CSB/{state}_CSB{years}_supplement_8_spatial.parquet", state=STATES, years=CSB_YEARS),
        csb1724_supplement_8_table = expand("data/edited/CSB/{state}_CSB{years}_supplement_8_table.parquet", state=STATES, years=CSB_YEARS)
    params:
        states = STATES,
        CSB_years = CSB_YEARS,
        soil_depth_cm = config["soil"]["soil_depth_cm"],
        CSB_input_dir = config["csb"]["edited_output_dir"],
        CSB_output_dir = config["csb_supplement"]["edited_output_dir"],
        soil_input_dir = config["soil"]["output_dir"],
        soil_output_dir = config["soil"]["edited_output_dir"],
        mukey_input_dir = config["soil"]["edited_output_dir"]
    script:
        "scripts/join_csb_supplement_8.py"


# =============================================================================
# Additional rules that might be useful
# =============================================================================

rule convert_csv_to_parquet:
    script: "scripts/convert_csv_to_parquet.py"