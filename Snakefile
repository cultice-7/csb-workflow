configfile: "config.yml"

########### For running all rules
rule all:
    input:
        tract = "data/Census/census_tract/cb_2023_us_tract_500k.shp",
        roads = "data/Census/road/prisecroads.shp",
        elevation = "data/Geo/elevation/elevation.tif",
        slope = "data/Geo/elevation/slope.tif",
        watershed = "data/Geo/watershed/watershed.shp",
        raw_cdl = expand("data/CDL/{year}_30m_cdls/{year}_30m_cdls.tif", year=range(2014,2025)),
        validated_regrow_shape = "data/edited/Regrow/regrow_with_cdl_validation.geojson",
        regrow_joined_supplemented_2_shape = "data/edited/Regrow/regrow_dises_supplemented_2.geojson",
        csb1623_joined_supplemented_2_shape = "data/edited/CSB/CSB1623_dises_supplemented_2.geojson"



###########



#---# Download datasets (excluding DISES and Regrow which have to be downloaded manually.)

# Download census tract
rule download_census_tract:
    input:
    output:
        tract = "data/Census/census_tract/cb_2023_us_tract_500k.shp"
    params:
        raw_dir = config["census_tract"]["raw_dir"],
        output_dir = config["census_tract"]["output_dir"],
        html = f"{config["census_tract"]["base_html"]}/cb_2023_us_tract_500k.zip"
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
        html = f"{config["state_bound"]["base_html"]}/cb_2023_us_state_500k.zip"
    script:
        "scripts/download_state_bound.py"

# Download primary & secondary roads
rule download_roads:
    input:
    output:
        roads = "data/Census/road/prisecroads.shp"
    params:
        prefixes = config["roads"]["prefixes"],
        raw_dir = config["roads"]["raw_dir"],
        output_dir = config["roads"]["output_dir"],
        html = config["roads"]["base_html"]
    script:
        "scripts/download_roads.py"

# Download 3DEP elevation file
rule download_3dep:
    input:
    output:
        elevation = "data/Geo/elevation/elevation.tif"
    params:
        y_range = config["elevation"]["y_range"],
        x_range = config["elevation"]["x_range"],
        raw_dir = config["elevation"]["raw_dir"],
        output_dir = config["elevation"]["output_dir"],
        html = config["elevation"]["base_html"]
    script:
        "scripts/download_3dep.py"

# Calculate slope from 3DEP elevation
rule calculate_slope:
    input: 
        elevation = "data/Geo/elevation/elevation.tif"
    output: 
        slope = "data/Geo/elevation/slope.tif"
    script:
        "scripts/calculate_slope.py"

# Download watershed dataset (National Hydrology Dataset)
rule download_watershed:
    input:
    output:
        watershed = "data/Geo/watershed/watershed.shp"
    params:
        states = config["watershed"]["states"],
        raw_dir = config["watershed"]["raw_dir"],
        output_dir = config["watershed"]["output_dir"],
        html = config["watershed"]["base_html"]
    script:
        "scripts/download_watershed.py"

# Download Crop Sequence Boundary (CSB) dataset
rule download_csb_1623:
    input: 
    output:
        raw_csb = directory("data/CSB/CSB1623.gdb"),
        metadata = "data/CSB/NationalCSB_2016-2023_rev23_metadata.htm"
    params:
        raw_dir = config["csb"]["raw_dir"],
        output_dir = config["csb"]["output_dir"],
        html = f"{config["csb"]["base_html"]}/NationalCSB_2016-2023_rev23.zip"
    script:
        "scripts/download_csb.py"

rule download_csb_1724:
    input: 
    output:
        raw_csb = directory("data/CSB/CSB1724.gdb"),
        matadata = "data/CSB/NationalCSB_2017-2024_rev23_metadata.htm"
    params:
        raw_dir = config["csb"]["raw_dir"],
        output_dir = config["csb"]["output_dir"],
        html = f"{config["csb"]["base_html"]}/NationalCSB_2017-2024_rev23.zip"
    script:
        "scripts/download_csb.py"

# Download Cropland Data Layer (CDL) dataset
rule download_cdl:
    input:
    output:
        raw_cdl = expand("data/CDL/{year}_30m_cdls/{year}_30m_cdls.tif", year=range(2014,2025))
    params:
        year_range = config["cdl"]["year"],
        raw_dir = config["cdl"]["raw_dir"],
        output_dir = config["cdl"]["output_dir"],
        html = config["cdl"]["base_html"]
    script:
        "scripts/download_cdl.py"



#---# Prepare main datasets

# Clean DISES data table    
rule clean_dises_table:
    input:
        dises_table = "data/DISES/combined_data_clean.csv"
    output:
        dises_table_short = "data/edited/DISES/combined_data_clean_short.csv"
    script:
        "scripts/clean_dises_table.py"

# Clean DISES tax parcels and join with data table
rule clean_dises_shape:
    input:
        dises_shape = "data/DISES/DISES_All_Parcels_05.15.25.shp",
        dises_table_short = "data/edited/DISES/combined_data_clean_short.csv"
    output:
        dises_shape_consolidated = "data/edited/DISES/dises_consolidated.gpkg"
    script:
        "scripts/clean_dises_shape.py"

# Clean Regrow data table
rule clean_regrow_table:
    input:
        regrow_main_crop_table = "data/Regrow/OH_main_crop_june24.csv"
    output:
        regrow_main_crop_table_wide = "data/edited/Regrow/OH_main_crop_wide_coded.csv"
    script:
        "scripts/clean_regrow_table.py"

# Clean Regrow field geojson and join with data table
rule clean_regrow_shape:
    input:
        regrow_boundary = "data/Regrow/OSU_field_boudaries.geojson",
        regrow_main_crop_table_wide = "data/edited/Regrow/OH_main_crop_wide_coded.csv"
    output:
        regrow_shape_clean = "data/edited/Regrow/regrow_clean.geojson"
    script:
        "scripts/clean_regrow_shape.py"

# Clip CDL raster for validation with Regrow
rule clip_cdl_rasters:
    input:
        cdl = expand("data/CDL/{year}_30m_cdls/{year}_30m_cdls.tif", year=range(2014,2025)),
        states = "data/Census/state_bound/cb_2023_us_state_500k.shp"
    output:
        clipped_cdl = expand("data/edited/CDL/{year}_30m_cdls_clipped.tif", year=range(2014, 2025))
    script:
        "scripts/clip_cdl_rasters.py"

# CDL validation for Regrow
rule validate_regrow_shape:
    input: 
        regrow_shape_clean = "data/edited/Regrow/regrow_clean.geojson",
        clipped_cdl_rasters = expand("data/edited/CDL/{year}_30m_cdls_clipped.tif", year=range(2014, 2025))
    output:
        validated_regrow_shape = "data/edited/Regrow/regrow_with_cdl_validation.geojson",
        summary_regrow_validation = "data/edited/Regrow/regrow_validity_summary_by_year.csv",
        summary_regrow_validation_cdl_1_5 = "data/edited/Regrow/regrow_validity_summary_cdl_1_5.csv"
    script:
        "scripts/validate_regrow_shape.py"

# Clip CSB shape
rule clip_csb_shape:
    input:
        csb1623 = "data/CSB/CSB1623.gdb",
        csb1724 = "data/CSB/CSB1724.gdb"
    output:
        csb1623_clipped = "data/edited/CSB/CSB1623_clipped.gpkg",
        csb1724_clipped = "data/edited/CSB/CSB1724_clipped.gpkg"
    script:
        "scripts/clip_csb_shape.py"



#---# Spatial join Regrow-DISES and CSB-DISES

# Join Regrow & DISES
rule join_regrow_dises:
    input:
        regrow_shape = "data/edited/Regrow/regrow_clean.geojson",
        dises_shape = "data/edited/DISES/dises_consolidated.gpkg"
    output:
        regrow_joined_shape = "data/edited/Regrow/regrow_dises_spatialjoin.geojson",
        regrow_joined_table = "data/edited/Regrow/regrow_dises_spatialjoin_table.csv"
    script:
        "scripts/join_regrow_dises.py"

# Join CSB & DISES
rule join_csb_dises:
    input:
        csb1623_clipped = "data/edited/CSB/CSB1623_clipped.gpkg",
        csb1724_clipped = "data/edited/CSB/CSB1724_clipped.gpkg"
    output:
        csb1623_joined_shape = "data/edited/CSB/CSB1623_dises_spatialjoin.geojson",
        csb1724_joined_shape = "data/edited/CSB/CSB1724_dises_spatialjoin.geojson",
        csb1623_joined_table = "data/edited/CSB/CSB1623_dises_spatialjoin_table.csv",
        csb1724_joined_table = "data/edited/CSB/CSB1724_dises_spatialjoin_table.csv"
    script:
        "scripts/join_csb_dises.py"




#------------# Add supplementary data

#---# Supplementary data 1: Slope & elevation (from 3DEP)

# Reproject slope & elevation to equal area projection (EPSG:5070)
rule reproject_supplement_1:
    input:
        elevation = "data/Geo/elevation/elevation.tif",
        slope = "data/Geo/elevation/slope.tif"
    output:
        elevation_proj = "data/Geo/elevation/elevation_reproj.tif",
        slope_proj = "data/Geo/elevation/slope_reproj.tif"
    script:
        "scripts/reproject_supplement_1.py"

# Add mean slope & mean elevation to Regrow-DISES spatial join dataset
rule join_regrow_dises_supplement_1:
    input:
        regrow_joined_shape = "data/edited/Regrow/regrow_dises_spatialjoin.geojson",
        elevation_proj = "data/Geo/elevation/elevation_reproj.tif",
        slope_proj = "data/Geo/elevation/slope_reproj.tif"
    output:
        regrow_joined_supplemented_1_shape = "data/edited/Regrow/regrow_dises_supplemented_1.geojson",
        regrow_joined_supplemented_1_table = "data/edited/Regrow/regrow_dises_supplemented_1_table.csv"   
    script:
        "scripts/join_regrow_dises_supplement_1.py"

# Add mean slope & mean elevation to CSB-DISES spatial join dataset
rule join_csb_dises_supplement_1:
    input:
        csb1623_joined_shape = "data/edited/CSB/CSB1623_dises_spatialjoin.geojson",
        elevation_proj = "data/Geo/elevation/elevation_reproj.tif",
        slope_proj = "data/Geo/elevation/slope_reproj.tif"
    output:
        csb1623_joined_supplemented_1_shape = "data/edited/CSB/CSB1623_dises_supplemented_1.geojson",
        csb1623_joined_supplemented_1_table = "data/edited/CSB/CSB1623_dises_supplemented_1_table.csv"
    script:
        "scripts/join_csb_dises_supplement_1.py"


#---# Supplementary data 2: Nearest distance to road

# Add nearest distance to road to Regrow-DISES spatial join supplemented 1,
# generate nearest point coordinates table, and generate points on road feature for future use.
rule join_regrow_dises_supplement_2:
    input:
        regrow_joined_supplemented_1_shape = "data/edited/Regrow/regrow_dises_supplemented_1.geojson",
        roads = "data/Census/road/prisecroads.shp"
    output:
        regrow_joined_supplemented_2_shape = "data/edited/Regrow/regrow_dises_supplemented_2.geojson",
        regrow_joined_supplemented_2_table = "data/edited/Regrow/regrow_dises_supplemented_2_table.csv", 
        regrow_points_on_road_table = "data/edited/Regrow/regrow_points_on_road_table.csv",
        regrow_points_on_road_shape = "data/edited/Regrow/regrow_points_on_road.geojson"
    script:
        "scripts/join_regrow_dises_supplement_2.py"

# Add nearest distance to road to CSB-DISES spatial join supplemented 1,
# generate nearest point coordinates table, and generate points on road feature for future use.
rule join_csb_dises_supplement_2:
    input:
        csb1623_joined_supplemented_1_shape = "data/edited/CSB/CSB1623_dises_supplemented_1.geojson",
        roads = "data/Census/road/prisecroads.shp"
    output:
        csb1623_joined_supplemented_2_shape = "data/edited/CSB/CSB1623_dises_supplemented_2.geojson",
        csb1623_joined_supplemented_2_table = "data/edited/CSB/CSB1623_dises_supplemented_2_table.csv",
        csb1623_points_on_road_table = "data/edited/CSB/CSB1623_points_on_road_table.csv",
        csb1623_points_on_road_shape = "data/edited/CSB/CSB1623_points_on_road.geojson"
    script:
        "scripts/join_csb_dises_supplement_2.py"


#---# Supplementary data 3: Census tract (incl. state & county fips) and watershed

