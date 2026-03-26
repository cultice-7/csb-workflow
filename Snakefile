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
        validated_regrow_shape = expand("data/edited/Regrow/validation/{state}_regrow_with_cdl_validation_cut.geojson", state=["OH", "MN", "WI", "IA", "IN", "IL"]),


#---# Download datasets (excluding DISES and Regrow which have to be downloaded manually)

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

# Download county boundary
rule download_county_bound:
    input:
    output:
        county_bound = "data/Census/county_bound/cb_2023_us_county_500k.shp"
    params:
        raw_dir = config["county_bound"]["raw_dir"],
        output_dir = config["county_bound"]["output_dir"],
        html = f"{config["county_bound"]["base_html"]}/cb_2023_us_county_500k.zip"
    script:
        "scripts/download_county_bound.py"

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
        state_bound = "data/Census/state_bound/cb_2023_us_state_500k.shp"
    output:
        elevation = "data/Geo/elevation/elevation.tif"
    params:
        y_range = config["elevation"]["y_range"],
        x_range = config["elevation"]["x_range"],
        raw_dir = config["elevation"]["raw_dir"],
        output_dir = config["elevation"]["output_dir"],
        html = config["elevation"]["base_html"],
        states = config["elevation"]["states"]
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
        subbasin = "data/Geo/watershed/subbasin.shp",
        watershed = "data/Geo/watershed/watershed.shp",
        subwatershed = "data/Geo/watershed/subwatershed.shp"
    params:
        states = config["watershed"]["states"],
        raw_dir = config["watershed"]["raw_dir"],
        output_dir = config["watershed"]["output_dir"],
        html = config["watershed"]["base_html"]
    script:
        "scripts/download_watershed.py"

# Download weather variables (PRISM, Oregon State University)
rule download_weather:
    input:
    output:
        ppt = expand("data/Weather/ppt/prism_ppt_us_30s_{year}{month}.tif", year=range(2014, 2026), month=[f"{m:02d}" for m in range(1, 13)]),
        tmax =  expand("data/Weather/tmax/prism_tmax_us_30s_{year}{month}.tif", year=range(2014, 2026), month=[f"{m:02d}" for m in range(1, 13)]),
        tmin =  expand("data/Weather/tmin/prism_tmin_us_30s_{year}{month}.tif", year=range(2014, 2026), month=[f"{m:02d}" for m in range(1, 13)]),
        tmean =  expand("data/Weather/tmean/prism_tmean_us_30s_{year}{month}.tif", year=range(2014, 2026), month=[f"{m:02d}" for m in range(1, 13)]),
        tdmean =  expand("data/Weather/tdmean/prism_tdmean_us_30s_{year}{month}.tif", year=range(2014, 2026), month=[f"{m:02d}" for m in range(1, 13)]),
        vpdmax =  expand("data/Weather/vpdmax/prism_vpdmax_us_30s_{year}{month}.tif", year=range(2014, 2026), month=[f"{m:02d}" for m in range(1, 13)]),
        vpdmin = expand("data/Weather/vpdmin/prism_vpdmin_us_30s_{year}{month}.tif", year=range(2014, 2026), month=[f"{m:02d}" for m in range(1, 13)]),
        soltotal = expand("data/Weather/soltotal/prism_soltotal_us_30s_{year}{month}.tif", year=range(2014, 2026), month=[f"{m:02d}" for m in range(1, 13)]),
        solslope = expand("data/Weather/solslope/prism_solslope_us_30s_{year}{month}.tif", year=range(2014, 2026), month=[f"{m:02d}" for m in range(1, 13)])
    params:
        weather_variables = config["weather"]["weather_variables"],
        states = config["weather"]["states"],
        year_range = config["weather"]["year"],
        raw_dir = config["weather"]["raw_dir"],
        output_dir = config["weather"]["output_dir"],
        html = config["weather"]["base_html"]
    script:
        "scripts/download_weather.py"

#---# Clean price data
rule clean_grain_price:
    input:
        corn_price = expand("data/Grain Price/{state}_corn_{level}_level.xlsx", state=["OH", "IN", "MI", "IA", "IL", "WI", "MN"], level=["elevator", "county"]),
        soybean_price = expand("data/Grain Price/{state}_soybeans_{level}_level.xlsx", state=["OH", "IN", "MI", "IA", "IL", "WI", "MN"], level=["elevator", "county"]),
        wheat_price = expand("data/Grain Price/{state}_wheat_{level}_level.xlsx", state=["OH", "IN", "MI", "IA", "IL", "WI", "MN"], level=["elevator"])
    output:
        elevator_location_geojson = expand("data/edited/Grain Price/{crop}_elevator_location.geojson", crop=['corn', 'soybeans', 'wheat']),
        elevator_location_xlsx = expand("data/edited/Grain Price/{crop}_elevator_location.xlsx", crop=['corn', 'soybeans', 'wheat']),
        index_county_location_geojson = expand("data/edited/Grain Price/{crop}_index_county_location.geojson", crop=['corn', 'soybeans']),
        index_county_location_xlsx = expand("data/edited/Grain Price/{crop}_index_county_location.xlsx", crop=['corn', 'soybeans'])
    params:
        states = config["price"]["states"],
        crops = config["price"]["crops"]
    script:
        "scripts/clean_grain_price.py"

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
        dises_table_short = "data/edited/DISES/combined_data_clean_all_columns.csv"
    script:
        "scripts/clean_dises_table.py"

# Clean DISES tax parcels and join with data table
rule clean_dises_shape:
    input:
        dises_shape = "data/DISES/DISES_All_Parcels_11.12.25.shp",
        dises_table_short = "data/edited/DISES/combined_data_clean_all_columns.csv"
    output:
        dises_shape_consolidated = "data/edited/DISES/dises_consolidated.gpkg"
    script:
        "scripts/clean_dises_shape.py"

# Join Regrow Updates (adding 2025 data)
rule join_regrow_updates:
    input:
        regrow_raw_geometry_2014_2024 = expand("data/Regrow/2014-2024/{state}_field_boundaries.geojson", state=["OH", "MN", "WI", "IA", "IN", "IL", "MI"]),
        regrow_raw_geometry_2025 = expand("data/Regrow/2025/{state}_2025_boundaries.geojson", state=["OH", "MN", "WI", "IA", "IN", "IL", "MI"]),
        regrow_raw_table_2014_2024 = expand("data/Regrow/2014-2024/Monitor_data_{state}.csv", state=["OH", "MN_WI_IA_IN_IL", "MI"]),
        regrow_raw_table_2025 = expand("data/Regrow/2025/{state}_2025_Data.csv", state=["OH", "MN", "WI", "IA", "IN", "IL", "MI"])
    output:
        regrow_merged_geometry = expand("data/Regrow/{state}_field_boundaries.geojson", state=["OH", "MN", "WI", "IA", "IN", "IL", "MI"]),
        regrow_concatenated_table = expand("data/Regrow/Monitor_data_{state}.csv", state=["OH", "MN_WI_IA_IN_IL", "MI"])
    params:
        input_2014_2024_dir = config["regrow"]["raw_2014_2024_dir"],
        input_2025_dir = config["regrow"]["raw_2025_dir"],
        output_dir = config["regrow"]["raw_dir"],
        states = config["regrow"]["states"]
    script:
        "scripts/join_regrow_updates.py"

# Clean Regrow data table (monitor data)
rule clean_regrow_table:
    input:
        regrow_concatenated_table = expand("data/Regrow/Monitor_data_{state}.csv", state=["OH", "MN_WI_IA_IN_IL", "MI"])
    output:
        regrow_wide_OH_MN_WI_IA_IN_IL_MI = "data/edited/Regrow/OH_MN_WI_IA_IN_IL_MI_regrow_wide_coded.parquet"
    params:
        input_dir = config["regrow"]["raw_dir"],
        output_dir = config["regrow"]["edited_output_dir"],
        states_monitor = config["regrow"]["states_monitor"]
    script:
        "scripts/clean_regrow_table.py"

# Join Regrow geometry with attribute tables
rule join_regrow_shape_table:
    input:
        regrow_boundary = expand("data/Regrow/{state}_field_boundaries.geojson", state=["OH", "MN", "WI", "IA", "IN", "IL", "MI"]),
        regrow_table_wide_OH_MN_WI_IA_IN_IL_MI = f"data/edited/Regrow/OH_MN_WI_IA_IN_IL_MI_regrow_wide_coded.parquet"
    output:
        regrow_shape_table = expand("data/edited/Regrow/{state}_regrow_shape_table.parquet", state=["OH", "MN", "WI", "IA", "IN", "IL", "MI"]),
        regrow_fieldID_geometry_parquet = expand("data/edited/Regrow/{state}_regrow_fieldID_geometry.parquet", state=["OH", "MN", "WI", "IA", "IN", "IL", "MI"]),
        regrow_fieldID_geometry_gpkg = expand("data/edited/Regrow/{state}_regrow_fieldID_geometry.gpkg", state=["OH", "MN", "WI", "IA", "IN", "IL", "MI"]),
        regrow_table = expand("data/edited/Regrow/{state}_regrow_table.parquet", state=["OH", "MN", "WI", "IA", "IN", "IL", "MI"])
    params:
        input_geometry_dir = config["regrow"]["raw_dir"],
        input_table_dir = config["regrow"]["edited_output_dir"],
        output_dir = config["regrow"]["edited_output_dir"],
        states = ["OH", "MN", "WI", "IA", "IN", "IL", "MI"]
    script:
        "scripts/join_regrow_shape_table.py"

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
        regrow_shape_clean = expand("data/edited/Regrow/{state}_regrow_shape_table.geojson", state=["OH", "MN", "WI", "IA", "IN", "IL", "MI"]),
        clipped_cdl_rasters = expand("data/edited/CDL/{year}_30m_cdls_clipped.tif", year=range(2014, 2025))
    output:
        validated_regrow_shape = expand("data/edited/Regrow/validation/{state}_regrow_with_cdl_validation.geojson", state=["OH", "MN", "WI", "IA", "IN", "IL", "MI"]),
        summary_regrow_validation = expand("data/edited/Regrow/validation/{state}_regrow_validity_summary_by_year.csv", state=["OH", "MN", "WI", "IA", "IN", "IL", "MI"]),
        summary_regrow_validation_cdl_1_5 = expand("data/edited/Regrow/validation/{state}_regrow_validity_summary_by_year_cdl_1_5.csv", state=["OH", "MN", "WI", "IA", "IN", "IL", "MI"])
    params:
        states = ["OH", "MN", "WI", "IA", "IN", "IL", "MI"]
    script:
        "scripts/validate_regrow_shape.py"

# Clip CSB shape
rule clip_csb_shape:
    input:
        csb1623 = "data/CSB/CSB1623.gdb",
        csb1724 = "data/CSB/CSB1724.gdb"
    output:
        csb1623_clipped = expand("data/edited/CSB/{state}_CSB1623_clipped.gpkg", state=["OH", "MN", "WI", "IA", "IN", "IL", "MI"]),
        csb1724_clipped = expand("data/edited/CSB/{state}_CSB1724_clipped.gpkg", state=["OH", "MN", "WI", "IA", "IN", "IL", "MI"])
    params:
        states = ["OH", "MN", "WI", "IA", "IN", "IL", "MI"]
    script:
        "scripts/clip_csb_shape.py"

# Split CSB shape
rule split_csb_shape:
    input:
        csb1623_clipped = expand("data/edited/CSB/{state}_CSB1623_clipped.gpkg", state=["OH", "MN", "WI", "IA", "IN", "IL", "MI"]),
        csb1724_clipped = expand("data/edited/CSB/{state}_CSB1724_clipped.gpkg", state=["OH", "MN", "WI", "IA", "IN", "IL", "MI"])
    output:
        csb_shape_table = expand("data/edited/CSB/{state}_CSB{years}_shape_table.parquet", state=["OH", "MN", "WI", "IA", "IN", "IL", "MI"], years=["1623", "1724"]),
        csb_CSBID_geometry_parquet = expand("data/edited/CSB/{state}_CSB{years}_CSBID_geometry.parquet", state=["OH", "MN", "WI", "IA", "IN", "IL", "MI"], years=["1623", "1724"]),
        csb_CSBID_geometry_gpkg = expand("data/edited/CSB/{state}_CSB{years}_CSBID_geometry.gpkg", state=["OH", "MN", "WI", "IA", "IN", "IL", "MI"], years=["1623", "1724"]),
        csb_table = expand("data/edited/CSB/{state}_CSB{years}_table.parquet", state=["OH", "MN", "WI", "IA", "IN", "IL", "MI"], years=["1623", "1724"])
    params:
        states = ["OH", "MN", "WI", "IA", "IN", "IL", "MI"],
        years = ["1623", "1724"]
    script:
        "scripts/split_csb_shape.py"

# Clip and reproject weather raster files
rule clip_reproject_weather_rasters:
    input:
        ppt = expand("data/Weather/ppt/prism_ppt_us_30s_{year}{month}.tif", year=range(2014, 2026), month=[f"{m:02d}" for m in range(1, 13)]),
        tmax =  expand("data/Weather/tmax/prism_tmax_us_30s_{year}{month}.tif", year=range(2014, 2026), month=[f"{m:02d}" for m in range(1, 13)]),
        tmin =  expand("data/Weather/tmin/prism_tmin_us_30s_{year}{month}.tif", year=range(2014, 2026), month=[f"{m:02d}" for m in range(1, 13)]),
        tmean =  expand("data/Weather/tmean/prism_tmean_us_30s_{year}{month}.tif", year=range(2014, 2026), month=[f"{m:02d}" for m in range(1, 13)]),
        tdmean =  expand("data/Weather/tdmean/prism_tdmean_us_30s_{year}{month}.tif", year=range(2014, 2026), month=[f"{m:02d}" for m in range(1, 13)]),
        vpdmax =  expand("data/Weather/vpdmax/prism_vpdmax_us_30s_{year}{month}.tif", year=range(2014, 2026), month=[f"{m:02d}" for m in range(1, 13)]),
        vpdmin = expand("data/Weather/vpdmin/prism_vpdmin_us_30s_{year}{month}.tif", year=range(2014, 2026), month=[f"{m:02d}" for m in range(1, 13)]),
        soltotal = expand("data/Weather/soltotal/prism_soltotal_us_30s_{year}{month}.tif", year=range(2014, 2026), month=[f"{m:02d}" for m in range(1, 13)]),
        solslope = expand("data/Weather/solslope/prism_solslope_us_30s_{year}{month}.tif", year=range(2014, 2026), month=[f"{m:02d}" for m in range(1, 13)]),
        states = "data/Census/state_bound/cb_2023_us_state_500k.shp"
    output:
        ppt_clipped = expand("data/edited/Weather/ppt/{state}_prism_ppt_us_30s_{year}{month}_clipped.tif", year=range(2014, 2026), month=[f"{m:02d}" for m in range(1, 13)], state=["OH", "MN", "WI", "IA", "IN", "IL", "MI"]),
        tdmean_clipped =  expand("data/edited/Weather/tdmean/{state}_prism_tdmean_us_30s_{year}{month}_clipped.tif", year=range(2014, 2026), month=[f"{m:02d}" for m in range(1, 13)], state=["OH", "MN", "WI", "IA", "IN", "IL", "MI"]),
        tmax_clipped =  expand("data/edited/Weather/tmax/{state}_prism_tmax_us_30s_{year}{month}_clipped.tif", year=range(2014, 2026), month=[f"{m:02d}" for m in range(1, 13)], state=["OH", "MN", "WI", "IA", "IN", "IL", "MI"]),
        tmin_clipped =  expand("data/edited/Weather/tmin/{state}_prism_tmin_us_30s_{year}{month}_clipped.tif", year=range(2014, 2026), month=[f"{m:02d}" for m in range(1, 13)], state=["OH", "MN", "WI", "IA", "IN", "IL", "MI"]),
        tmean_clipped =  expand("data/edited/Weather/tmean/{state}_prism_tmean_us_30s_{year}{month}_clipped.tif", year=range(2014, 2026), month=[f"{m:02d}" for m in range(1, 13)], state=["OH", "MN", "WI", "IA", "IN", "IL", "MI"]),
        vpdmax_clipped =  expand("data/edited/Weather/vpdmax/{state}_prism_vpdmax_us_30s_{year}{month}_clipped.tif", year=range(2014, 2026), month=[f"{m:02d}" for m in range(1, 13)], state=["OH", "MN", "WI", "IA", "IN", "IL", "MI"]),
        vpdmin_clipped = expand("data/edited/Weather/vpdmin/{state}_prism_vpdmin_us_30s_{year}{month}_clipped.tif", year=range(2014, 2026), month=[f"{m:02d}" for m in range(1, 13)], state=["OH", "MN", "WI", "IA", "IN", "IL", "MI"]),
        soltotal_clipped = expand("data/edited/Weather/soltotal/{state}_prism_soltotal_us_30s_{year}{month}_clipped.tif", year=range(2014, 2026), month=[f"{m:02d}" for m in range(1, 13)], state=["OH", "MN", "WI", "IA", "IN", "IL", "MI"]),
        solslope_clipped = expand("data/edited/Weather/solslope/{state}_prism_solslope_us_30s_{year}{month}_clipped.tif", year=range(2014, 2026), month=[f"{m:02d}" for m in range(1, 13)], state=["OH", "MN", "WI", "IA", "IN", "IL", "MI"])
    params:
        states = config["weather"]["states"],
        weather_variables = config["weather"]["weather_variables"]
    script:
        "scripts/clip_reproject_weather_rasters.py"

# Clip gSSURGO mukey raster files
rule clip_gSSURGO_rasters:
    input:
        gSSURGO_mukey_grid = "data/gSSURGO/FY2026_gSSURGO_mukey_grid\MURASTER_30m.tif",
        states = "data/Census/state_bound/cb_2023_us_state_500k.shp"
    output:
        mukey_clipped = expand("data/edited/Soil/gSSURGO Mukey Grid/{state}_MURASTER_30m.tif", state=["OH", "MN", "WI", "IA", "IN", "IL", "MI"])
    params:
        states = ['IN', 'IL', 'MN', 'IA', 'OH', 'WI', 'MI']
    script:
        "scripts/clip_gSSURGO_rasters.py"

# Rasterize Regrow to gSSURGO mukey raster format
rule rasterize_regrow_to_gSSURGO_grid:
    input:
        mukey_clipped = expand("data/edited/Soil/gSSURGO Mukey Grid/{state}_MURASTER_30m.tif", state=["OH", "MN", "WI", "IA", "IN", "IL", "MI"]),
        regrow_geometry = expand("data/edited/Regrow/{state}_regrow_fieldID_geometry.parquet", state=["OH", "MN", "WI", "IA", "IN", "IL", "MI"])
    output:
        regrow_raster_to_gSSURGO_grid = expand("data/edited/Regrow/{state}_regrow_raster_to_gSSURGO_grid.tif", state=["OH", "MN", "WI", "IA", "IN", "IL", "MI"]),
        regrow_duplicating_fields = expand("data/edited/Regrow/{state}_regrow_duplicating_fields.parquet", state=["OH", "MN", "WI", "IA", "IN", "IL", "MI"])
    params:
        states = config["soil"]["states"],
        regrow_input_dir = config["regrow"]["edited_output_dir"],
        regrow_output_dir = config["regrow"]["edited_output_dir"],
        soil_output_dir = config["soil"]["edited_output_dir"]
    script:
        "scripts/rasterize_regrow_to_gSSURGO_grid.py"

# Rasterize CSB to gSSURGO mukey raster format
rule rasterize_CSB_to_gSSURGO_grid:
    input:
        mukey_clipped = expand("data/edited/Soil/gSSURGO Mukey Grid/{state}_MURASTER_30m.tif", state=["OH", "MN", "WI", "IA", "IN", "IL", "MI"]),
        csb1724_geometry = expand("data/edited/CSB/{state}_CSB1724_CSBID_geometry.parquet", state=['IN', 'IL', 'MN', 'IA', 'OH', 'WI', 'MI'])
    output:
        CSB_raster_to_gSSURGO_grid = expand("data/edited/CSB/{state}_CSB1724_raster_to_gSSURGO_grid.tif", state=["OH", "MN", "WI", "IA", "IN", "IL", "MI"]),
        CSB_duplicating_fields = expand("data/edited/CSB/{state}_CSB1724_duplicating_fields.parquet", state=["OH", "MN", "WI", "IA", "IN", "IL", "MI"])
    params:
        states = config["soil"]["states"],
        CSB_input_dir = config["csb"]["edited_output_dir"],
        CSB_output_dir = config["csb"]["edited_output_dir"],
        soil_output_dir = config["soil"]["edited_output_dir"]
    script:
        "scripts/rasterize_CSB_to_gSSURGO_grid.py"




#---# Spatial join Regrow-DISES and CSB-DISES

# Join Regrow & DISES
rule join_regrow_dises:
    input:
        regrow_shape = expand("data/edited/Regrow/{state}_regrow_shape_table.parquet", state=["OH", "MN", "WI", "IA", "IN", "IL", "MI"]),
        dises_shape = "data/edited/DISES/dises_consolidated.gpkg"
    output:
        #regrow_joined_shape = expand("data/edited/Regrow/{state}_regrow_dises_spatial.parquet", state=["OH", "IN", "MI"]),
        regrow_joined_table = expand("data/edited/Regrow/{state}_regrow_dises_table.parquet", state=["OH", "IN", "MI"])
    params:
        states = ["OH", "IN", "MI"]
    script:
        "scripts/join_regrow_dises.py"

# Join CSB & DISES
rule join_csb_dises:
    input:
        csb_clipped = expand("data/edited/CSB/{state}_CSB{years}_shape_table.parquet", state=["OH", "MN", "WI", "IA", "IN", "IL", "MI"], years=["1623", "1724"]),
        dises_shape = "data/edited/DISES/dises_consolidated.gpkg"
    output:
        #csb_joined_shape = expand("data/edited/CSB/{state}_CSB{years}_dises_spatial.geojson", state=["OH", "IN", "MI"], years=["1623", "1724"]),
        csb_joined_table = expand("data/edited/CSB/{state}_CSB{years}_dises_table.parquet", state=["OH", "IN", "MI"], years=["1623", "1724"])
    params:
        states = ["OH", "IN", "MI"],
        years = ["1623", "1724"]
    script:
        "scripts/join_csb_dises.py"




#------------# Add supplementary data

# Reproject slope & elevation to equal area projection (EPSG:5070)
rule reproject_slope_elevation:
    input:
        elevation = "data/Geo/elevation/elevation.tif",
        slope = "data/Geo/elevation/slope.tif"
    output:
        elevation_proj = "data/Geo/elevation/elevation_reproj.tif",
        slope_proj = "data/Geo/elevation/slope_reproj.tif"
    script:
        "scripts/reproject_slope_elevation.py"

#---# Supplementary data 1: Census tract-level data
rule join_regrow_supplement_1:
    input:
        regrow_joined_shape = expand("data/edited/Regrow/{state}_regrow_shape_table.geojson", state=['IN', 'IL', 'MN', 'IA', 'OH', 'WI', 'MI']),
        tract_boundaries = "data/Census/census_tract/cb_2023_us_tract_500k.shp"
    output:
        regrow_supplement_1_shape = expand("data/edited/Regrow/{state}_regrow_supplement_1_spatial.geojson", state=['IN', 'IL', 'MN', 'IA', 'OH', 'WI', 'MI']),
        regrow_supplement_1_table = expand("data/edited/Regrow/{state}_regrow_supplement_1_table.csv", state=['IN', 'IL', 'MN', 'IA', 'OH', 'WI', 'MI'])
    params:
        states = ['IN', 'IL', 'MN', 'IA', 'OH', 'WI', 'MI']
    script:
        "scripts/join_regrow_supplement_1.py"

rule join_csb_supplement_1:
    input:
        csb1724_shape = expand("data/edited/CSB/{state}_CSB1724_clipped.gpkg", state=['IN', 'IL', 'MN', 'IA', 'OH', 'WI', 'MI']),
        tract_boundaries = "data/Census/census_tract/cb_2023_us_tract_500k.shp"
    output:
        csb1724_supplement_1_shape = expand("data/edited/CSB/{state}_CSB1724_supplement_1_spatial.geojson", state=['IN', 'IL', 'MN', 'IA', 'OH', 'WI', 'MI']),
        csb1724_supplement_1_table = expand("data/edited/CSB/{state}_CSB1724_supplement_1_table.csv", state=['IN', 'IL', 'MN', 'IA', 'OH', 'WI', 'MI'])
    params:
        states = ['IN', 'IL', 'MN', 'IA', 'OH', 'WI', 'MI']
    script:
        "scripts/join_csb_supplement_1.py"

#---# Supplementary data 2: Slope & elevation (from 3DEP)
rule join_regrow_supplement_2:
    input:
        regrow_joined_shape = expand("data/edited/Regrow/{state}_regrow_shape_table.geojson", state=['IN', 'IL', 'MN', 'IA', 'OH', 'WI', 'MI']),
        elevation_proj = "data/Geo/elevation/elevation_reproj.tif",
        slope_proj = "data/Geo/elevation/slope_reproj.tif"
    output:
        #regrow_supplement_2_shape = expand("data/edited/Regrow/{state}_regrow_supplement_2_spatial.geojson", state=['IN', 'IL', 'MN', 'IA', 'OH', 'WI', 'MI']),
        regrow_supplement_2_table = expand("data/edited/Regrow/{state}_regrow_supplement_2_table.csv", state=['IN', 'IL', 'MN', 'IA', 'OH', 'WI', 'MI'])
    params:
        states = ['IN', 'IL', 'MN', 'IA', 'OH', 'WI', 'MI']
    script:
        "scripts/join_regrow_supplement_2.py"

rule join_csb_supplement_2:
    input:
        csb1724_shape = expand("data/edited/CSB/{state}_CSB1724_clipped.gpkg", state=['IN', 'IL', 'MN', 'IA', 'OH', 'WI', 'MI']),
        elevation_proj = "data/Geo/elevation/elevation_reproj.tif",
        slope_proj = "data/Geo/elevation/slope_reproj.tif"
    output:
        #csb1724_supplement_2_shape = expand("data/edited/CSB/{state}_CSB1724_supplement_2_spatial.geojson", state=['IN', 'IL', 'MN', 'IA', 'OH', 'WI', 'MI']),
        csb1724_supplement_2_table = expand("data/edited/CSB/{state}_CSB1724_supplement_2_table.csv", state=['IN', 'IL', 'MN', 'IA', 'OH', 'WI', 'MI'])
    params:
        states = ['IN', 'IL', 'MN', 'IA', 'OH', 'WI', 'MI']
    script:
        "scripts/join_csb_supplement_2.py"

#---# Supplementary data 3: Watershed data
rule join_regrow_supplement_3:
    input:
        regrow_joined_shape = expand("data/edited/Regrow/{state}_regrow_shape_table.geojson", state=['IN', 'IL', 'MN', 'IA', 'OH', 'WI', 'MI']),
        subbasin = "data/Geo/watershed/subbasin.shp",
        watershed = "data/Geo/watershed/watershed.shp",
        subwatershed = "data/Geo/watershed/subwatershed.shp"
    output:
        #regrow_supplement_3_shape = expand("data/edited/Regrow/{state}_regrow_supplement_3_spatial.geojson", state=['IN', 'IL', 'MN', 'IA', 'OH', 'WI', 'MI']),
        regrow_supplement_3_table = expand("data/edited/Regrow/{state}_regrow_supplement_3_table.csv", state=['IN', 'IL', 'MN', 'IA', 'OH', 'WI', 'MI'])
    params:
        states = ['IN', 'IL', 'MN', 'IA', 'OH', 'WI', 'MI']
    script:
        "scripts/join_regrow_supplement_3.py"

rule join_csb_supplement_3:
    input:
        csb1724_shape = expand("data/edited/CSB/{state}_CSB1724_clipped.gpkg", state=['IN', 'IL', 'MN', 'IA', 'OH', 'WI', 'MI']),
        subbasin = "data/Geo/watershed/subbasin.shp",
        watershed = "data/Geo/watershed/watershed.shp",
        subwatershed = "data/Geo/watershed/subwatershed.shp"
    output:
        #csb1724_supplement_3_shape = expand("data/edited/CSB/{state}_CSB1724_supplement_3_spatial.geojson", state=['IN', 'IL', 'MN', 'IA', 'OH', 'WI', 'MI']),
        csb1724_supplement_3_table = expand("data/edited/CSB/{state}_CSB1724_supplement_3_table.csv", state=['IN', 'IL', 'MN', 'IA', 'OH', 'WI', 'MI'])
    params:
        states = ['IN', 'IL', 'MN', 'IA', 'OH', 'WI', 'MI']
    script:
        "scripts/join_csb_supplement_3.py"

#---# Supplementary data 4: Nearest distance to road

# Join nearest distance to road with Regrow geospatial dataset. Generate nearest point coordinates table, and generate points on road feature for future use
rule join_regrow_supplement_4:
    input:
        regrow_joined_shape = expand("data/edited/Regrow/{state}_regrow_shape_table.geojson", state=['IN', 'IL', 'MN', 'IA', 'OH', 'WI', 'MI']),
        roads = "data/Census/road/prisecroads.shp"
    output:
        #regrow_supplement_4_shape = expand("data/edited/Regrow/{state}_regrow_supplement_4_spatial.geojson", state=['IN', 'IL', 'MN', 'IA', 'OH', 'WI', 'MI']),
        regrow_supplement_4_table = expand("data/edited/Regrow/{state}_regrow_supplement_4_table.csv", state=['IN', 'IL', 'MN', 'IA', 'OH', 'WI', 'MI']), 
        regrow_points_on_road_shape = expand("data/edited/road/{state}_regrow_points_on_road.geojson", state=['IN', 'IL', 'MN', 'IA', 'OH', 'WI', 'MI']),
        regrow_points_on_road_table = expand("data/edited/road/{state}_regrow_points_on_road_table.csv", state=['IN', 'IL', 'MN', 'IA', 'OH', 'WI', 'MI'])
    params:
        states = ['IN', 'IL', 'MN', 'IA', 'OH', 'WI', 'MI']
    script:
        "scripts/join_regrow_supplement_4.py"

# Join nearest distance to road with CSB geospatial dataset. Generate nearest point coordinates table, and generate points on road feature for future use
rule join_csb_supplement_4:
    input:
        csb1724_shape = expand("data/edited/CSB/{state}_CSB1724_clipped.gpkg", state=['IN', 'IL', 'MN', 'IA', 'OH', 'WI', 'MI']),
        roads = "data/Census/road/prisecroads.shp"
    output:
        #csb1724_supplement_4_shape = expand("data/edited/CSB/{state}_CSB1724_supplement_4_spatial.geojson", state=['IN', 'IL', 'MN', 'IA', 'OH', 'WI', 'MI']),
        csb1724_supplement_4_table = expand("data/edited/CSB/{state}_CSB1724_supplement_4_table.csv", state=['IN', 'IL', 'MN', 'IA', 'OH', 'WI', 'MI']),
        csb1724_points_on_road_shape = expand("data/edited/road/{state}_CSB1724_points_on_road.geojson", state=['IN', 'IL', 'MN', 'IA', 'OH', 'WI', 'MI']),
        csb1724_points_on_road_table = expand("data/edited/road/{state}_CSB1724_points_on_road_table.csv", state=['IN', 'IL', 'MN', 'IA', 'OH', 'WI', 'MI'])
    params:
        states = ['IN', 'IL', 'MN', 'IA', 'OH', 'WI', 'MI']
    script:
        "scripts/join_csb_supplement_4.py"


#---# Supplementary data 5: land management activities on neighboring fields (either intersecting or sharing a boundary)
rule join_regrow_supplement_5:
    input:
        regrow_joined_shape = expand("data/edited/Regrow/{state}_regrow_shape_table.geojson", state=['IN', 'IL', 'MN', 'IA', 'OH', 'WI', 'MI'])
    output:
        #regrow_supplement_5_shape = expand("data/edited/Regrow/{state}_regrow_supplement_5_spatial.geojson", state=['IN', 'IL', 'MN', 'IA', 'OH', 'WI', 'MI']),
        regrow_supplement_5_table = expand("data/edited/Regrow/{state}_regrow_supplement_5_table.csv", state=['IN', 'IL', 'MN', 'IA', 'OH', 'WI', 'MI'])
    params:
        states = ['IN', 'IL', 'MN', 'IA', 'OH', 'WI', 'MI']
    script:
        "scripts/join_regrow_supplement_5.py"

rule join_csb_supplement_5:
    input:
        csb1724_shape = expand("data/edited/CSB/{state}_CSB1724_clipped.gpkg", state=['IN', 'IL', 'MN', 'IA', 'OH', 'WI', 'MI'])
    output:
        #csb1724_supplement_5_shape = expand("data/edited/CSB/{state}_CSB1724_supplement_5_spatial.geojson", state=['IN', 'IL', 'MN', 'IA', 'OH', 'WI', 'MI']),
        csb1724_supplement_5_table = expand("data/edited/CSB/{state}_CSB1724_supplement_5_table.csv", state=['IN', 'IL', 'MN', 'IA', 'OH', 'WI', 'MI'])
    params:
        states = ['IN', 'IL', 'MN', 'IA', 'OH', 'WI', 'MI']
    script:
        "scripts/join_csb_supplement_5.py"

        
#---# Supplementary data 6: weather data 
rule join_regrow_supplement_6:
    input:
        regrow_shape = expand("data/edited/Regrow/{state}_regrow_fieldID_geometry.parquet", state=['IN', 'IL', 'MN', 'IA', 'OH', 'WI', 'MI']),
        #ppt_clipped = expand("data/edited/Weather/ppt/{state}_prism_ppt_us_30s_{year}{month}_clipped.tif", year=range(2014, 2026), month=[f"{m:02d}" for m in range(1, 13)], state=["OH", "MN", "WI", "IA", "IN", "IL", "MI"]),
        #tdmean_clipped =  expand("data/edited/Weather/tdmean/{state}_prism_tdmean_us_30s_{year}{month}_clipped.tif", year=range(2014, 2026), month=[f"{m:02d}" for m in range(1, 13)], state=["OH", "MN", "WI", "IA", "IN", "IL", "MI"]),
        #tmax_clipped =  expand("data/edited/Weather/tmax/{state}_prism_tmax_us_30s_{year}{month}_clipped.tif", year=range(2014, 2026), month=[f"{m:02d}" for m in range(1, 13)], state=["OH", "MN", "WI", "IA", "IN", "IL", "MI"]),
        #tmin_clipped =  expand("data/edited/Weather/tmin/{state}_prism_tmin_us_30s_{year}{month}_clipped.tif", year=range(2014, 2026), month=[f"{m:02d}" for m in range(1, 13)], state=["OH", "MN", "WI", "IA", "IN", "IL", "MI"]),
        #tmean_clipped =  expand("data/edited/Weather/tmean/{state}_prism_tmean_us_30s_{year}{month}_clipped.tif", year=range(2014, 2026), month=[f"{m:02d}" for m in range(1, 13)], state=["OH", "MN", "WI", "IA", "IN", "IL", "MI"]),
        #vpdmax_clipped =  expand("data/edited/Weather/vpdmax/{state}_prism_vpdmax_us_30s_{year}{month}_clipped.tif", year=range(2014, 2026), month=[f"{m:02d}" for m in range(1, 13)], state=["OH", "MN", "WI", "IA", "IN", "IL", "MI"]),
        #vpdmin_clipped = expand("data/edited/Weather/vpdmin/{state}_prism_vpdmin_us_30s_{year}{month}_clipped.tif", year=range(2014, 2026), month=[f"{m:02d}" for m in range(1, 13)], state=["OH", "MN", "WI", "IA", "IN", "IL", "MI"]),
        #soltotal_clipped = expand("data/edited/Weather/soltotal/{state}_prism_soltotal_us_30s_{year}{month}_clipped.tif", year=range(2014, 2026), month=[f"{m:02d}" for m in range(1, 13)], state=["OH", "MN", "WI", "IA", "IN", "IL", "MI"]),
        #solslope_clipped = expand("data/edited/Weather/solslope/{state}_prism_solslope_us_30s_{year}{month}_clipped.tif", year=range(2014, 2026), month=[f"{m:02d}" for m in range(1, 13)], state=["OH", "MN", "WI", "IA", "IN", "IL", "MI"])
    output:
        #regrow_supplement_6_shape = expand("data/edited/Regrow/{state}_regrow_supplement_6_spatial.geojson", state=['IN', 'IL', 'MN', 'IA', 'OH', 'WI', 'MI']),
        regrow_supplement_6_table = expand("data/edited/Regrow/{state}_regrow_supplement_6_table.parquet", state=['IN', 'IL', 'MN', 'IA', 'OH', 'WI', 'MI'])
    params:
        states = config["weather"]["states"],
        weather_variables = config["weather"]["weather_variables"],
        regrow_input_dir = config["regrow"]["edited_output_dir"],
        regrow_output_dir = config["regrow"]["edited_output_dir"]
    script:
        "scripts/join_regrow_supplement_6.py"

rule join_csb_supplement_6:
    input:
        csb1724_shape = expand("data/edited/CSB/{state}_CSB1724_CSBID_geometry.parquet", state=['IN', 'IL', 'MN', 'IA', 'OH', 'WI', 'MI']),
        ppt_clipped = expand("data/edited/Weather/ppt/{state}_prism_ppt_us_30s_{year}{month}_clipped.tif", year=range(2014, 2026), month=[f"{m:02d}" for m in range(1, 13)], state=["OH", "MN", "WI", "IA", "IN", "IL", "MI"]),
        tdmean_clipped =  expand("data/edited/Weather/tdmean/{state}_prism_tdmean_us_30s_{year}{month}_clipped.tif", year=range(2014, 2026), month=[f"{m:02d}" for m in range(1, 13)], state=["OH", "MN", "WI", "IA", "IN", "IL", "MI"]),
        tmax_clipped =  expand("data/edited/Weather/tmax/{state}_prism_tmax_us_30s_{year}{month}_clipped.tif", year=range(2014, 2026), month=[f"{m:02d}" for m in range(1, 13)], state=["OH", "MN", "WI", "IA", "IN", "IL", "MI"]),
        tmin_clipped =  expand("data/edited/Weather/tmin/{state}_prism_tmin_us_30s_{year}{month}_clipped.tif", year=range(2014, 2026), month=[f"{m:02d}" for m in range(1, 13)], state=["OH", "MN", "WI", "IA", "IN", "IL", "MI"]),
        tmean_clipped =  expand("data/edited/Weather/tmean/{state}_prism_tmean_us_30s_{year}{month}_clipped.tif", year=range(2014, 2026), month=[f"{m:02d}" for m in range(1, 13)], state=["OH", "MN", "WI", "IA", "IN", "IL", "MI"]),
        vpdmax_clipped =  expand("data/edited/Weather/vpdmax/{state}_prism_vpdmax_us_30s_{year}{month}_clipped.tif", year=range(2014, 2026), month=[f"{m:02d}" for m in range(1, 13)], state=["OH", "MN", "WI", "IA", "IN", "IL", "MI"]),
        vpdmin_clipped = expand("data/edited/Weather/vpdmin/{state}_prism_vpdmin_us_30s_{year}{month}_clipped.tif", year=range(2014, 2026), month=[f"{m:02d}" for m in range(1, 13)], state=["OH", "MN", "WI", "IA", "IN", "IL", "MI"]),
        soltotal_clipped = expand("data/edited/Weather/soltotal/{state}_prism_soltotal_us_30s_{year}{month}_clipped.tif", year=range(2014, 2026), month=[f"{m:02d}" for m in range(1, 13)], state=["OH", "MN", "WI", "IA", "IN", "IL", "MI"]),
        solslope_clipped = expand("data/edited/Weather/solslope/{state}_prism_solslope_us_30s_{year}{month}_clipped.tif", year=range(2014, 2026), month=[f"{m:02d}" for m in range(1, 13)], state=["OH", "MN", "WI", "IA", "IN", "IL", "MI"])
    output:
        #csb1724_supplement_6_shape = expand("data/edited/CSB/{state}_CSB1724_supplement_6_spatial.geojson", state=['IN', 'IL', 'MN', 'IA', 'OH', 'WI', 'MI']),
        csb1724_supplement_6_table = expand("data/edited/CSB/{state}_CSB1724_supplement_6_table.parquet", state=['IN', 'IL', 'MN', 'IA', 'OH', 'WI', 'MI'])
    params:
        states = config["weather"]["states"],
        weather_variables = config["weather"]["weather_variables"],
        CSB_input_dir = config["csb"]["edited_output_dir"],
        CSB_output_dir = config["csb"]["edited_output_dir"]
    script:
        "scripts/join_csb_supplement_6.py"

#---# Supplementary data 7: crop price data
rule join_regrow_supplement_7:
    input:
        regrow_joined_shape = expand("data/edited/Regrow/{state}_regrow_supplement_1_spatial.geojson", state=['IN', 'IL', 'MN', 'IA', 'OH', 'WI', 'MI']),
        corn_price = expand("data/Grain Price/{state}_corn_{level}_level.xlsx", state=["OH", "IN", "MI", "IA", "IL", "WI", "MN"], level=["elevator", "county"]),
        soybean_price = expand("data/Grain Price/{state}_soybeans_{level}_level.xlsx", state=["OH", "IN", "MI", "IA", "IL", "WI", "MN"], level=["elevator", "county"]),
        wheat_price = expand("data/Grain Price/{state}_wheat_{level}_level.xlsx", state=["OH", "IN", "MI", "IA", "IL", "WI", "MN"], level=["elevator"]),
        elevator_location_geojson = expand("data/edited/Grain Price/{crop}_elevator_location.geojson", crop=['corn', 'soybeans', 'wheat']),
        index_county_location_geojson = expand("data/edited/Grain Price/{crop}_index_county_location.geojson", crop=['corn', 'soybeans'])
    output:
        #regrow_supplement_7_shape = expand("data/edited/Regrow/{state}_regrow_supplement_7_spatial.geojson", state=['IN', 'IL', 'MN', 'IA', 'OH', 'WI', 'MI']),
        regrow_supplement_7_table = expand("data/edited/Regrow/{state}_regrow_supplement_7_table.csv", state=['IN', 'IL', 'MN', 'IA', 'OH', 'WI', 'MI'])
    params:
        states = config["price"]["states"],
        crops = config["price"]["crops"],
        number_of_neighbors = config["price"]["N_nearest"]
    script:
        "scripts/join_regrow_supplement_7.py"

rule join_csb_supplement_7:
    input:
        csb1724_shape = expand("data/edited/CSB/{state}_CSB1724_supplement_1_spatial.geojson", state=['IN', 'IL', 'MN', 'IA', 'OH', 'WI', 'MI']),
        corn_price = expand("data/Grain Price/{state}_corn_{level}_level.xlsx", state=["OH", "IN", "MI", "IA", "IL", "WI", "MN"], level=["elevator", "county"]),
        soybean_price = expand("data/Grain Price/{state}_soybeans_{level}_level.xlsx", state=["OH", "IN", "MI", "IA", "IL", "WI", "MN"], level=["elevator", "county"]),
        wheat_price = expand("data/Grain Price/{state}_wheat_{level}_level.xlsx", state=["OH", "IN", "MI", "IA", "IL", "WI", "MN"], level=["elevator"]),
        elevator_location_geojson = expand("data/edited/Grain Price/{crop}_elevator_location.geojson", crop=['corn', 'soybeans', 'wheat']),
        index_county_location_geojson = expand("data/edited/Grain Price/{crop}_index_county_location.geojson", crop=['corn', 'soybeans'])
    output:
        #csb1724_supplement_7_shape = expand("data/edited/CSB/{state}_CSB1724_supplement_7_spatial.geojson", state=['IN', 'IL', 'MN', 'IA', 'OH', 'WI', 'MI']),
        csb1724_supplement_7_table = expand("data/edited/CSB/{state}_CSB1724_supplement_7_table.csv", state=['IN', 'IL', 'MN', 'IA', 'OH', 'WI', 'MI'])
    params:
        states = config["price"]["states"],
        crops = config["price"]["crops"],
        number_of_neighbors = config["price"]["N_nearest"]
    script:
        "scripts/join_csb_supplement_7.py"

rule convert_csv_to_parquet:
    script: "scripts/convert_csv_to_parquet.py"

#---# Supplementary data 8: soil composition data
rule join_regrow_supplement_8:
    input:
        regrow_geometry = expand("data/edited/Regrow/{state}_regrow_fieldID_geometry.parquet", state=['IN', 'IL', 'MN', 'IA', 'OH', 'WI', 'MI']),
        regrow_raster_to_gSSURGO_grid = expand("data/edited/Regrow/{state}_regrow_raster_to_gSSURGO_grid.tif", state=["OH", "MN", "WI", "IA", "IN", "IL", "MI"]),
        regrow_duplicating_fields = expand("data/edited/Regrow/{state}_regrow_duplicating_fields.parquet", state=["OH", "MN", "WI", "IA", "IN", "IL", "MI"]),
        gSSURGO_tabular = expand("data/gSSURGO/gSSURGO_{state}/gSSURGO_{state}.gdb", state=["OH", "IN", "MI", "IA", "IL", "WI", "MN"])
    output:
        #regrow_supplement_8_shape = expand("data/edited/Regrow/{state}_regrow_supplement_8_spatial.parquet", state=['IN', 'IL', 'MN', 'IA', 'OH', 'WI', 'MI']),
        regrow_supplement_8_table = expand("data/edited/Regrow/{state}_regrow_supplement_8_table.parquet", state=['IN', 'IL', 'MN', 'IA', 'OH', 'WI', 'MI'])
    params:
        states = config["soil"]["states"],
        soil_depth_cm = config["soil"]["soil_depth_cm"],
        regrow_input_dir = config["regrow"]["edited_output_dir"],
        regrow_output_dir = config["regrow"]["edited_output_dir"],
        soil_output_dir = config["soil"]["edited_output_dir"]
    script:
        "scripts/join_regrow_supplement_8.py"

rule join_csb_supplement_8:
    input:
        csb1724_geometry = expand("data/edited/CSB/{state}_CSB1724_CSBID_geometry.parquet", state=['IN', 'IL', 'MN', 'IA', 'OH', 'WI', 'MI']),
        CSB_raster_to_gSSURGO_grid = expand("data/edited/CSB/{state}_CSB1724_raster_to_gSSURGO_grid.tif", state=["OH", "MN", "WI", "IA", "IN", "IL", "MI"]),
        CSB_duplicating_fields = expand("data/edited/CSB/{state}_CSB1724_duplicating_fields.parquet", state=["OH", "MN", "WI", "IA", "IN", "IL", "MI"]),
        gSSURGO_tabular = expand("data/gSSURGO/gSSURGO_{state}/gSSURGO_{state}.gdb", state=["OH", "IN", "MI", "IA", "IL", "WI", "MN"])
    output:
        #csb1724_supplement_8_shape = expand("data/edited/CSB/{state}_CSB1724_supplement_8_spatial.parquet", state=['IN', 'IL', 'MN', 'IA', 'OH', 'WI', 'MI']),
        csb1724_supplement_8_table = expand("data/edited/CSB/{state}_CSB1724_supplement_8_table.parquet", state=['IN', 'IL', 'MN', 'IA', 'OH', 'WI', 'MI'])
    params:
        states = config["soil"]["states"],
        soil_depth_cm = config["soil"]["soil_depth_cm"],
        CSB_input_dir = config["csb"]["edited_output_dir"],
        CSB_output_dir = config["csb"]["edited_output_dir"],
        soil_output_dir = config["soil"]["edited_output_dir"]
    script:
        "scripts/join_csb_supplement_8.py"
