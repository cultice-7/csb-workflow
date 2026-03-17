# =============================================================================
# DOWNLOAD RULES
# Download external datasets (excluding DISES and Regrow which require manual download)
# =============================================================================


rule download_census_tract:
    output:
        tract = "data/Census/census_tract/cb_2023_us_tract_500k.shp"
    params:
        raw_dir = config["census_tract"]["raw_dir"],
        output_dir = config["census_tract"]["output_dir"],
        html = f'{config["census_tract"]["base_html"]}/cb_2023_us_tract_500k.zip'
    log: "logs/download_census_tract.log"
    script:
        "../scripts/download_census_tract.py"


rule download_state_bound:
    output:
        state_bound = "data/Census/state_bound/cb_2023_us_state_500k.shp"
    params:
        raw_dir = config["state_bound"]["raw_dir"],
        output_dir = config["state_bound"]["output_dir"],
        html = f'{config["state_bound"]["base_html"]}/cb_2023_us_state_500k.zip'
    log: "logs/download_state_bound.log"
    script:
        "../scripts/download_state_bound.py"


rule download_county_bound:
    output:
        county_bound = "data/Census/county_bound/cb_2023_us_county_500k.shp"
    params:
        raw_dir = config["county_bound"]["raw_dir"],
        output_dir = config["county_bound"]["output_dir"],
        html = f'{config["county_bound"]["base_html"]}/cb_2023_us_county_500k.zip'
    log: "logs/download_county_bound.log"
    script:
        "../scripts/download_county_bound.py"


rule download_roads:
    output:
        roads = "data/Census/road/prisecroads.shp"
    params:
        prefixes = config["roads"]["prefixes"],
        raw_dir = config["roads"]["raw_dir"],
        output_dir = config["roads"]["output_dir"],
        html = config["roads"]["base_html"]
    log: "logs/download_roads.log"
    script:
        "../scripts/download_roads.py"


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
    log: "logs/download_3dep.log"
    threads: 2
    resources:
        mem_mb = 8000
    script:
        "../scripts/download_3dep.py"


rule download_watershed:
    output:
        subbasin = "data/Geo/watershed/subbasin.shp",
        watershed = "data/Geo/watershed/watershed.shp",
        subwatershed = "data/Geo/watershed/subwatershed.shp"
    params:
        states = config["watershed"]["states"],
        raw_dir = config["watershed"]["raw_dir"],
        output_dir = config["watershed"]["output_dir"],
        html = config["watershed"]["base_html"]
    log: "logs/download_watershed.log"
    script:
        "../scripts/download_watershed.py"


rule download_weather:
    output:
        ppt     = expand("data/Weather/ppt/prism_ppt_us_30s_{year}{month}.tif",         year=YEARS, month=MONTHS),
        tmax    = expand("data/Weather/tmax/prism_tmax_us_30s_{year}{month}.tif",        year=YEARS, month=MONTHS),
        tmin    = expand("data/Weather/tmin/prism_tmin_us_30s_{year}{month}.tif",        year=YEARS, month=MONTHS),
        tmean   = expand("data/Weather/tmean/prism_tmean_us_30s_{year}{month}.tif",      year=YEARS, month=MONTHS),
        tdmean  = expand("data/Weather/tdmean/prism_tdmean_us_30s_{year}{month}.tif",    year=YEARS, month=MONTHS),
        vpdmax  = expand("data/Weather/vpdmax/prism_vpdmax_us_30s_{year}{month}.tif",    year=YEARS, month=MONTHS),
        vpdmin  = expand("data/Weather/vpdmin/prism_vpdmin_us_30s_{year}{month}.tif",    year=YEARS, month=MONTHS),
        soltotal = expand("data/Weather/soltotal/prism_soltotal_us_30s_{year}{month}.tif", year=YEARS, month=MONTHS),
        solslope = expand("data/Weather/solslope/prism_solslope_us_30s_{year}{month}.tif", year=YEARS, month=MONTHS),
    params:
        weather_variables = WEATHER_VARS,
        states = config["weather"]["states"],
        year_range = config["weather"]["year"],
        raw_dir = config["weather"]["raw_dir"],
        output_dir = config["weather"]["output_dir"],
        html = config["weather"]["base_html"]
    log: "logs/download_weather.log"
    threads: 4
    resources:
        mem_mb = 4000
    script:
        "../scripts/download_weather.py"


rule download_csb_1623:
    output:
        raw_csb = directory("data/CSB/CSB1623.gdb"),
        metadata = "data/CSB/NationalCSB_2016-2023_rev23_metadata.htm"
    params:
        raw_dir = config["csb"]["raw_dir"],
        output_dir = config["csb"]["output_dir"],
        html = f'{config["csb"]["base_html"]}/NationalCSB_2016-2023_rev23.zip'
    log: "logs/download_csb_1623.log"
    script:
        "../scripts/download_csb.py"


rule download_csb_1724:
    output:
        raw_csb = directory("data/CSB/CSB1724.gdb"),
        metadata = "data/CSB/NationalCSB_2017-2024_rev23_metadata.htm"
    params:
        raw_dir = config["csb"]["raw_dir"],
        output_dir = config["csb"]["output_dir"],
        html = f'{config["csb"]["base_html"]}/NationalCSB_2017-2024_rev23.zip'
    log: "logs/download_csb_1724.log"
    script:
        "../scripts/download_csb.py"


rule download_cdl:
    output:
        raw_cdl = expand("data/CDL/{year}_30m_cdls/{year}_30m_cdls.tif", year=YEARS)
    params:
        year_range = config["cdl"]["year"],
        raw_dir = config["cdl"]["raw_dir"],
        output_dir = config["cdl"]["output_dir"],
        html = config["cdl"]["base_html"]
    log: "logs/download_cdl.log"
    script:
        "../scripts/download_cdl.py"
