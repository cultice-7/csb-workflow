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