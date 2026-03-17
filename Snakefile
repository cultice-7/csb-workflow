# =============================================================================
# CSB Workflow - Crop Sequence Boundary Pipeline
# =============================================================================
#
# This Snakemake pipeline automates the download, cleaning, clipping, and
# spatial joining of agricultural geospatial data for 7 US Midwestern states.
#
# Two main datasets are processed in parallel:
#   - Regrow: Field-level agricultural monitoring data
#   - CSB: USDA Crop Sequence Boundaries
#
# Each dataset is enriched with 8 supplementary data layers:
#   1. Census tract boundaries
#   2. Slope & elevation (3DEP)
#   3. Watershed boundaries (NHD)
#   4. Nearest road distance (TIGER)
#   5. Neighboring field activities
#   6. Weather data (PRISM)
#   7. Crop prices (grain elevators & county indices)
#   8. Soil composition (gSSURGO)
#
# Prerequisites:
#   - DISES and Regrow datasets must be downloaded manually (not automated)
#   - gSSURGO geodatabases must be downloaded manually per state
#   - Grain price Excel files must be obtained manually
#   - conda env: conda env create -f envs/env.yml
#
# Usage:
#   snakemake --cores <N>             # Run full pipeline
#   snakemake <rule_name> --cores <N> # Run a specific rule
#   snakemake --dag | dot -Tpng > dag.png  # Generate DAG image
#
# =============================================================================

configfile: "config.yml"


# =============================================================================
# Constants
# =============================================================================

STATES = config["elevation"]["states"]          # ['IN', 'IL', 'MN', 'IA', 'OH', 'WI', 'MI']
DISES_STATES = ["OH", "IN", "MI"]               # States with DISES data
YEARS = range(2014, 2025)                        # CDL/weather year range (2014-2024)
CSB_YEARS = ["1623", "1724"]                     # CSB dataset year ranges
WEATHER_VARS = config["weather"]["weather_variables"]  # 9 PRISM variables
MONTHS = [f"{m:02d}" for m in range(1, 13)]      # 01-12


# =============================================================================
# Target rule: all final outputs
# =============================================================================

rule all:
    input:
        # Regrow CDL validation
        expand("data/edited/Regrow/validation/{state}_regrow_with_cdl_validation.geojson",
               state=STATES),
        expand("data/edited/Regrow/validation/{state}_regrow_validity_summary_by_year.csv",
               state=STATES),
        # Regrow-DISES spatial join
        expand("data/edited/Regrow/{state}_regrow_dises_table.parquet",
               state=DISES_STATES),
        # CSB-DISES spatial join
        expand("data/edited/CSB/{state}_CSB{years}_dises_table.parquet",
               state=DISES_STATES, years=CSB_YEARS),
        # Regrow supplement tables 1-7 (CSV)
        expand("data/edited/Regrow/{state}_regrow_supplement_{n}_table.csv",
               state=STATES, n=range(1, 8)),
        # Regrow supplement table 8 (Parquet - soil)
        expand("data/edited/Regrow/{state}_regrow_supplement_8_table.parquet",
               state=STATES),
        # CSB supplement tables 1-7 (CSV)
        expand("data/edited/CSB/{state}_CSB1724_supplement_{n}_table.csv",
               state=STATES, n=range(1, 8)),
        # CSB supplement table 8 (Parquet - soil)
        expand("data/edited/CSB/{state}_CSB1724_supplement_8_table.parquet",
               state=STATES),


# =============================================================================
# Include rule modules
# =============================================================================

include: "rules/download.smk"
include: "rules/clean.smk"
include: "rules/clip.smk"
include: "rules/joins.smk"
include: "rules/supplements.smk"
