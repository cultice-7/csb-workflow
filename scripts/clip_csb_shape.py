import geopandas as gp
from pathlib import Path
import gc

### Note: Two CSB files are processed separately due to their large size ###

# FIPS code to state abbreviation mapping
FIPS_TO_STATE = {
    "17": "IL",
    "18": "IN",
    "19": "IA",
    "26": "MI",
    "27": "MN",
    "39": "OH",
    "55": "WI",
}

# Get states from snakemake params
states = snakemake.params.states

# Build FIPS list from requested states (reverse lookup)
STATE_TO_FIPS = {v: k for k, v in FIPS_TO_STATE.items()}
target_fips = [STATE_TO_FIPS[s] for s in states]

output_dir = Path("data/edited/CSB")
output_dir.mkdir(parents=True, exist_ok=True)

# Load and clip CSB 2016-2023
csb1623 = gp.read_file("data/CSB/CSB1623.gdb")

for fips in target_fips:
    state_name = FIPS_TO_STATE[fips]
    csb1623_clipped = csb1623[csb1623["STATEFIPS"] == fips]
    csb1623_clipped.to_file(
        output_dir / f"{state_name}_CSB1623_clipped.gpkg", driver="GPKG"
    )
    print(f"CSB1623 for {state_name} is clipped")

# Delete to free memory and force garbage collection
del csb1623
gc.collect()


# Load and clip CSB 2017-2024
csb1724 = gp.read_file("data/CSB/CSB1724.gdb")

for fips in target_fips:
    state_name = FIPS_TO_STATE[fips]
    csb1724_clipped = csb1724[csb1724["STATEFIPS"] == fips]
    csb1724_clipped.to_file(
        output_dir / f"{state_name}_CSB1724_clipped.gpkg", driver="GPKG"
    )
    print(f"CSB1724 for {state_name} is clipped")

# Delete to free memory and force garbage collection
del csb1724
gc.collect()
