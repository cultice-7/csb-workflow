import geopandas as gpd
import os

# Import parameters from Snakemake
CSB_input_folder = snakemake.params.CSB_input_dir
CSB_output_folder = snakemake.params.CSB_output_dir
states = snakemake.params.states
states_codes = snakemake.params.states_codes
CSB_years = snakemake.params.CSB_years


# Clip each 8‑year CSB file by state
def clip_raw_CSB(states, state_code, CSB_year, CSB_input_folder, CSB_output_folder):
    
    # Paths to input CSB files
    csb_input_path =  os.path.join(CSB_input_folder, f"CSB{CSB_year}.gdb")
    
    # Create an ouput directory
    os.makedirs(CSB_output_folder, exist_ok=True)

    # Clip CSB to the selected state
    csb_clipped = gpd.read_file(
        csb_input_path,
        where=f"STATEFIPS = '{state_code}'"
    )

    # Save clipped CSB file
    pos = states_codes.index(state_code)
    csb_output_path =  os.path.join(CSB_output_folder, f"{states[pos]}_CSB{CSB_year}_clipped.parquet")
    csb_clipped.to_parquet(csb_output_path, compression="zstd")
    
    print(f"CSB{CSB_year} for {states[pos]} ({pos}) is clipped")


# Main code
for CSB_year in CSB_years:
    for state_code in states_codes:
        clip_raw_CSB(states, state_code, CSB_year, CSB_input_folder, CSB_output_folder)
