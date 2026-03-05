import geopandas as gp
import os
import gc

### Note: Two CSB files are processed separately due to their large size ###

# Load CSB 2016-2023
csb1623 = gp.read_file("data/CSB/CSB1623.gdb")

# Define FIPS codes for states
# IL - 17, IN - 18, IA - 19, MI - 26, MN - 27, OH -  39, WI - 55
target_fips = ['17', '18', '19', '26', '27', '39', '55']
state_name =  ['IL', 'IN', 'IA', 'MI', 'MN', 'OH', 'WI']

for state_num in target_fips:
    # Clip CSBs to the selected states
    csb1623_clipped = csb1623[csb1623['STATEFIPS'] == state_num]

    # Save to new files
    os.makedirs("data/edited/CSB", exist_ok=True)
    pos = target_fips.index(state_num)
    csb1623_clipped.to_file(f"data/edited/CSB/{state_name[pos]}_CSB1623_clipped.gpkg", driver="GPKG")
    print(f"CSB1623 for {state_name[pos]} is clipped")

# Delete the files to free memory and force garbage collection
del csb1623
del csb1623_clipped
gc.collect()



# Load CSB 2017-2024
csb1724 = gp.read_file("data/CSB/CSB1724.gdb")

# Define FIPS codes for states
# IL - 17, IN - 18, IA - 19, MI - 26, MN - 27, OH -  39, WI - 55
target_fips = ['17', '18', '19', '26', '27', '39', '55']
state_name =  ['IL', 'IN', 'IA', 'MI', 'MN', 'OH', 'WI']

for state_num in target_fips:
    # Clip CSBs to the selected states
    csb1724_clipped = csb1724[csb1724['STATEFIPS'] == state_num]

    # Save to the existimg directory
    pos = target_fips.index(state_num)
    csb1724_clipped.to_file(f"data/edited/CSB/{state_name[pos]}_CSB1724_clipped.gpkg", driver="GPKG")
    print(f"CSB1724 for {state_name[pos]} is clipped")

# Delete the files to free memory and force garbage collection
del csb1724
del csb1724_clipped
gc.collect()