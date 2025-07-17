#untested

import geopandas as gp

# Load CSBs
csb1623 = gp.read_file("data/CSB/CSB1623.gdb")
csb1724 = gp.read_file("data/CSB/CSB1724.gdb")

# Define FIPS codes for OH, MI, IN
target_fips = ['39', '26', '18']

# Clip CSBs to the selected states
csb1623_clipped = csb1623[csb1623['STATEFIPS'].isin(target_fips)]
csb1724_clipped = csb1724[csb1724['STATEFIPS'].isin(target_fips)]

# Save to new GDB files
csb1623_clipped.to_file("data/edited/CSB/CSB1623_clipped.gdb", driver="FileGDB")
csb1724_clipped.to_file("data/edited/Regrow/CSB1724_clipped.gdb", driver="FileGDB")
