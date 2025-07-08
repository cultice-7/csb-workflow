import geopandas as gp
import pandas as pd

gdf = gp.read_file("data/DISES_All_Parcels_05.15.25.shp")

# Delete comprehensive ID = 0
gdf_wip = gdf[gdf['Comprehe_1'] != 0]

# Rename comprehensive ID column from Comprehe_1 to Comprehensive_ID
gdf_wip.rename(columns={'Comprehe_1': 'Comprehensive_ID'}, inplace=True)

# Consolidate multiple tax parcels into one for each owner
gdf_consol = gdf_wip.dissolve(by='Comprehensive_ID')

# Load dises data table
csv = pd.read_csv("combined_data_clean_short.csv")

# Join shape file and table using comprehensive ID
gdf_joined = gdf_consol.merge(csv, on='Comprehensive_ID', how='left')

# Save shape file
gdf_joined.to_file("dises_consolidated.shp")