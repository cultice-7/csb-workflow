import geopandas as gp
import pandas as pd

gdf = gp.read_file("data/DISES/DISES_All_Parcels_05.15.25.shp")

# Delete comprehensive ID = 0
gdf_wip = gdf[gdf['Comprehe_1'] != 0].copy()

# Rename comprehensive ID column from Comprehe_1 to Comp_ID
gdf_wip.rename(columns={'Comprehe_1': 'comp_id'}, inplace=True)

# Keep only comp_id and geometry
gdf_minimal = gdf_wip[['comp_id', 'geometry']]

# Consolidate multiple tax parcels into one for each owner
gdf_consol = gdf_minimal.dissolve(by='comp_id')

# Reset index so comp_id becomes a column again
gdf_consol.reset_index(inplace=True)

# Load dises data table
csv = pd.read_csv("data/edited/DISES/combined_data_clean_short.csv")

# Add indication that this comp_id exist in data table, to be shown in joined shape file
csv['from_data_table'] = "Y"
cols = ['comp_id', 'from_data_table', 'field_name', 'field_size', 'field_crop']
csv = csv[cols]

# Replace field_crop values to match CDL classification (We'll do the same with Regrow)
csv['field_crop'] = csv['field_crop'].replace(2,5)

# Join shape file and table using comprehensive ID
gdf_joined = gdf_consol.merge(csv, on='comp_id', how='left')

# Save shape file to GeoPackage (preserve blanks in numeric fields) & save attribute table as csv
gdf_joined.to_file("data/edited/DISES/dises_consolidated.gpkg", layer='dises_consolidated', driver="GPKG")
gdf_joined.drop(columns='geometry').to_csv("data/edited/DISES/dises_consolidated_attributes.csv", index=False)