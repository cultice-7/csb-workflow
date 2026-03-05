import geopandas as gp
import pandas as pd

gdf = gp.read_file("data/DISES/DISES_All_Parcels_11.12.25.shp")

# Delete comprehensive ID = 0
gdf_wip = gdf[gdf['comprehens'] != 0].copy()

# Rename comprehensive ID column from Comprehe_1 to Comp_ID
gdf_wip.rename(columns={'comprehens': 'comp_id'}, inplace=True)

# Keep only comp_id and geometry
gdf_minimal = gdf_wip[['comp_id', 'geometry']]

# Consolidate multiple tax parcels into one for each owner
gdf_consol = gdf_minimal.dissolve(by='comp_id')

# Reset index so comp_id becomes a column again
gdf_consol.reset_index(inplace=True)

# Reproject to NAD83/Conus Albers
gdf_consol = gdf_consol.to_crs(epsg=5070)

#  Calculate area in US survey acres (1 acre = 4046.8564224 sq m)
gdf_consol["area_acre"] = gdf_consol.geometry.area/4046.8564224

# Load dises data table
df_csv = pd.read_csv("data/edited/DISES/combined_data_clean_all_columns.csv")

# Add indication that this comp_id exist in data table, to be shown in joined shape file
df_csv['survey_responded'] = "Y"

# Replace field_crop values to match CDL classification (We'll do the same with Regrow)
df_csv['field_crop_23'] = df_csv['field_crop_23'].replace(2,5)

# Replace field_CC values to match Regrow
df_csv['field_cover_23'] = df_csv['field_cover_23'].replace(1,3) #1-Yes to 3-confirmed CC
df_csv['field_cover_23'] = df_csv['field_cover_23'].replace(2,1) #2-No to 1-No cC

# Join shape file and table using comprehensive ID
gdf_joined = gdf_consol.merge(df_csv, on='comp_id', how='left')

# Save shape file to GeoPackage (preserve blanks in numeric fields) & save attribute table as csv
gdf_joined.to_file("data/edited/DISES/dises_consolidated.gpkg", layer='dises_consolidated', driver="GPKG")
gdf_joined.drop(columns='geometry').to_csv("data/edited/DISES/dises_consolidated_attributes.csv", index=False)