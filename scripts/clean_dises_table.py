import pandas as pd
import os

df = pd.read_csv("data/DISES/combined_data_clean.csv")

# Delete comprehensive ID = NA
df_wip = df[df['Comprehensive_ID'].notna()].copy()

# Keep only columns comprehensive ID, field name, field size, field crop
df_final = df_wip[['Comprehensive_ID', 'field_name', 'field_size', 'field_crop']].copy()

# Rename column Comprehensive_ID to comp_id
df_final.rename(columns={'Comprehensive_ID': 'comp_id'}, inplace=True)

# Save dises data table for future join with dises shape
os.makedirs("data/edited/DISES", exist_ok=True)
df_final.to_csv("data/edited/DISES/combined_data_clean_short.csv", index=False)