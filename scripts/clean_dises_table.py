import pandas as pd

df = pd.read_csv("../data/DISES/combined_data_clean.csv")

# Delete comprehensive ID = NA
df_wip = df[df['Comprehensive_ID'] != 'NA']

# Keep only columns comprehensive ID, field name, field size, field crop
df_final = df_wip[['Comprehensive_ID', 'field_name', 'field_size' , 'field_crop']]

# Save dises data table for future join with dises shape
df_final.to_csv("../data/edited/DISES/combined_data_clean_short.csv")