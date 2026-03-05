import pandas as pd
import os

df = pd.read_csv("data/DISES/combined_data_clean.csv")

# Delete comprehensive ID = NA
df_wip = df[df['Comprehensive_ID'].notna()].copy()

# Keep only specific columns 
#column_names = ['Comprehensive_ID', 'tot_acres', 'tot_acres_owned', 'tot_acres_rented',
#               'succession', 'livestock', 'field_size', 'profit_good', 'profit_bad',
#               'field_crop', 'field_tillage', 'field_CC',
#               'FI_1', 'FI_2', 'FI_3', 'FI_4', 'FI_5', 'FI_6', 'FI_7', 'FI_8', 'FI_9',
#                'FI_10', 'FI_11', 'FI_12', 'FI_13', 'FI_14', 'FI_15', 'FI_16',
#                'Sex', 'Age', 'Education', 'NR_education', 'OFI', 'involved_org', 'involved_farm'] 
#df_final = df_wip[column_names].copy()

# Keep all columns 
df_final = df_wip.copy()

# Rename column Comprehensive_ID to comp_id, field_crop to field_crop_23, field_tillage to field_till_23, field_CC to field_cover_23  
df_final.rename(columns={'Comprehensive_ID': 'comp_id',
                         'field_crop': 'field_crop_23',
                         'field_tillage': 'field_till_23',
                         'field_CC': 'field_cover_23'}, inplace = True)

# Compute farmer types using FI answers from the survey
# Productivism index
df_final['productivism_index'] = df_final[['FI_1', 'FI_4', 'FI_7', 'FI_9', 'FI_10']].mean(axis=1)
# Conservationism index
df_final['conservationism_index'] = df_final[['FI_2', 'FI_3', 'FI_5', 'FI_6', 'FI_8']].mean(axis=1)
# Civic index
df_final['civic_index'] = df_final[['FI_11', 'FI_12', 'FI_13', 'FI_14', 'FI_15', 'FI_16']].mean(axis=1)

# Save dises data table for future join with dises shape
os.makedirs("data/edited/DISES", exist_ok = True)
#df_final.to_csv("data/edited/DISES/combined_data_clean_short.csv", index = False)
df_final.to_csv("data/edited/DISES/combined_data_clean_all_columns.csv", index = False)