import pandas as pd
import os

df = pd.read_csv("data/Regrow/OH_green_cover_crop_july7.csv")

# Change the format of date columns to datetime
cols_to_convert = ['start_date', 'end_date']
df[cols_to_convert] = df[cols_to_convert].apply(pd.to_datetime, errors = 'coerce')

# Change the format of practice_determination to string
# df['practice_determination'] = df['practice_determination'].astype(str) 

# Extract year from date columns
df['start_year'] = df['start_date'].dt.year
df['end_year'] = df['end_date'].dt.year

# Compute planting length (days)
df['planting_days'] = (df['end_date'] - df['start_date']).dt.days

# Sort dataset by boundary_id, start_year, start_date (so that earlier plantings within each year come first)
df.sort_values(by = ['boundary_id', 'start_year', 'start_date'], inplace = True)
df.reset_index(drop = True, inplace = True)

# Multicropping detection
df['multicrop'] = df.groupby(['boundary_id', 'start_year']).cumcount() + 1
#df['maxmulticrop'] = df.groupby(['boundary_id', 'plant_year'])['multicrop'].transform('max')
#df.drop(columns=['maxmulticrop'], inplace=True)

# Create year_cropset identifier
df['year_cropset'] = df['start_year'].astype(str).str[2:] + '_' + df['multicrop'].astype(str)

# Drop unnecessary columns
#df.drop(columns = ['start_date', 'end_date', 'overall_confidence', 'start_year', 'end_year', 'multicrop'], inplace=True)

# Reshape from long to wide
df_wide = df.pivot(index = 'boundary_id', columns = 'year_cropset', values = 'practice_determination')
df_wide.reset_index(inplace = True)
df_wide.columns = ['field_id'] + ['cover' + str(col) for col in df_wide.columns[1:]]

# Recode crop names to CDL codes 
# Note: recode non_cropland as 999 which doesn't exist in CDL, 
# recode grass_perennial as 176 which is Grassland/Pasture in CDL,
# recode hay and ryegrass as 37 which is Other Hay/Non Alfalfa
crop_map = {
    "GREEN_COVER_CLASS_COVER_CROP_NOT_TRACKED": "no_crop", 
    "GREEN_COVER_CLASS_POTENTIAL_COVER_CROP": "potential", 
    "GREEN_COVER_CLASS_COVER_CROP": "cover_crop",
    "GREEN_COVER_CLASS_NO_DATA": "NA",
    "GREEN_COVER_CLASS_NOT_APPLICABLE": "NA"
}

df_wide.replace(crop_map, inplace = True)

# Save to csv
os.makedirs("data/edited/Regrow", exist_ok = True)
df_wide.to_csv("data/edited/Regrow/OH_cover_crop_wide_coded.csv", index = False)