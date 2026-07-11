import pandas as pd
import os

df = pd.read_csv("data/Regrow/OH_tillage_june24.csv")

# Change the format of date columns to datetime
cols_to_convert = ['start_date', 'end_date', 'observation_window_start_date', 'observation_window_end_date']
df[cols_to_convert] = df[cols_to_convert].apply(pd.to_datetime, errors = 'coerce')

# Change the format of practice_determination to str
# df['practice_determination'] = df['practice_determination'].astype(str) 

# Replace missing values in "start_date" with "observation_window_start_date"
# Replace gaps in "end_date" with values from "end_window_start_date"
df['start_date'] = df['start_date'].fillna(df['observation_window_start_date'])
df['end_date'] = df['end_date'].fillna(df['observation_window_end_date'])

# Extract year from date columns
df['start_year'] = df['start_date'].dt.year
df['end_year'] = df['end_date'].dt.year

# Compute tillage length (days)
df['tillage_days'] = (df['end_date'] - df['start_date']).dt.days

# Sort dataset by boundary_id, start_year, plant_date (so that earlier tillage within each year comes first)
df.sort_values(by = ['boundary_id', 'start_year', 'start_date'], inplace = True)
df.reset_index(drop = True, inplace = True)

# Compute the number of tillages for each field within each year
df['tillage_count'] = df.groupby(['boundary_id', 'start_year']).cumcount() + 1

# Create year_tillage_count identifier
df['year_tillage_count'] = df['start_year'].astype(str).str[2:] + '_' + df['tillage_count'].astype(str)

# Reshape from long to wide
df_wide = df.pivot(index = 'boundary_id', columns = 'year_tillage_count', values = 'practice_determination')
df_wide.reset_index(inplace = True)
df_wide.columns = ['field_id'] + ['till' + str(col) for col in df_wide.columns[1:]]

# Recode tillage names to numerical codes 
tillage_map = {
    "TILLAGE_INTENSITY_NO_TILLAGE": "no_till", 
    "TILLAGE_INTENSITY_REDUCED_TILLAGE": "reduced", 
    "TILLAGE_INTENSITY_CONVENTIONAL_TILLAGE": "conventional",
    "TILLAGE_INTENSITY_NO_TILLAGE_DATA": "NA"
}

df_wide.replace(tillage_map, inplace = True)

# Save to csv
os.makedirs("data/edited/Regrow", exist_ok = True)
df_wide.to_csv("data/edited/Regrow/OH_tillage_wide_coded.csv", index = False)