import pandas as pd

df = pd.read_csv("data/Regrow/OH_main_crop_june24.csv")

# Extract year from date strings
df['plantyear'] = df['plant_date'].str[:4]
df['harvestyear'] = df['harvest_date'].str[:4]

# Multicropping detection
df['multicrop'] = df.groupby(['boundary_id', 'plantyear']).cumcount() + 1
df['maxmulticrop'] = df.groupby(['boundary_id', 'plantyear'])['multicrop'].transform('max')
df.drop(columns=['maxmulticrop'], inplace=True)

# Create year_cropset identifier
df['year_cropset'] = df['plantyear'].astype(str) + '_' + df['multicrop'].astype(str)

# Drop unnecessary columns
df.drop(columns=['plant_date', 'harvest_date', 'crop_type_confidence', 'plantyear', 'multicrop'], inplace=True)

# Reshape from long to wide
df_wide = df.pivot(index='boundary_id', columns='year_cropset', values='practice_determination')
df_wide.reset_index(inplace=True)
df_wide.columns = ['boundary_id'] + ['crop' + str(col) for col in df_wide.columns[1:]]

# Recode crop names to CDL codes 
# Note: recode non_cropland as 999 which doesn't exist in CDL, 
# recode grass_perennial as 176 which is Grassland/Pasture in CDL,
# recode hay and ryegrass as 37 which is Other Hay/Non Alfalfa
crop_map = {
    "corn": "1", "cotton": "2", "rice": "3", "sorghum": "4", "soybean": "5", "sunflower": "6",
    "peanut": "10", "tobacco": "11", "sweet_corn": "12", "pop_or_orn_corn": "13", "mint": "14",
    "barley": "21", "wheat_durum": "22", "wheat_spring": "23", "wheat_winter": "24", "rye": "27", "oat": "28",
    "millet": "29", "speltz": "30", "canola": "31", "flaxseed": "32", "safflower": "33",
    "mustard": "35", "alfalfa": "36", "hay": "37", "ryegrass": "37", "camelina": "38", "buckwheat": "39",
    "sugar_beet": "41", "dry_bean": "42", "potato": "43", "other": "44", "sugarcane": "45",
    "vegetable": "47", "cucumber": "50", "pea": "53", "tomato": "54", "herb": "57",
    "clover": "58", "sod_grass": "59", "fallow": "61",
    "cherry": "66", "peach": "67", "apple": "68", "grape": "69", "walnut": "76",
    "forest_deciduous": "141", "forest_evergreen": "142", "evergreen": "142", "shrub": "152",
    "pistachio": "204", "triticale": "205", "carrot": "206", "pepper": "216", "greens": "219", "squash": "222",
    "pumpkin": "229", "cabbage": "243", "grass_perennial": "176", "turnip": "247", "non_cropland": "999"
}

df_wide.replace(crop_map, inplace=True)

# Save to csv
df_wide.to_csv("data/edited/Regrow/OH_main_crop_wide_coded.csv", index=False)

