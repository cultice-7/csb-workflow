import pandas as pd
import os

# Import parameters from Snakemake
CSB_input_folder = snakemake.params.CSB_input_dir
CSB_output_folder = snakemake.params.CSB_output_dir
ag_census_input_folder = snakemake.params.ag_census_input_dir
states = snakemake.params.states
CSB_years = snakemake.params.CSB_years


def clean_ag_census_column_name(name):
    # Drop bracket characters, keeping the words inside as-is
    name = name.replace('(', '').replace(')', '')

    # Normalize remaining separators
    name = name.replace(', ', '_')
    name = name.replace(': ', '_ ')
    name = name.replace(' - ', '_')

    return name


def add_county_state_name(df, state):
    # Title-case all string columns (raw USDA NASS data is upper case)
    str_cols = df.select_dtypes(include=['object', 'string']).columns
    df[str_cols] = df[str_cols].apply(lambda col: col.str.title())
    # Build the county_state_name key, matching the CSB-side key (county_name + state)
    df['county_state_name'] = df['County Name'] + '_' + state
    return df


# Extract county-level tillage and cover-crop practice acreage from the Conservation Practices census file
def assign_census_conservation_practices(ag_census_input_folder, state):
    census = pd.read_csv(os.path.join(ag_census_input_folder, f"{state}_Census_Economics_Farm Operations_Farm Ownership_Conservation Practices_2017&2022.csv"))

    # Clean numeric values (strip thousands separators) and build a readable county name
    census['Value'] = pd.to_numeric(
        census['Value'].astype(str).str.replace(',', ''),
        errors='coerce'
    )
    census['County Name'] = census['County'] + ' ' + census['Geo Level']

    # Identify the rows for each practice (no-till, reduced till, conventional till, cover crop)
    measure = 'ACRES'  # 'ACRES' or 'NUMBER OF OPERATIONS' or 'MEASURED IN ACRES / OPERATION'
    mask_2017 = census['Year'] == 2017
    mask_2022 = census['Year'] == 2022
    mask_no_till = (census['Data Item'] == f'PRACTICES, LAND USE, CROPLAND, CONSERVATION TILLAGE, NO-TILL - {measure}'
                   ) & (census['Domain'] == 'TOTAL') & (census['Domain Category'] == 'NOT SPECIFIED')
    mask_reduced_till = (census['Data Item'] == f'PRACTICES, LAND USE, CROPLAND, CONSERVATION TILLAGE, (EXCL NO-TILL) - {measure}'
                   ) & (census['Domain'] == 'TOTAL') & (census['Domain Category'] == 'NOT SPECIFIED')
    mask_conventional_till = (census['Data Item'] == f'PRACTICES, LAND USE, CROPLAND, CONVENTIONAL TILLAGE - {measure}'
                   ) & (census['Domain'] == 'TOTAL') & (census['Domain Category'] == 'NOT SPECIFIED')
    mask_cover_crop = (census['Data Item'] == f'PRACTICES, LAND USE, CROPLAND, COVER CROP PLANTED, (EXCL CRP) - {measure}'
                   ) & (census['Domain'] == 'TOTAL') & (census['Domain Category'] == 'NOT SPECIFIED')

    # Filter and sort each practice's rows separately, then stack them into one table
    no_till = census[(mask_2017 | mask_2022) & mask_no_till].sort_values(by=['Year', 'County']).reset_index(drop=True)
    reduced_till = census[(mask_2017 | mask_2022) & mask_reduced_till].sort_values(by=['Year', 'County']).reset_index(drop=True)
    conventional_till = census[(mask_2017 | mask_2022) & mask_conventional_till].sort_values(by=['Year', 'County']).reset_index(drop=True)
    cover_crop = census[(mask_2017 | mask_2022) & mask_cover_crop].sort_values(by=['Year', 'County']).reset_index(drop=True)

    conservation_practices = pd.concat([no_till, reduced_till, conventional_till, cover_crop], ignore_index=True)

    return add_county_state_name(conservation_practices, state)


# Extract county-level harvested acreage for corn, soybeans, and winter wheat from the Field Crops census file
def assign_census_field_crops(ag_census_input_folder, state):
    census = pd.read_csv(os.path.join(ag_census_input_folder, f"{state}_Census_Crops_Field Crops_2017&2022.csv"))

    # Clean numeric values and build a readable county name
    census['Value'] = pd.to_numeric(
        census['Value'].astype(str).str.replace(',', ''),
        errors='coerce'
    )
    census['County Name'] = census['County'] + ' ' + census['Geo Level']

    # Identify the rows for each crop
    measure = 'ACRES HARVESTED'  # 'ACRES HARVESTED' or 'OPERATIONS WITH AREA HARVESTED' or 'PRODUCTION, MEASURED IN BU'
    mask_2017 = census['Year'] == 2017
    mask_2022 = census['Year'] == 2022
    mask_corn_grain = (census['Data Item'] == f'CORN, GRAIN - {measure}'
                   ) & (census['Domain'] == 'TOTAL') & (census['Domain Category'] == 'NOT SPECIFIED')
    mask_corn_silage = (census['Data Item'] == f'CORN, SILAGE - {measure}'
                   ) & (census['Domain'] == 'TOTAL') & (census['Domain Category'] == 'NOT SPECIFIED')
    mask_soybeans = (census['Data Item'] == f'SOYBEANS - {measure}'
                   ) & (census['Domain'] == 'TOTAL') & (census['Domain Category'] == 'NOT SPECIFIED')
    mask_winter_wheat = (census['Data Item'] == f'WHEAT, WINTER - {measure}'
                   ) & (census['Domain'] == 'TOTAL') & (census['Domain Category'] == 'NOT SPECIFIED')

    # Filter and sort each crop's rows separately, then stack them into one table
    corn_grain = census[(mask_2017 | mask_2022) & mask_corn_grain].sort_values(by=['Year', 'County']).reset_index(drop=True)
    corn_silage = census[(mask_2017 | mask_2022) & mask_corn_silage].sort_values(by=['Year', 'County']).reset_index(drop=True)
    soybeans = census[(mask_2017 | mask_2022) & mask_soybeans].sort_values(by=['Year', 'County']).reset_index(drop=True)
    winter_wheat = census[(mask_2017 | mask_2022) & mask_winter_wheat].sort_values(by=['Year', 'County']).reset_index(drop=True)

    field_crops = pd.concat([corn_grain, corn_silage, soybeans, winter_wheat], ignore_index=True)

    return add_county_state_name(field_crops, state)


# Extract county-level farmland/cropland/pastureland/woodland acreage from the Farm Operations census file
def assign_census_land_use(ag_census_input_folder, state):
    census = pd.read_csv(os.path.join(ag_census_input_folder, f"{state}_Census_Economics_Farm Operations_Farm Ownership_Conservation Practices_2017&2022.csv"))

    # Clean numeric values and build a readable county name
    census['Value'] = pd.to_numeric(
        census['Value'].astype(str).str.replace(',', ''),
        errors='coerce'
    )
    census['County Name'] = census['County'] + ' ' + census['Geo Level']

    mask_2017 = census['Year'] == 2017
    mask_2022 = census['Year'] == 2022

    # Total farm operations acreage
    measure = 'ACRES OPERATED'  # 'ACRES OPERATED' or 'NUMBER OF OPERATIONS'
    mask_farmland = (census['Data Item'] == f'FARM OPERATIONS - {measure}'
                   ) & (census['Domain'] == 'TOTAL') & (census['Domain Category'] == 'NOT SPECIFIED')

    # Land-use breakdown: cropland, cropland harvested, pastureland, woodland, other land
    measure = 'ACRES'  # 'ACRES' or 'NUMBER OF OPERATIONS'
    mask_cropland = (census['Data Item'] == f'AG LAND, CROPLAND - {measure}'
                   ) & (census['Domain'] == 'TOTAL') & (census['Domain Category'] == 'NOT SPECIFIED')
    mask_cropland_harvested = (census['Data Item'] == f'AG LAND, CROPLAND, HARVESTED - {measure}'
                   ) & (census['Domain'] == 'TOTAL') & (census['Domain Category'] == 'NOT SPECIFIED')
    mask_pastureland = (census['Data Item'] == f'AG LAND, PASTURELAND - {measure}'
                   ) & (census['Domain'] == 'TOTAL') & (census['Domain Category'] == 'NOT SPECIFIED')
    mask_pastureland_excl_cropland_woodland = (census['Data Item'] == f'AG LAND, PASTURELAND, (EXCL CROPLAND & WOODLAND) - {measure}'
                   ) & (census['Domain'] == 'TOTAL') & (census['Domain Category'] == 'NOT SPECIFIED')
    mask_woodland = (census['Data Item'] == f'AG LAND, WOODLAND - {measure}'
                   ) & (census['Domain'] == 'TOTAL') & (census['Domain Category'] == 'NOT SPECIFIED')
    mask_otherland = (census['Data Item'] == f'AG LAND, (EXCL CROPLAND & PASTURELAND & WOODLAND) - {measure}'
                   ) & (census['Domain'] == 'TOTAL') & (census['Domain Category'] == 'NOT SPECIFIED')

    # Farmland ownership split: owned vs. rented
    measure = 'ACRES'  # 'ACRES' or 'NUMBER OF OPERATIONS'
    mask_farmland_owned = (census['Data Item'] == f'AG LAND, OWNED, IN FARMS - {measure}'
                   ) & (census['Domain'] == 'TOTAL') & (census['Domain Category'] == 'NOT SPECIFIED')
    mask_farmland_rented = (census['Data Item'] == f'AG LAND, RENTED FROM OTHERS, IN FARMS - {measure}'
                   ) & (census['Domain'] == 'TOTAL') & (census['Domain Category'] == 'NOT SPECIFIED')

    # Filter and sort each land-use category's rows separately, then stack them into one table
    farmland = census[(mask_2017 | mask_2022) & mask_farmland].sort_values(by=['Year', 'County']).reset_index(drop=True)
    cropland = census[(mask_2017 | mask_2022) & mask_cropland].sort_values(by=['Year', 'County']).reset_index(drop=True)
    cropland_harvested = census[(mask_2017 | mask_2022) & mask_cropland_harvested].sort_values(by=['Year', 'County']).reset_index(drop=True)
    pastureland = census[(mask_2017 | mask_2022) & mask_pastureland].sort_values(by=['Year', 'County']).reset_index(drop=True)
    pastureland_excl_cropland_woodland = census[(mask_2017 | mask_2022) & mask_pastureland_excl_cropland_woodland].sort_values(by=['Year', 'County']).reset_index(drop=True)
    woodland = census[(mask_2017 | mask_2022) & mask_woodland].sort_values(by=['Year', 'County']).reset_index(drop=True)
    otherland = census[(mask_2017 | mask_2022) & mask_otherland].sort_values(by=['Year', 'County']).reset_index(drop=True)

    farmland_owned = census[(mask_2017 | mask_2022) & mask_farmland_owned].sort_values(by=['Year', 'County']).reset_index(drop=True)
    farmland_rented = census[(mask_2017 | mask_2022) & mask_farmland_rented].sort_values(by=['Year', 'County']).reset_index(drop=True)

    land_use = pd.concat([farmland, cropland, cropland_harvested,
                          pastureland, pastureland_excl_cropland_woodland,
                          woodland, otherland,
                          farmland_owned, farmland_rented], ignore_index=True)

    return add_county_state_name(land_use, state)


# Extract county-level farmland and cropland acreage by ownership tenure from the Land Ownership census file
def assign_census_land_ownership(ag_census_input_folder, state):
    census = pd.read_csv(os.path.join(ag_census_input_folder, f"{state}_Census_Demographics_Farms_Land Ownership_Assets_2017&2022.csv"))

    # Clean numeric values and build a readable county name
    census['Value'] = pd.to_numeric(
        census['Value'].astype(str).str.replace(',', ''),
        errors='coerce'
    )
    census['County Name'] = census['County'] + ' ' + census['Geo Level']

    mask_2017 = census['Year'] == 2017
    mask_2022 = census['Year'] == 2022

    # Tenure breakdown (full owner / part owner / tenant, etc.) for farmland and harvested cropland
    measure_farmland = 'ACRES OPERATED'  # 'ACRES OPERATED' or 'NUMBER OF OPERATIONS'
    measure_cropland = 'ACRES'  # 'ACRES' or 'NUMBER OF OPERATIONS'

    mask_farmland_ownership = (census['Data Item'] == f'FARM OPERATIONS - {measure_farmland}'
                   ) & (census['Domain'] == 'TENURE')
    mask_cropland_ownership = (census['Data Item'] == f'AG LAND, CROPLAND, HARVESTED - {measure_cropland}'
                   ) & (census['Domain'] == 'TENURE')

    # Filter and sort each tenure category's rows separately (kept by Domain Category, e.g. full owner/tenant, not collapsed)
    farmland_ownership = census[(mask_2017 | mask_2022) & mask_farmland_ownership].sort_values(by=['Year', 'County', 'Domain Category']).reset_index(drop=True)
    cropland_ownership = census[(mask_2017 | mask_2022) & mask_cropland_ownership].sort_values(by=['Year', 'County', 'Domain Category']).reset_index(drop=True)

    land_ownership = pd.concat([farmland_ownership, cropland_ownership], ignore_index=True)
    # Fold the tenure category into the Data Item label so each tenure group gets its own output column after pivoting
    land_ownership['Data Item'] = land_ownership['Data Item'] + ' ' + land_ownership['Domain Category']

    return add_county_state_name(land_ownership, state)


# Extract county-level producer demographic and farm-decisionmaking characteristics from the Producer Characteristics census file
def assign_census_farmer_characteristics(ag_census_input_folder, state):
    census = pd.read_csv(os.path.join(ag_census_input_folder, f"{state}_Census_Demographics_Producer Characteristics_2017&2022.csv"))

    # Clean numeric values and build a readable county name
    census['Value'] = pd.to_numeric(
        census['Value'].astype(str).str.replace(',', ''),
        errors='coerce'
    )
    census['County Name'] = census['County'] + ' ' + census['Geo Level']

    mask_2017 = census['Year'] == 2017
    mask_2022 = census['Year'] == 2022

    measure = 'NUMBER OF PRODUCERS'  # 'NUMBER OF PRODUCERS' or 'ACRES OPERATED' or 'NUMBER OF OPERATIONS'

    # Overall producer count
    mask_producer_number = (census['Data Item'] == f'PRODUCERS - {measure}'
                        ) & (census['Domain'] == 'TOTAL') & (census['Domain Category'] == 'NOT SPECIFIED')

    # Producer sex
    mask_producer_sex = ((census['Data Item'] == f'PRODUCERS, FEMALE - {measure}') |
                         (census['Data Item'] == f'PRODUCERS, MALE - {measure}')
                        ) & (census['Domain'] == 'TOTAL') & (census['Domain Category'] == 'NOT SPECIFIED')

    # Producer age brackets plus average age
    mask_producer_age = ((census['Data Item'] == f'PRODUCERS, AGE LT 25 - {measure}') |
                         (census['Data Item'] == f'PRODUCERS, AGE 25 TO 34 - {measure}') |
                         (census['Data Item'] == f'PRODUCERS, AGE 35 TO 44 - {measure}') |
                         (census['Data Item'] == f'PRODUCERS, AGE 45 TO 54 - {measure}') |
                         (census['Data Item'] == f'PRODUCERS, AGE 55 TO 64 - {measure}') |
                         (census['Data Item'] == f'PRODUCERS, AGE 65 TO 74 - {measure}') |
                         (census['Data Item'] == f'PRODUCERS, AGE GE 75 - {measure}') |
                         (census['Data Item'] == 'PRODUCERS - AGE, AVG, MEASURED IN YEARS')
                        ) & (census['Domain'] == 'TOTAL') & (census['Domain Category'] == 'NOT SPECIFIED')

    # Producer race
    mask_producer_race = ((census['Data Item'] == f'PRODUCERS, AMERICAN INDIAN OR ALASKA NATIVE - {measure}') |
                          (census['Data Item'] == f'PRODUCERS, ASIAN - {measure}') |
                          (census['Data Item'] == f'PRODUCERS, BLACK OR AFRICAN AMERICAN - {measure}') |
                          (census['Data Item'] == f'PRODUCERS, NATIVE HAWAIIAN OR OTHER PACIFIC ISLANDER - {measure}') |
                          (census['Data Item'] == f'PRODUCERS, WHITE - {measure}') |
                          (census['Data Item'] == f'PRODUCERS, MULTI-RACE - {measure}')
                         ) & (census['Domain'] == 'TOTAL') & (census['Domain Category'] == 'NOT SPECIFIED')

    # Producer ethnicity (Hispanic)
    mask_producer_hispanic = (census['Data Item'] == f'PRODUCERS, HISPANIC - {measure}'
                         ) & (census['Domain'] == 'TOTAL') & (census['Domain Category'] == 'NOT SPECIFIED')

    # Household size
    mask_producer_household = (census['Data Item'] == 'PRODUCERS - PERSONS IN HOUSEHOLD, MEASURED IN PERSONS'
                         ) & (census['Domain'] == 'TOTAL') & (census['Domain Category'] == 'NOT SPECIFIED')

    # Farming experience (years on any operation / present operation)
    mask_producer_experience = ((census['Data Item'] == 'PRODUCERS - YEARS ON ANY OPERATION, AVG, MEASURED IN YEARS') |
                                (census['Data Item'] == f'PRODUCERS, YEARS ON ANY OPERATION, LT 6 YEARS - {measure}') |
                                (census['Data Item'] == f'PRODUCERS, YEARS ON ANY OPERATION, 6 TO 10 YEARS - {measure}') |
                                (census['Data Item'] == f'PRODUCERS, YEARS ON ANY OPERATION, GE 11 YEARS - {measure}') |
                                (census['Data Item'] == 'PRODUCERS - YEARS ON PRESENT OPERATION, AVG, MEASURED IN YEARS')
                               ) & (census['Domain'] == 'TOTAL') & (census['Domain Category'] == 'NOT SPECIFIED')

    # Which farm decisions the producer is involved in
    mask_producer_decisionmaking = ((census['Data Item'] == f'PRODUCERS, DAY TO DAY DECISIONMAKING - {measure}') |
                          (census['Data Item'] == f'PRODUCERS, LAND USE OR CROP DECISIONMAKING - {measure}') |
                          (census['Data Item'] == f'PRODUCERS, LIVESTOCK DECISIONMAKING - {measure}') |
                          (census['Data Item'] == f'PRODUCERS, MARKETING DECISIONMAKING - {measure}') |
                          (census['Data Item'] == f'PRODUCERS, RECORD KEEPING OR FINANCIAL MGMT DECISIONMAKING - {measure}') |
                          (census['Data Item'] == f'PRODUCERS, ESTATE OR SUCCESSION PLANNING DECISIONMAKING - {measure}')
                         ) & (census['Domain'] == 'TOTAL') & (census['Domain Category'] == 'NOT SPECIFIED')

    # Whether the producer lives on the operation
    mask_producer_residence = ((census['Data Item'] == f'PRODUCERS, RESIDENCE, ON OPERATION - {measure}') |
                                (census['Data Item'] == f'PRODUCERS, RESIDENCE, NOT ON OPERATION - {measure}')
                               ) & (census['Domain'] == 'TOTAL') & (census['Domain Category'] == 'NOT SPECIFIED')

    # Military service history
    mask_producer_military = ((census['Data Item'] == f'PRODUCERS, MILITARY SERVICE, ACTIVE DUTY NOW OR IN THE PAST - {measure}') |
                                (census['Data Item'] == f'PRODUCERS, MILITARY SERVICE, NEVER SERVED OR ONLY ON ACTIVE DUTY FOR TRAINING IN RESERVES OR NATIONAL GUARD - {measure}')
                               ) & (census['Domain'] == 'TOTAL') & (census['Domain Category'] == 'NOT SPECIFIED')

    # Primary occupation (farming vs. other)
    mask_producer_occupation = ((census['Data Item'] == f'PRODUCERS, PRIMARY OCCUPATION, FARMING - {measure}') |
                                (census['Data Item'] == f'PRODUCERS, PRIMARY OCCUPATION, (EXCL FARMING) - {measure}')
                               ) & (census['Domain'] == 'TOTAL') & (census['Domain Category'] == 'NOT SPECIFIED')

    # Filter and sort each characteristic's rows separately, then stack them into one table
    producer_number = census[(mask_2017 | mask_2022) & mask_producer_number].sort_values(by=['Year', 'County', 'Data Item']).reset_index(drop=True)
    producer_sex = census[(mask_2017 | mask_2022) & mask_producer_sex].sort_values(by=['Year', 'County', 'Data Item']).reset_index(drop=True)
    producer_age = census[(mask_2017 | mask_2022) & mask_producer_age].sort_values(by=['Year', 'County', 'Data Item']).reset_index(drop=True)
    producer_race = census[(mask_2017 | mask_2022) & mask_producer_race].sort_values(by=['Year', 'County', 'Data Item']).reset_index(drop=True)
    producer_hispanic = census[(mask_2017 | mask_2022) & mask_producer_hispanic].sort_values(by=['Year', 'County', 'Data Item']).reset_index(drop=True)
    producer_household = census[(mask_2017 | mask_2022) & mask_producer_household].sort_values(by=['Year', 'County', 'Data Item']).reset_index(drop=True)
    producer_experience = census[(mask_2017 | mask_2022) & mask_producer_experience].sort_values(by=['Year', 'County', 'Data Item']).reset_index(drop=True)
    producer_decisionmaking = census[(mask_2017 | mask_2022) & mask_producer_decisionmaking].sort_values(by=['Year', 'County', 'Data Item']).reset_index(drop=True)
    producer_residence = census[(mask_2017 | mask_2022) & mask_producer_residence].sort_values(by=['Year', 'County', 'Data Item']).reset_index(drop=True)
    producer_military = census[(mask_2017 | mask_2022) & mask_producer_military].sort_values(by=['Year', 'County', 'Data Item']).reset_index(drop=True)
    producer_occupation = census[(mask_2017 | mask_2022) & mask_producer_occupation].sort_values(by=['Year', 'County', 'Data Item']).reset_index(drop=True)

    farmer_characteristics = pd.concat([producer_number, producer_sex, producer_age, producer_race,
                                       producer_hispanic, producer_household, producer_experience, producer_decisionmaking,
                                       producer_residence, producer_military, producer_occupation], ignore_index=True)

    return add_county_state_name(farmer_characteristics, state)


# =============================================================================
# Main script: for each CSB year-window and state, build the county-level USDA
# Ag Census table, reshape it to wide format, and attach it to the sorted list
# of counties present among that state's CSB fields (county-level output,
# keyed on county_state_name -- not merged onto individual fields)
# =============================================================================
for year in CSB_years:
    for state in states:

        # Build the 5 category tables for this state
        land_use = assign_census_land_use(ag_census_input_folder, state)
        field_crops = assign_census_field_crops(ag_census_input_folder, state)
        conservation_practices = assign_census_conservation_practices(ag_census_input_folder, state)
        land_ownership = assign_census_land_ownership(ag_census_input_folder, state)
        farmer_characteristics = assign_census_farmer_characteristics(ag_census_input_folder, state)

        # Stack all 5 county-level census categories into one long table
        census_county_data = pd.concat([land_use, land_ownership, field_crops,
                                conservation_practices, farmer_characteristics], ignore_index=True)

        # Reshape from long to wide: one row per county_state_name, one column per Data Item x Year
        census_county_data_wide = census_county_data.pivot(index='county_state_name',
                                                           columns=['Year', 'Data Item'],
                                                           values='Value')

        census_county_data_wide.columns = [
            clean_ag_census_column_name(f'{var}_{yr}') for yr, var in census_county_data_wide.columns
        ]

        # Build the unique, alphabetically-sorted list of counties present among this state's CSB fields
        # Supplement_9 is a county-level dataset keyed on county_state_name
        CSB_input_path = os.path.join(CSB_input_folder, f"{state}_CSB{year}_census_tract_table.parquet")
        CSB_supplement_1 = pd.read_parquet(CSB_input_path)

        # Drop counties belonging to a different state (can happen when a border field's nearest-tract
        # fallback in supplement_1 matches a neighboring state's tract instead of its own state's)
        counties = CSB_supplement_1[['county_name', 'state_name', 'county_state_name']].drop_duplicates()
        counties = counties[counties['state_name'] == state]
        counties = counties.sort_values('county_state_name').reset_index(drop=True)

        # Merge county-level census data onto the sorted county list
        CSB_supplement_9 = counties.merge(census_county_data_wide, on='county_state_name', how='left')

        # Convert float64 columns to float32 to save memory
        float64_cols = CSB_supplement_9.select_dtypes(include=["float64"]).columns
        CSB_supplement_9[float64_cols] = CSB_supplement_9[float64_cols].astype("float32")

        # Save attribute table
        output_path_table = os.path.join(CSB_output_folder, f"{state}_CSB{year}_ag_census_table.parquet")
        CSB_supplement_9.to_parquet(output_path_table, index=False, compression="zstd")
        print(f'Supplementary data 9 for {state} CSB{year} is created and saved')
