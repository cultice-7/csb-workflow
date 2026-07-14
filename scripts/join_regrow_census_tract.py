import geopandas as gpd
import os

# Import parameters from Snakemake
regrow_input_folder = snakemake.params.regrow_input_dir
regrow_output_folder = snakemake.params.regrow_output_dir
census_tract_input_folder = snakemake.params.census_tract_input_dir
states = snakemake.params.states
target_CRS = snakemake.params.target_CRS


#---# Load required datasets
# Census tract boundaries
tract_boundaries_input_path = os.path.join(census_tract_input_folder, "cb_2023_us_tract_500k.shp")
tract_boundaries = gpd.read_file(tract_boundaries_input_path)
tract_boundaries = tract_boundaries[['STATEFP', 'STUSPS', 'COUNTYFP', 'NAMELSADCO', 'GEOID', 'ALAND', 'AWATER', 'geometry']]
tract_boundaries.rename(columns={'STATEFP': 'state_id', 'STUSPS': 'state_name', 
                                 'COUNTYFP':'county_id', 'NAMELSADCO': 'county_name', 
                                 'GEOID': 'census_tract_id', 'ALAND':'tract_land_area', 
                                 'AWATER':'tract_water_area'}, inplace=True)

for state in states:
    
    regrow_input_path = os.path.join(regrow_input_folder, f"{state}_regrow_fieldID_geometry.parquet")
    
    # Load Regrow data
    regrow_shape = gpd.read_parquet(regrow_input_path)

    # Setting active geometry column
    regrow_shape = regrow_shape.set_geometry('geometry')
    
    # Keep only the columns necessary for the spatial join
    cols_to_keep = ['field_id', 'geometry']
    regrow_shape = regrow_shape[cols_to_keep]

    #---# Add Census data: state, county, tract
    # Set paths to output files
    output_path_spatial = os.path.join(regrow_output_folder, f"{state}_regrow_census_tract_spatial.parquet")
    output_path_table = os.path.join(regrow_output_folder, f"{state}_regrow_census_tract_table.parquet")
    try:
        #---# Polygon–polygon intersections (too time-consuming)
        #tract_boundaries = tract_boundaries.to_crs(epsg=5070)
        #intersections_location = gpd.overlay(regrow_shape, tract_boundaries, how='intersection')
        #intersections_location['overlap_area_temp'] = intersections_location.geometry.area
        #largest_overlap_location = intersections_location.sort_values('overlap_area_temp', ascending=False).drop_duplicates('field_id')
        #columns_to_merge_location = ['field_id', 'state_id', 'state_name', 'county_id', 'county_name', 'census_tract_id', 'tract_land_area', 'tract_water_area']
        #regrow_shape = regrow_shape.merge(largest_overlap_location[columns_to_merge_location], on='field_id', how='left', suffixes=('', '_temp'))
        #cols_to_drop = [col for col in regrow_shape.columns if col.endswith('_temp')]
        #regrow_shape.drop(columns=cols_to_drop, inplace=True)
        
        #---# Centroid point-in-polygon spatial joins (faster tool)
        # Ensure both datasets are in the same projected CRS
        tract_boundaries = tract_boundaries.to_crs(target_CRS)

        # Build the county_state_name variable, used as a key to join county-level data (e.g. supplement 9/USDA Ag Census)
        tract_boundaries['county_state_name'] = tract_boundaries['county_name'] + '_' + tract_boundaries['state_name']
        
        # Compute centroids for fields (faster than polygon intersections)
        field_centroids = regrow_shape.copy()
        field_centroids["geometry"] = field_centroids.geometry.centroid

        # Spatial join: assign each field centroid to the tract it falls within
        joined = gpd.sjoin(field_centroids, tract_boundaries, how="left", predicate="within")

        # Guard against overlapping/invalid tract polygons producing more than one match per centroid
        assert not joined['field_id'].duplicated().any(), (
            "A field centroid matched more than one Census tract — check tract_boundaries for overlapping polygons"
        )

        # Identify unmatched fields
        unmatched_mask = joined['census_tract_id'].isna()

        # For unmatched fields, find the nearest existing Census tract
        nearest_tract = gpd.sjoin_nearest(joined.loc[unmatched_mask, ['field_id', 'geometry']], tract_boundaries, how='left', distance_col='distance')

        # Break ties deterministically if a centroid is exactly equidistant from more than one tract
        nearest_tract = nearest_tract.drop_duplicates(subset='field_id', keep='first')

        # Replace rows without assigned tract data in joined
        # Index both frames on field_id so .update() aligns by key, not row position
        joined = joined.set_index('field_id', drop=False)
        joined.update(nearest_tract.set_index('field_id', drop=False), overwrite=False)
        joined = joined.reset_index(drop=True)

        # Columns you want to bring back from the join
        cols_to_merge = ["field_id", "state_id", "state_name", "county_id", "county_name", "county_state_name", "census_tract_id", "tract_land_area", "tract_water_area"]

        # Keep only existing columns
        cols_available = [col for col in cols_to_merge if col in joined.columns]

        # Merge the tract attributes back to the original polygons
        regrow_shape = regrow_shape.merge(joined[cols_available], on="field_id", how="left", suffixes=("", "_temp"), validate="one_to_one")

        # Drop temp duplicates if any (from suffixes) and the spatial join index
        cols_to_drop = [col for col in regrow_shape.columns if col.endswith("_temp")]
        regrow_shape.drop(columns=cols_to_drop, inplace=True)
        print("Tract data added to attribute table")

    except Exception as e:
        print(f"Error processing slope: {e}")
        raise
    

    # Check whether there are any missing values
    cols_to_check = ['field_id', 'state_id', 'county_id', 'census_tract_id']
    if regrow_shape[cols_to_check].isna().any().any():
        print(regrow_shape[regrow_shape[cols_to_check].isna().any(axis=1)])
        
    # Convert float64 columns to float32 to save memory
    float64_cols = regrow_shape.select_dtypes(include=["float64"]).columns
    regrow_shape[float64_cols] = regrow_shape[float64_cols].astype("float32")
    
    #---# Save files w/ and w/o geometry
    #regrow_shape.to_parquet(output_path_spatial, compression="zstd")
    attribute_table = regrow_shape.drop(columns='geometry')
    print(attribute_table.shape) #check df shape
    attribute_table.to_parquet(output_path_table, index=False, compression="zstd")
    print(f'Supplementary data 1 for {state} is created and saved')