import geopandas as gpd
import os

#---# Load required datasets
# State and country boundaries
tract_boundaries = gpd.read_file("data/Census/census_tract/cb_2023_us_tract_500k.shp")
tract_boundaries = tract_boundaries[['STATEFP', 'STUSPS', 'COUNTYFP', 'NAMELSADCO', 'GEOID', 'ALAND', 'AWATER', 'geometry']]
tract_boundaries.rename(columns={
    'STATEFP': 'state_id', 'STUSPS': 'state_name', 'COUNTYFP':'county_id', 'NAMELSADCO': 'county_name', 'GEOID': 'census_tract_id', 'ALAND':'tract_land_area', 'AWATER':'tract_water_area'}, inplace=True)

# Input and output folders for CSB
input_folder_CSB = "data/edited/CSB/"
output_folder_CSB = "data/edited/CSB/"

# List of states
states = snakemake.params.states


for state in states:
    
    input_path_CSB = os.path.join(input_folder_CSB, f"{state}_CSB1724_clipped.gpkg")

    # Load CSB data
    CSB_shape = gpd.read_file(input_path_CSB)

    # Setting active geometry column
    CSB_shape = CSB_shape.set_geometry('geometry')
    
    # Reproject geometry to an equal-area CRS (NAD83/CONUS Albers)
    CSB_shape = CSB_shape.to_crs(epsg=5070)
    
    # Keep only the columns necessary for the spatial join
    cols_to_keep = ['CSBID', 'geometry']
    CSB_shape = CSB_shape[cols_to_keep]


    #---#  Add Census data: state, county, tract
    # Set paths to output files
    output_path_geojson = os.path.join(output_folder_CSB, f"{state}_CSB1724_supplement_1_spatial.geojson")
    output_path_csv = os.path.join(output_folder_CSB, f"{state}_CSB1724_supplement_1_table.csv")
    try:
        #---# Polygon–polygon intersections (too time-consuming)
        #tract_boundaries = tract_boundaries.to_crs(epsg=5070)
        #intersections_location = gpd.overlay(CSB_shape, tract_boundaries, how='intersection')
        #intersections_location['overlap_area_temp'] = intersections_location.geometry.area
        #largest_overlap_location = intersections_location.sort_values('overlap_area_temp', ascending=False).drop_duplicates('CSBID')
        #columns_to_merge_location = ['CSBID', 'state_id', 'state_name', 'county_id', 'county_name', 'census_tract_id', 'tract_land_area', 'tract_water_area']
        #CSB_shape = CSB_shape.merge(largest_overlap_location[columns_to_merge_location], on='CSBID', how='left', suffixes=('', '_temp'))
        #cols_to_drop = [col for col in CSB_shape.columns if col.endswith('_temp')]
        #CSB_shape.drop(columns=cols_to_drop, inplace=True)
        
        #---# Centroid point-in-polygon spatial joins (faster tool)
        # Ensure dataset is in the same projected CRS
        tract_boundaries = tract_boundaries.to_crs(epsg=5070)
        
        # Compute centroids for fields (faster than polygon intersections)
        parcel_centroids = CSB_shape.copy()
        parcel_centroids["geometry"] = parcel_centroids.geometry.centroid

        # Spatial join: assign each field centroid to the tract it falls within
        joined = gpd.sjoin(parcel_centroids, tract_boundaries, how="left", predicate="within")
        
        # Identify unmatched fields
        unmatched_mask = joined['census_tract_id'].isna()
        
        # For unmatched fields, find the nearest existing Census tract
        nearest_tract = gpd.sjoin_nearest(joined.loc[unmatched_mask, ['CSBID', 'geometry']], tract_boundaries, how='left', distance_col='distance')
        
        # Align the indices of the joined and nearest-tract dataframes to fill rows only for unmatched fields
        nearest_tract.index = joined.loc[unmatched_mask].index

        # Replace rows without assigned tract data in joined
        joined.update(nearest_tract, overwrite=False)

        # Columns you want to bring back from the join
        cols_to_merge = ["CSBID", "state_id", "state_name", "county_id", "county_name", "census_tract_id", "tract_land_area", "tract_water_area"]

        # Keep only existing columns
        cols_available = [col for col in cols_to_merge if col in joined.columns]

        # Merge the tract attributes back to the original polygons
        CSB_shape = CSB_shape.merge(joined[cols_available], on="CSBID", how="left", suffixes=("", "_temp"))

        # Drop temp duplicates if any (from suffixes) and the spatial join index
        cols_to_drop = [col for col in CSB_shape.columns if col.endswith("_temp")]
        CSB_shape.drop(columns=cols_to_drop, inplace=True)
        print("Tract data added to attribute table")

    except Exception as e:
        print(f"Error processing slope: {e}")
        raise
    
    # Check whether there are any missing values
    cols_to_check = ['CSBID', 'state_id', 'county_id', 'census_tract_id']
    if CSB_shape[cols_to_check].isna().any().any():
        print(CSB_shape[CSB_shape[cols_to_check].isna().any(axis=1)])
    
    # Save geojson and csv files
    CSB_shape.to_file(output_path_geojson, driver="GeoJSON")
    attribute_table = CSB_shape.drop(columns='geometry')
    print(attribute_table.shape) #check df shape
    attribute_table.to_csv(output_path_csv, index=False)
    print(f'Supplementary data 1 for {state} is created and saved')