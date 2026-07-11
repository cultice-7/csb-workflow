import geopandas as gpd
from rasterstats import zonal_stats
import os

#---# Load required datasets
# Projected raster paths
elevation_proj_path = "data/Geo/elevation/elevation_reproj.tif"
slope_proj_path = "data/Geo/elevation/slope_reproj.tif"

# State and country boundaries
tract_boundaries = gpd.read_file("data/Census/census_tract/cb_2023_us_tract_500k.shp")
tract_boundaries = tract_boundaries[['STATEFP', 'STUSPS', 'COUNTYFP', 'NAMELSADCO', 'GEOID', 'ALAND', 'AWATER', 'geometry']]
tract_boundaries.rename(columns={
    'STATEFP': 'state_id', 'STUSPS': 'state_name', 'COUNTYFP':'county_id', 'NAMELSADCO': 'county_name', 'GEOID': 'census_tract_id', 'ALAND':'tract_land_area', 'AWATER':'tract_water_area'}, inplace=True)

# Watershed data
subbasin = gpd.read_file("data/Geo/watershed/subbasin.shp")
watershed = gpd.read_file("data/Geo/watershed/watershed.shp")
subwatershed = gpd.read_file("data/Geo/watershed/subwatershed.shp")

# Input and output folders for Regrow
input_folder_Regrow = "data/edited/Regrow/"
output_folder_Regrow = "data/edited/Regrow/"

states = snakemake.params.states


for state in states:
    
    input_path_Regrow = os.path.join(input_folder_Regrow, f"{state}_regrow_shape_table.geojson")
    
    # Load Regrow data
    regrow_shape = gpd.read_file(input_path_Regrow)

    # Setting active geometry column
    regrow_shape = regrow_shape.set_geometry('geometry')
    
    # Reproject geometry to an equal-area CRS (NAD83/CONUS Albers)
    regrow_shape = regrow_shape.to_crs(epsg=5070)
    
    # Keep only the columns necessary for the spatial join
    cols_to_keep = ['field_id', 'geometry']
    regrow_shape = regrow_shape[cols_to_keep]


    #---# Join supplementary data to main dataset
    # Add Census data: state, county, tract
    regrow_shape_copy = regrow_shape.copy()
    output_path_geojson = os.path.join(output_folder_Regrow, f"{state}_regrow_supplement_1_spatial.geojson")
    output_path_csv = os.path.join(output_folder_Regrow, f"{state}_regrow_supplement_1_table.csv")
    try:
        # Polygon–polygon intersections (too time-consuming)
        #tract_boundaries = tract_boundaries.to_crs(epsg=5070)
        #intersections_location = gpd.overlay(regrow_shape_copy, tract_boundaries, how='intersection')
        #intersections_location['overlap_area_temp'] = intersections_location.geometry.area
        #largest_overlap_location = intersections_location.sort_values('overlap_area_temp', ascending=False).drop_duplicates('field_id')
        #columns_to_merge_location = ['field_id', 'state_id', 'state_name', 'county_id', 'county_name', 'census_tract_id', 'tract_land_area', 'tract_water_area']
        #regrow_shape_copy = regrow_shape_copy.merge(largest_overlap_location[columns_to_merge_location], on='field_id', how='left', suffixes=('', '_temp'))
        #cols_to_drop = [col for col in regrow_shape_copy.columns if col.endswith('_temp')]
        #regrow_shape_copy.drop(columns=cols_to_drop, inplace=True)
        
        # Centroid point-in-polygon spatial joins (faster tool)
        # Ensure both datasets are in the same projected CRS
        tract_boundaries = tract_boundaries.to_crs(epsg=5070)
        
        # Compute centroids for parcels (faster than polygon intersections)
        parcel_centroids = regrow_shape_copy.copy()
        parcel_centroids["geometry"] = parcel_centroids.geometry.centroid

        # Spatial join: assign each parcel centroid to the tract it falls within
        joined = gpd.sjoin(parcel_centroids, tract_boundaries, how="left", predicate="within")

        # Columns you want to bring back from the join
        cols_to_merge = ["field_id", "state_id", "state_name", "county_id", "county_name", "census_tract_id", "tract_land_area", "tract_water_area"]

        # Keep only existing columns
        cols_available = [col for col in cols_to_merge if col in joined.columns]

        # Merge the tract attributes back to the original polygons
        regrow_shape_copy = regrow_shape_copy.merge(joined[cols_available], on="field_id", how="left", suffixes=("", "_temp"))

        # Drop temp duplicates if any (from suffixes) and the spatial join index
        cols_to_drop = [col for col in regrow_shape_copy.columns if col.endswith("_temp")]
        regrow_shape_copy.drop(columns=cols_to_drop, inplace=True)
        print("Tract data added to attribute table")

    except Exception as e:
        print(f"Error processing slope: {e}")
        raise
    
    # Save geojson and csv files
    regrow_shape_copy.to_file(output_path_geojson, driver="GeoJSON")
    attribute_table = regrow_shape_copy.drop(columns='geometry')
    print(attribute_table.shape) #check df shape
    attribute_table.to_csv(output_path_csv, index=False)
    print(f'Supplementary data 1 for {state} is created and saved')
    


    # Add zonal statistics for elevation and slope
    regrow_shape_copy = regrow_shape.copy()
    output_path_geojson = os.path.join(output_folder_Regrow, f"{state}_regrow_supplement_2_spatial.geojson")
    output_path_csv = os.path.join(output_folder_Regrow, f"{state}_regrow_supplement_2_table.csv")
    # Zonal statistics for elevation
    try:
        print("Calculating and adding mean elevation...")
        elevation_stats = zonal_stats(regrow_shape_copy, elevation_proj_path, stats="mean")
        elevation_means = [stat['mean'] for stat in elevation_stats]
        regrow_shape_copy['elevation_mean'] = elevation_means
        print("Mean elevation added to attribute table.")
    except Exception as e:
        print(f"Error processing elevation: {e}")
        raise

    # Zonal statistics for slope
    try:
        print("Calculating and adding mean slope...")
        slope_stats = zonal_stats(regrow_shape_copy, slope_proj_path, stats="mean")
        slope_means = [stat['mean'] for stat in slope_stats]
        regrow_shape_copy['slope_mean'] = slope_means
        print("Mean slope added to attribute table.")
    except Exception as e:
        print(f"Error processing slope: {e}")
        raise
    
    # Save geojson and csv files
    #regrow_shape_copy.to_file(output_path_geojson, driver="GeoJSON")
    attribute_table = regrow_shape_copy.drop(columns='geometry')
    print(attribute_table.shape) #check df shape
    attribute_table.to_csv(output_path_csv, index=False)
    print(f'Supplementary data 2 for {state} is created and saved')



    # Add subbasin, watershed and subwatershed
    regrow_shape_copy = regrow_shape.copy()
    output_path_geojson = os.path.join(output_folder_Regrow, f"{state}_regrow_supplement_3_spatial.geojson")
    output_path_csv = os.path.join(output_folder_Regrow, f"{state}_regrow_supplement_3_table.csv")
    print("Adding hydrography data...")
    hu_codes = ['8', '10', '12']
    for code in hu_codes:
        if (code == '8'): 
            df_hu = subbasin.copy()
            df_hu = df_hu.drop_duplicates(subset=[f'huc{code}', 'geometry']).reset_index(drop=True)
            cols = [f'huc{code}', 'name', 'geometry']
            df_hu = df_hu[cols]
            df_hu.rename(columns={f'huc{code}':'subbasin_id', 'name': 'subbasin_name'}, inplace = True)
        elif (code == '10'): 
            df_hu = watershed.copy()
            df_hu = df_hu.drop_duplicates(subset=[f'huc{code}', 'geometry']).reset_index(drop=True)
            cols = [f'huc{code}', 'name', 'hutype', 'geometry']
            df_hu = df_hu[cols]
            df_hu.rename(columns={f'huc{code}':'watershed_id', 'name': 'watershed_name', 'hutype': 'watershed_type'}, inplace = True)
        elif (code == '12'): 
            df_hu = subwatershed.copy()
            df_hu = df_hu.drop_duplicates(subset=[f'huc{code}', 'geometry']).reset_index(drop=True)
            cols = [f'huc{code}', 'name', 'hutype', 'geometry']
            df_hu = df_hu[cols]
            df_hu.rename(columns={f'huc{code}':'subwatershed_id', 'name': 'subwatershed_name', 'hutype': 'subwatershed_type'}, inplace = True)
        
        try:
            # Polygon–polygon intersections (too time-consuming)
            #df_hu = df_hu.to_crs(epsg=5070) 
            #intersections_hydro = gpd.overlay(regrow_shape_copy, df_hu, how='intersection')
            #intersections_hydro['overlap_area_temp'] = intersections_hydro.geometry.area
            #largest_overlap_hydro = intersections_hydro.sort_values('overlap_area_temp', ascending=False).drop_duplicates('field_id')
            #columns_to_merge_hydro = largest_overlap_hydro.columns[largest_overlap_hydro.columns.str.contains('field_id|_id|_name|_type', regex=True)]
            #regrow_shape_copy = regrow_shape_copy.merge(largest_overlap_hydro[columns_to_merge_hydro], on='field_id', how='left', suffixes=('', '_temp'))
            #cols_to_drop = [col for col in regrow_shape_copy.columns if col.endswith('_temp')]
            #regrow_shape_copy.drop(columns=cols_to_drop, inplace=True)
            #print("Hydrography data added to attribute table.")
            
            # Centroid point-in-polygon spatial joins (faster tool)
            # Reproject watershed df to EPSG:5070
            df_hu = df_hu.to_crs(epsg=5070)

            # Compute centroids for each parcel (much faster than polygon intersections)
            parcel_centroids = regrow_shape_copy.copy()
            parcel_centroids["geometry"] = parcel_centroids.geometry.centroid

            # Spatial join: assign each parcel centroid to its watershed
            joined = gpd.sjoin(parcel_centroids, df_hu, how="left", predicate="within")

            # Select the watershed attribute columns you care about
            columns_to_keep = joined.filter(regex="field_id|_id|_name|_type").columns

            # Merge watershed attributes back to the original parcels
            regrow_shape_copy = regrow_shape_copy.merge(joined[list(columns_to_keep)], on="field_id", how="left", suffixes=("", "_temp"))

            # Clean up any temporary duplicate columns
            cols_to_drop = [col for col in regrow_shape_copy.columns if col.endswith("_temp")]
            regrow_shape_copy.drop(columns=cols_to_drop, inplace=True)
            print("Hydrography data added to attribute table.")

        except Exception as e:
            print(f"Error processing slope: {e}")
            raise
        
    # Save geojson and csv files
    #regrow_shape_copy.to_file(output_path_geojson, driver="GeoJSON")
    attribute_table = regrow_shape_copy.drop(columns='geometry')
    print(attribute_table.shape) #check df shape
    attribute_table.to_csv(output_path_csv, index=False)
    print(f'Supplementary data 3 for {state} is created and saved')

