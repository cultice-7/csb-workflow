import geopandas as gpd
import os

#---# Load required datasets

# Watershed data
subbasin = gpd.read_file("data/Geo/watershed/subbasin.shp")
watershed = gpd.read_file("data/Geo/watershed/watershed.shp")
subwatershed = gpd.read_file("data/Geo/watershed/subwatershed.shp")

# Input and output folders for Regrow
input_folder_Regrow = "data/edited/Regrow/"
output_folder_Regrow = "data/edited/Regrow/"

# List of states
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


    #---# Add subbasin, watershed and subwatershed
    # Set paths to output files
    output_path_geojson = os.path.join(output_folder_Regrow, f"{state}_regrow_supplement_3_spatial.geojson")
    output_path_csv = os.path.join(output_folder_Regrow, f"{state}_regrow_supplement_3_table.csv")
    
    print("Adding hydrography data...")
    # Each hydrological unit has a specific dataset structure, so we process them separately
    # Codes of hydrological units: 8 - subbasin, 10 - watershed, 12 - subwatershed
    # Each dataset transformation steps: remove duplicating rows -> keep only needed columns -> rename columns
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
            #---# Polygon–polygon intersections (too time-consuming)
            #df_hu = df_hu.to_crs(epsg=5070) 
            #intersections_hydro = gpd.overlay(regrow_shape, df_hu, how='intersection')
            #intersections_hydro['overlap_area_temp'] = intersections_hydro.geometry.area
            #largest_overlap_hydro = intersections_hydro.sort_values('overlap_area_temp', ascending=False).drop_duplicates('field_id')
            #columns_to_merge_hydro = largest_overlap_hydro.columns[largest_overlap_hydro.columns.str.contains('field_id|_id|_name|_type', regex=True)]
            #regrow_shape = regrow_shape.merge(largest_overlap_hydro[columns_to_merge_hydro], on='field_id', how='left', suffixes=('', '_temp'))
            #cols_to_drop = [col for col in regrow_shape.columns if col.endswith('_temp')]
            #regrow_shape.drop(columns=cols_to_drop, inplace=True)
            #print("Hydrography data added to attribute table.")
            
            #---# Centroid point-in-polygon spatial joins (faster tool)
            # Reproject watershed df to EPSG:5070
            df_hu = df_hu.to_crs(epsg=5070)

            # Compute centroids for each parcel (much faster than polygon intersections)
            parcel_centroids = regrow_shape.copy()
            parcel_centroids["geometry"] = parcel_centroids.geometry.centroid

            # Spatial join: assign each parcel centroid to its watershed
            joined = gpd.sjoin(parcel_centroids, df_hu, how="left", predicate="within")

            # Select the watershed attribute columns you care about
            columns_to_keep = joined.filter(regex="field_id|_id|_name|_type").columns

            # Merge watershed attributes back to the original parcels
            regrow_shape = regrow_shape.merge(joined[list(columns_to_keep)], on="field_id", how="left", suffixes=("", "_temp"))

            # Clean up any temporary duplicate columns
            cols_to_drop = [col for col in regrow_shape.columns if col.endswith("_temp")]
            regrow_shape.drop(columns=cols_to_drop, inplace=True)
            print("Hydrography data added to attribute table.")

        except Exception as e:
            print(f"Error processing slope: {e}")
            raise
    
    # Check whether there are any missing values
    cols_to_check = ['field_id', 'subbasin_id', 'watershed_id', 'subwatershed_id']
    if regrow_shape[cols_to_check].isna().any().any():
        print(regrow_shape[regrow_shape[cols_to_check].isna().any(axis=1)])

    # Save geojson and csv files
    #regrow_shape.to_file(output_path_geojson, driver="GeoJSON")
    attribute_table = regrow_shape.drop(columns='geometry')
    print(attribute_table.shape) #check df shape
    attribute_table.to_csv(output_path_csv, index=False)
    print(f'Supplementary data 3 for {state} is created and saved')