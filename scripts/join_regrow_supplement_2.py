import geopandas as gpd
from rasterstats import zonal_stats
import os
from pathlib import Path

# Import parameters from Snakemake
regrow_input_folder = snakemake.params.regrow_input_dir
regrow_output_folder = snakemake.params.regrow_output_dir
elevation_input_folder = snakemake.params.elevation_input_dir
slope_input_folder = snakemake.params.slope_input_dir
states = snakemake.params.states


for state in states:

    # Paths to the reprojected slope and elevation raster files
    elevation_proj_path = os.path.join(elevation_input_folder, f"{state}_elevation_clipped.tif")
    slope_proj_path = os.path.join(slope_input_folder, f"{state}_slope_clipped.tif")

    input_path_Regrow = os.path.join(regrow_input_folder, f"{state}_regrow_fieldID_geometry.parquet")
    
    # Load Regrow data
    regrow_shape = gpd.read_parquet(input_path_Regrow)

    # Setting active geometry column
    regrow_shape = regrow_shape.set_geometry('geometry')
    
    # Keep only the columns necessary for the spatial join
    cols_to_keep = ['field_id', 'geometry']
    regrow_shape = regrow_shape[cols_to_keep]
     
    
    #---# Add zonal statistics for elevation and slope
    # Set paths to output files
    output_path_spatial = os.path.join(regrow_output_folder, f"{state}_regrow_supplement_2_spatial.parquet")
    output_path_table = os.path.join(regrow_output_folder, f"{state}_regrow_supplement_2_table.parquet")
    
    # Zonal statistics for elevation
    try:
        print("Calculating and adding mean elevation...")
        elevation_stats = zonal_stats(regrow_shape, elevation_proj_path, stats="mean")
        elevation_means = [stat['mean'] for stat in elevation_stats]
        regrow_shape['elevation_mean'] = elevation_means
        print("Mean elevation added to attribute table.")
    except Exception as e:
        print(f"Error processing elevation: {e}")
        raise

    # Zonal statistics for slope
    try:
        print("Calculating and adding mean slope...")
        slope_stats = zonal_stats(regrow_shape, slope_proj_path, stats="mean")
        slope_means = [stat['mean'] for stat in slope_stats]
        regrow_shape['slope_mean'] = slope_means
        print("Mean slope added to attribute table.")
    except Exception as e:
        print(f"Error processing slope: {e}")
        raise
    
    # Check whether there are any missing values
    cols_to_check = ['field_id', 'elevation_mean', 'slope_mean']
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
    print(f'Supplementary data 2 for {state} is created and saved')
    

