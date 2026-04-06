import geopandas as gpd
from rasterstats import zonal_stats
import os

# Import parameters from Snakemake
CSB_input_folder = snakemake.params.CSB_input_dir
CSB_output_folder = snakemake.params.CSB_output_dir
elevation_input_folder = snakemake.params.elevation_input_dir
slope_input_folder = snakemake.params.slope_input_dir
states = snakemake.params.states
CSB_years = snakemake.params.CSB_years


for year in CSB_years:
    for state in states:
        
        # Paths to the reprojected slope and elevation raster files
        elevation_proj_path = os.path.join(elevation_input_folder, f"{state}_elevation_clipped.tif")
        slope_proj_path = os.path.join(slope_input_folder, f"{state}_slope_clipped.tif")
        
        input_path_CSB = os.path.join(CSB_input_folder, f"{state}_CSB{year}_CSBID_geometry.parquet")

        # Load CSB data
        CSB_shape = gpd.read_parquet(input_path_CSB)

        # Setting active geometry column
        CSB_shape = CSB_shape.set_geometry('geometry')
        
        # Keep only the columns necessary for the spatial join
        cols_to_keep = ['CSBID', 'geometry']
        CSB_shape = CSB_shape[cols_to_keep]
        
        
        #---# Add zonal statistics for elevation and slope
        # Set paths to output files
        output_path_spatial = os.path.join(CSB_output_folder, f"{state}_CSB{year}_supplement_2_spatial.parquet")
        output_path_table = os.path.join(CSB_output_folder, f"{state}_CSB{year}_supplement_2_table.parquet")
        
        # Zonal statistics for elevation
        try:
            print("Calculating mean elevation...")
            elevation_stats = zonal_stats(CSB_shape, elevation_proj_path, stats="mean")
            elevation_means = [stat['mean'] for stat in elevation_stats]
            CSB_shape['elevation_mean'] = elevation_means
            print("Mean elevation added to attribute table.")
        except Exception as e:
            print(f"Error processing elevation: {e}")
            raise

        # Zonal statistics for slope
        try:
            print("Calculating mean slope...")
            slope_stats = zonal_stats(CSB_shape, slope_proj_path, stats="mean")
            slope_means = [stat['mean'] for stat in slope_stats]
            CSB_shape['slope_mean'] = slope_means
            print("Mean slope added to attribute table.")
        except Exception as e:
            print(f"Error processing slope: {e}")
        
        # Check whether there are any missing values
        cols_to_check = ['CSBID', 'elevation_mean', 'slope_mean']
        if CSB_shape[cols_to_check].isna().any().any():
            print(CSB_shape[CSB_shape[cols_to_check].isna().any(axis=1)])
            
        # Convert float64 columns to float32 to save memory
        float64_cols = CSB_shape.select_dtypes(include=["float64"]).columns
        CSB_shape[float64_cols] = CSB_shape[float64_cols].astype("float32")
        
        #---# Save files w/ and w/o geometry
        #CSB_shape.to_parquet(output_path_spatial, compression="zstd")
        attribute_table = CSB_shape.drop(columns='geometry')
        print(attribute_table.shape) #check df shape
        attribute_table.to_parquet(output_path_table, index=False, compression="zstd")
        print(f'Supplementary data 2 for {state} is created and saved')
        