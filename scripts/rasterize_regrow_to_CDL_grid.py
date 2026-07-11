import pandas as pd
import geopandas as gpd
import glob
import rasterio
from rasterio.features import rasterize
import os
from pathlib import Path

# Import parameters from Snakemake
regrow_input_folder = snakemake.params.regrow_input_dir
regrow_checks_folder = snakemake.params.regrow_checks_dir
CDL_input_folder = snakemake.params.CDL_input_dir
regrow_raster_output_folder = snakemake.params.regrow_raster_output_dir
states = snakemake.params.states
rasterization_year = snakemake.params.rasterization_year
target_CRS = snakemake.params.target_CRS


def preprocess_regrow(state, regrow_input_folder, regrow_checks_folder, regrow_raster_output_folder):
    
    # Path to input and output fields
    regrow_input_path = os.path.join(regrow_input_folder, f"{state}_regrow_fieldID_geometry.parquet")
    overlapping_fields_input_path = os.path.join(regrow_checks_folder, f"{state}_regrow_overlapping_fields.parquet")
    fieldID_pid_output_path = os.path.join(regrow_raster_output_folder, f"{state}_regrow_fieldID_pid_correspondence.parquet")
    
    # Load regrow shape (geometry) data
    regrow_geometry = gpd.read_parquet(regrow_input_path)
    
    ### Clean Regrow: Regrow geometries can overlap which harms rasterization, we need to remove overlaps ###
    regrow_overlaps = pd.read_parquet(overlapping_fields_input_path)
    # Extract all smaller parcel IDs from the pairs
    overlap_smaller_field_ids = regrow_overlaps["field_id_smaller"].unique()
    # Drop rows whose field_id is in smaller_ids and reset index
    regrow_geometry = (regrow_geometry[~regrow_geometry["field_id"].isin(overlap_smaller_field_ids)]).reset_index(drop=True)
  
    # In rasterio, if multiple polygons overlap a single pixel, the value from the last polygon in the input sequence is assigned to that pixel in the final raster
    # If polygons overlap, we want overlapping areas to take the value of the larger polygon
    # Sort polygons by area in ascending order so that larger polygons are processed last
    regrow_geometry = regrow_geometry.iloc[regrow_geometry.geometry.area.argsort()]
    
    
    # Mapping: Regrow field_id: unique field integers (pid)
    id_map = {s: i for i, s in enumerate(regrow_geometry["field_id"].unique())}
    regrow_geometry["pid"] = regrow_geometry["field_id"].map(id_map)
    # Save Regrow-pid mapping
    regrow_geometry[['field_id', 'pid']].to_parquet(fieldID_pid_output_path, compression="zstd")
    
    return regrow_geometry


def extract_raster_info(path):
    with rasterio.open(path) as src:
        return {
            "file": path,
            "crs": src.crs,
            "transform": src.transform,
            "res_x": src.res[0],
            "res_y": src.res[1],
            "width": src.width,
            "height": src.height,
            "bounds": src.bounds,
            "dtype": src.dtypes[0],
            "nodata": src.nodata
        }

def compare_to_template(row, template):
    return {
        "crs_match": row["crs"] == template["crs"],
        "res_match": (row["res_x"], row["res_y"]) == (template["res_x"], template["res_y"]),
        "size_match": (row["width"], row["height"]) == (template["width"], template["height"]),
        "transform_match": row["transform"] == template["transform"],
        "dtype_match": row["dtype"] == template["dtype"]
    }


def regrow_rasterization(regrow_geometry, state, rasterization_year, CDL_input_folder, regrow_raster_output_folder):
    
    # Path to input and output fields
    CDL_input_path = os.path.join(CDL_input_folder, f"{state}_{rasterization_year}_30m_cdls_clipped.tif")
    regrow_output_path = os.path.join(regrow_raster_output_folder, f"{state}_regrow_raster_to_CDL_grid.tif")

    ### Process CDL raster file to create mapping CDL pixel → Regrow pid ###
    with rasterio.open(CDL_input_path) as src:
        # Match CRS between vector Regrow and CDL raster
        regrow_geometry = regrow_geometry.to_crs(src.crs)

        # Properties of CDL raster file
        transform = src.transform
        height = src.height
        width = src.width

        # Rasterize fields: assign unique field integers (pid) to pixels whose centers lie inside a given polygon
        shapes = ((geom, pid) for geom, pid in zip(regrow_geometry.geometry, regrow_geometry.pid))
        
        parcel_raster = rasterize(
            shapes=shapes,
            out_shape=(height, width),
            transform=transform,
            fill=-1,            # pixels not belonging to any parcel
            all_touched=False,  # only pixels whose center is inside polygon
            dtype="int32",
            skip_invalid = False
        )
        
        # Save raster to disk
        out_meta = {
            "driver": "GTiff",
            "height": height,
            "width": width,
            "count": 1,
            "dtype": parcel_raster.dtype,
            "crs": src.crs,
            "transform": transform
        }

        with rasterio.open(regrow_output_path, "w", **out_meta) as dst:
            dst.write(parcel_raster, 1)
            
        print(f"Rasterization for {state} is complete and output files are saved.")


# Main code
# Make sure output folder exists
Path(regrow_raster_output_folder).mkdir(parents=True, exist_ok=True)

for state in states:
    # Adress the issue with overlapping fields
    regrow_geometry = preprocess_regrow(state, regrow_input_folder, regrow_checks_folder, regrow_raster_output_folder)
    
    # Find all raster files for the given state (e.g., Ohio_*.tif)
    raster_file_paths = glob.glob(f"{CDL_input_folder}/{state}_*.tif")

    # Extract key spatial parameters (CRS, resolution, transform, etc.) for each raster
    records = [extract_raster_info(f) for f in raster_file_paths]

    # Convert the list of dictionaries into a DataFrame (one row per raster)
    raster_parameter_df = pd.DataFrame(records)

    # Locate the template raster for the specified year (used as reference grid)
    template_file = glob.glob(f"{CDL_input_folder}/{state}_{rasterization_year}_*.tif")

    # Ensure exactly one template raster is found (avoid ambiguity)
    if len(template_file) != 1:
        raise ValueError(f"Expected 1 template raster, found {len(template_file)}: {template_file}")
    else:
        # Extract spatial parameters of the template raster
        template = extract_raster_info(template_file[0])

    # Compare each raster’s parameters to the template
    # Returns a DataFrame of boolean checks (e.g., crs_match, res_match, transform_match)
    qa_results = raster_parameter_df.apply(
        lambda r: compare_to_template(r, template),
        axis=1,
        result_type="expand"
    )

    # Combine file names with QA results for easier inspection
    qa_df = pd.concat([raster_parameter_df["file"], qa_results], axis=1)

    # Validate that ALL rasters match the template across ALL parameters
    # - First .all(): checks each column (parameter)
    # - Second .all(): checks overall result
    if not qa_df.drop(columns=["file"]).all().all():
        # Print which parameters failed across the dataset
        print(qa_df.drop(columns=["file"]).all())

        # Stop execution if any mismatch is found
        raise ValueError("Raster alignment check FAILED: at least one parameter does not match the template.")
    
    regrow_rasterization(regrow_geometry, state, rasterization_year, CDL_input_folder, regrow_raster_output_folder)