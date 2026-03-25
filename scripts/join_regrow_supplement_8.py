import numpy as np
import pandas as pd
import geopandas as gpd
import rasterio
from rasterio.features import rasterize
import os
from pathlib import Path
import pickle

# Input and output folders for Regrow
input_folder_Regrow = snakemake.params.regrow_input_dir
output_folder_Regrow = snakemake.params.regrow_output_dir
output_folder_soil = snakemake.params.soil_output_dir

# Pull list of states for running the code
states = snakemake.params.states
# Specify the maximum soil depth (in cm) to be used in variable calculations (weighting)
soil_depth_cm = snakemake.params.soil_depth_cm

for state in states:
    
    input_path_Regrow_geometry = os.path.join(input_folder_Regrow, f"{state}_regrow_fieldID_geometry.parquet")
    input_path_Regrow_raster = os.path.join(input_folder_Regrow, f"{state}_regrow_raster_to_gSSURGO_grid.tif")
    input_path_duplicating_fields = os.path.join(input_folder_Regrow, f"{state}_regrow_duplicating_fields.parquet")
    output_path_spatial = os.path.join(output_folder_Regrow, f"{state}_regrow_supplement_8_spatial.parquet")
    output_path_table = os.path.join(output_folder_Regrow, f"{state}_regrow_supplement_8_table.parquet")
    output_path_integrated_soil = os.path.join(output_folder_soil, f"{state}_integrated_soil_variables.parquet")
    
    # Load regrow_dises joined datasets
    regrow_geometry = gpd.read_parquet(input_path_Regrow_geometry)
    
    # Reproject geometry to the same CRS (NAD83/CONUS Albers)
    regrow_geometry = regrow_geometry.to_crs(epsg=5070)
    
    # Create soil dataset
    regrow_soil = regrow_geometry[['field_id']]
    
    
    
    ### Clean Regrow: Regrow geometries contain overlaps which harms rasterization, we need to remove overlaps ###
    regrow_overlaps = pd.read_parquet(input_path_duplicating_fields)
    
    # Clean regrow fields to remove overlapping fields
    # Extract all smaller parcel IDs from the pairs
    overlap_smaller_field_ids = regrow_overlaps["field_id_smaller"].unique()
    # Drop rows whose field_id is in smaller_ids and reset index
    regrow_geometry = (regrow_geometry[~regrow_geometry["field_id"].isin(overlap_smaller_field_ids)]).reset_index(drop=True)
        
    
    
    ### Process mukey raster file to create pairs of Regrow field_id: mukeys (pixel values) ###
    # Regrow field_id → unique field integers
    id_map = {s: i for i, s in enumerate(regrow_geometry["field_id"].unique())}
    # Unique field integers → Regrow field_id
    reverse_id_map = {v: k for k, v in id_map.items()} 
    regrow_geometry["pid"] = regrow_geometry["field_id"].map(id_map)

    # Open gSSURGO mukey raster
    with rasterio.open(f"data/edited/Soil/gSSURGO Mukey Grid/{state}_MURASTER_30m.tif") as src_gSSURGO:
        # Read gSSURGO mukey raster
        mukey_raster = src_gSSURGO.read(1)

    # Open rasterized Regrow
    with rasterio.open(input_path_Regrow_raster) as src_regrow:
        # Read rasterized Regrow
        regrow_raster = src_regrow.read(1)

    # Extract pixel values for each field
    # Indices of all pixels that belong to some field
    rows, cols = np.where(regrow_raster >= 0)
    parcel_ids = regrow_raster[rows, cols]
    values = mukey_raster[rows, cols].astype(str)

    # Build dictionary of field_id → pixel values
    fieldID_mukeys = {field_id: [] for field_id in regrow_geometry.field_id}
    for pid, val in zip(parcel_ids, values):
        fieldID_mukeys[reverse_id_map[pid]].append(val)
        
    # Find and print fields with no pixels assigned
    print(f"Regrow fields in {state} with no pixels assigned")
    for k, v in fieldID_mukeys.items():
        if len(v) == 0:
            print(k, v)
    
    
    
    ### Process gSSURGO tabular data ###
    
    ## Uploade tabular datasets and keep specific soil variables
    # Path to gSSURGO dataset for a specific state
    gdb_path = f"data/gSSURGO/gSSURGO_{state}/gSSURGO_{state}.gdb"

    # Load a layer "mapunit" and keep specific columns from the list of soil variables
    mapunit = gpd.read_file(gdb_path, layer="mapunit")
    mapunit = mapunit[['mukey', 'lkey']]
    # Drop any duplicates from the shrinked mapunit dataset
    mapunit.drop_duplicates(inplace=True)

    # Load a layer "component" and keep specific columns from the list of soil variables
    component = gpd.read_file(gdb_path, layer="component")
    component = component[['mukey', 'cokey', 'majcompflag', 'compname', 'comppct_r', 'drainagecl', 'cropprodindex']]
    # Drop any duplicates from the shrinked component dataset
    component.drop_duplicates(inplace=True)

    # Load a layer "chorizon" and keep specific columns from the list of soil variables
    chorizon = gpd.read_file(gdb_path, layer="chorizon")
    chorizon = chorizon[['cokey', 'hzdept_r', 'hzdepb_r', 'sandtotal_r', 'claytotal_r', 'ph1to1h2o_r', 'om_r']]
    # Drop any duplicates from the shrinked chorizon dataset
    chorizon.drop_duplicates(inplace=True)

    # Load a layer "component" and keep specific columns from the list of soil variables
    corestrictions = gpd.read_file(gdb_path, layer="corestrictions")
    corestrictions = corestrictions[['cokey', 'resdept_r']]
    # Drop any duplicates from the shrinked corestrictions dataset
    corestrictions.drop_duplicates(inplace=True)

    # Load a layer "muaggatt" and keep specific columns from the list of soil variables
    muaggatt = gpd.read_file(gdb_path, layer="muaggatt")
    muaggatt = muaggatt[['mukey', 'slopegraddcp']]
    # Drop any duplicates from the shrinked muaggatt dataset
    muaggatt.drop_duplicates(inplace=True)

    # Load a layer "legend" and keep specific columns
    legend = gpd.read_file(gdb_path, layer="legend")
    legend = legend[["lkey", "areasymbol"]]
    # Drop any duplicates from the shrinked legend dataset
    legend.drop_duplicates(inplace=True)

    # Load a layer "sacatalog" and keep specific columns
    sacatalog = gpd.read_file(gdb_path, layer="sacatalog")
    sacatalog = sacatalog[["areasymbol", "saverest"]]
    # Delete time from "saverest" datetime column
    sacatalog["saverest"] = pd.to_datetime(sacatalog["saverest"]).dt.date
    # Drop any duplicates from the shrinked sacatalog dataset
    sacatalog.drop_duplicates(inplace=True)
    
    ## Dataset cleaning
    # Check duplicates in all datasets
    if (mapunit.duplicated(keep=False)).sum() != 0:
        print(f"WARNING DUPLICATING ROWS: {mapunit[mapunit.duplicated(keep=False)]}")

    if (component.duplicated(subset=['mukey', 'cokey'], keep=False)).sum() != 0:
        print(f"WARNING DUPLICATING ROWS: {component[component.duplicated(subset=['mukey', 'cokey'], keep=False)]}")

    if (chorizon.duplicated(subset=['cokey', 'hzdept_r'], keep=False)).sum() != 0:
        print(f"WARNING DUPLICATING ROWS: {chorizon[chorizon.duplicated(subset=['cokey', 'hzdept_r'], keep=False)]}")

    if (corestrictions.duplicated(subset=['cokey'], keep=False)).sum() != 0:
        print(f"WARNING DUPLICATING ROWS: {corestrictions[corestrictions.duplicated(subset=['cokey'], keep=False)]}")

    if (muaggatt.duplicated(subset=['mukey'], keep=False)).sum() != 0:
        print(f"WARNING DUPLICATING ROWS: {muaggatt[muaggatt.duplicated(subset=['mukey'], keep=False)]}")

    if (legend.duplicated(subset=['lkey'], keep=False)).sum() != 0:
        print(f"WARNING DUPLICATING ROWS: {legend[legend.duplicated(subset=['lkey'], keep=False)]}")

    if (sacatalog.duplicated(subset=['areasymbol'], keep=False)).sum() != 0:
        print(f"WARNING DUPLICATING ROWS: {sacatalog[sacatalog.duplicated(subset=['areasymbol'], keep=False)]}")
    
    # Keep only minimum soil restrictive layer for each cokey
    corestrictions = corestrictions.groupby("cokey")["resdept_r"].min().reset_index()
    
    ## Dataset merging
    # Merge mapunit and component on the key "mukey" first. Then, merge corestrictions, legend and sacatalog
    mapunit_component = mapunit.merge(component, on="mukey", how="left")
    mapunit_component_corestr_legend_sacatalog = mapunit_component.merge(corestrictions, on="cokey", how="left").merge(legend, on="lkey", how="left").merge(sacatalog, on="areasymbol", how="left")
    
    # Check duplicates
    if (mapunit_component_corestr_legend_sacatalog.duplicated(subset=['mukey', 'cokey'], keep=False)).sum() != 0:
        print(f"WARNING DUPLICATING ROWS: {mapunit_component_corestr_legend_sacatalog[mapunit_component_corestr_legend_sacatalog.duplicated(subset=['mukey', 'cokey'], keep=False)]}")
    
    # Create integrated soil dataset by merging all datasets with each other in the following order: 
    # 1) component, chorizon and corerestictions on the key "cokey"
    # 2) the mapunit, merged dataset from (1), and muaggatt on the key "mukey"
    # 3) the merged dataset from (2), legend and sacatalog on the keys "lkey" and "areasymbol"
    component_chorizon = component.merge(chorizon, on="cokey", how="left")
    component_chorizon_corerestictions = component_chorizon.merge(corestrictions, on="cokey", how="left")

    mapunit_component_chorizon = mapunit.merge(component_chorizon, on="mukey", how="left")
    mapunit_component_chorizon_corerestictions = mapunit.merge(component_chorizon_corerestictions, on="mukey", how="left")
    mapunit_component_chorizon_corerestictions_muaggatt = mapunit_component_chorizon_corerestictions.merge(muaggatt, on="mukey", how="left")

    mapunit_component_chorizon_corerestictions_muaggatt_legend = mapunit_component_chorizon_corerestictions_muaggatt.merge(legend, on="lkey", how="left")
    mapunit_component_chorizon_corerestictions_muaggatt_legend_sacatalog = mapunit_component_chorizon_corerestictions_muaggatt_legend.merge(sacatalog, on="areasymbol", how="left")
    
    # Check duplicates
    if (mapunit_component_chorizon_corerestictions_muaggatt_legend_sacatalog.duplicated(subset=['mukey', 'cokey', 'hzdept_r'], keep=False)).sum() != 0:
        print(f"WARNING DUPLICATING ROWS: {mapunit_component_chorizon_corerestictions_muaggatt_legend_sacatalog[mapunit_component_chorizon_corerestictions_muaggatt_legend_sacatalog.duplicated(subset=['mukey', 'cokey', 'hzdept_r'], keep=False)]}") 
    
    # Save the integrated soil dataset
    mapunit_component_chorizon_corerestictions_muaggatt_legend_sacatalog.to_parquet(output_path_integrated_soil, compression="zstd")
    
    
    ## Creating soil variables (the list of variables is provided by the soil team)
    top = 0
    bottom = soil_depth_cm
    
    # Set an integrated mukey dataset for spatial join
    mukey_soil_variables = mapunit.copy()
    
    # Variables which values are taken only from the dominant component
    dominant_component = (
        mapunit_component_corestr_legend_sacatalog[mapunit_component_corestr_legend_sacatalog["majcompflag"] == "Yes"]
        .sort_values(["mukey", "comppct_r"], ascending=[True, False])
        .groupby("mukey", as_index=False)
        .first()
    )
    
    # Rename variables
    dominant_component.rename(columns={"compname": "compname_dominant",
                                    "comppct_r": "comppct_r_dominant",
                                    "drainagecl": "drainagecl_dominant",
                                    "cropprodindex": "cropprodindex_dominant",
                                    "resdept_r":"resdept_r_dominant",
                                    "areasymbol":"areasymbol_dominant",
                                    "saverest":"saverest_dominant"}, inplace=True)

    # Add dominant component variables to the integrated dataset
    mukey_soil_variables = mukey_soil_variables.merge(dominant_component.filter(regex="dominant|mukey"), on="mukey", how="left")
    
    # Add variables that need to be weighted by soil layer composition: clay, sand, and pH
    mapunit_component_chorizon_copy = mapunit_component_chorizon.copy()
    
    # Compute horizon overlap with 0–30 cm
    mapunit_component_chorizon_copy["overlap"] = np.maximum(
        0,
        np.minimum(mapunit_component_chorizon_copy["hzdepb_r"], bottom) - np.maximum(mapunit_component_chorizon_copy["hzdept_r"], top))

    # Function to compute horizon-weighted mean per component (depth-weighted mean for the first 30 cm using horizon thickness overlap)
    def horizon_weighted_30(x, col):
        vals = pd.to_numeric(x[col], errors="coerce")
        weights = x["overlap"]

        mask = weights > 0
        vals = vals[mask]
        weights = weights[mask]

        if weights.sum() == 0:
            return np.nan

        return np.sum(vals * weights) / np.sum(weights)
    
    # Compute component-level values
    component_weighted_30 = (
        mapunit_component_chorizon_copy.groupby(["mukey", "cokey", "comppct_r"])
        .apply(lambda x: pd.Series({
            "claytotal_r_30cm": horizon_weighted_30(x, "claytotal_r"),
            "sandtotal_r_30cm": horizon_weighted_30(x, "sandtotal_r"),
            "ph1to1h2o_r_30cm": horizon_weighted_30(x, "ph1to1h2o_r"),
            "om_r_30cm": horizon_weighted_30(x, "om_r")
        }))
        .reset_index())
    
    # Component-weighted mean for each mukey (component-weighted mean using comppct_r)
    def component_weighted_mean(x, col):
        vals = x[col]
        weights = x["comppct_r"]

        mask = np.isfinite(vals) & np.isfinite(weights)
        vals = vals[mask]
        weights = weights[mask]

        if weights.sum() == 0:
            return np.nan

        return np.sum(vals * weights) / np.sum(weights)
    
    # Final mukey-level soil values
    mukey_weighted = (
        component_weighted_30.groupby("mukey")
        .apply(lambda x: pd.Series({
            "claytotal_r_30cm_weighted": component_weighted_mean(x, "claytotal_r_30cm"),
            "sandtotal_r_30cm_weighted": component_weighted_mean(x, "sandtotal_r_30cm"),
            "ph1to1h2o_r_30cm_weighted": component_weighted_mean(x, "ph1to1h2o_r_30cm"),
            "om_r_30cm_weighted": component_weighted_mean(x, "om_r_30cm")
        }))
        .reset_index()
    )

    # Add weighted variables variables to the integrated dataset
    mukey_soil_variables = mukey_soil_variables.merge(mukey_weighted, on="mukey", how="left")

    # Add slope data to the integrated dataset
    mukey_soil_variables = mukey_soil_variables.merge(muaggatt, on="mukey", how="left")



    ### Compute soil variables for each Regrow field ###
    # Split numerical and categorical columns
    num_cols = mukey_soil_variables.select_dtypes(include='number').columns.drop(['mukey', 'lkey'], errors='ignore')
    cat_cols = mukey_soil_variables.select_dtypes(exclude='number').columns.drop(['mukey', 'lkey'], errors='ignore')

    # Index soil and target datasets to speed up searching
    soil_indexed = mukey_soil_variables.set_index("mukey")
    regrow_soil = regrow_soil.set_index("field_id")

    # Compute soil variables for each Regrow field by taking average (mode) across all mukeys within a given field
    for fieldID, mukeys in fieldID_mukeys.items():
        if mukeys:
            sub = soil_indexed.reindex(mukeys).dropna(how="all")
            if sub.empty:
                continue

            num_mean = sub[num_cols].mean()
            cat_mode = sub[cat_cols].apply(lambda x: x.dropna().value_counts().idxmax() if not x.dropna().empty else None)

            row = pd.concat([num_mean, cat_mode])

            regrow_soil.loc[fieldID, row.index] = row.values

    # Map each smaller_field_id to its corresponding larger_field_id and then copy the values from the larger fields back into the target dataset
    regrow_soil.loc[regrow_overlaps["field_id_smaller"]] = regrow_soil.loc[regrow_overlaps["field_id_larger"]].values
    regrow_soil = regrow_soil.reset_index()
    
    ### Save output files ###
    # Convert all float64 to float32
    float64_cols = regrow_soil.select_dtypes(include="float64").columns
    regrow_soil[float64_cols] = regrow_soil[float64_cols].astype("float32")
    
    regrow_soil.to_parquet(output_path_table, compression="zstd")
    print(f"Creating and saving soil dataset for {state} is complete")