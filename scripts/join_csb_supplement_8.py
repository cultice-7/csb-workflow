import numpy as np
import pandas as pd
import geopandas as gpd
from collections import defaultdict
import rasterio
from rasterio.features import rasterize
import os
from pathlib import Path
import gc

# Import parameters from Snakemake
CSB_input_folder = snakemake.params.CSB_input_dir
rasterized_CSB_input_folder = snakemake.params.rasterized_CSB_input_dir
CSB_checks_folder = snakemake.params.CSB_checks_dir
CSB_output_folder = snakemake.params.CSB_output_dir
soil_input_folder =snakemake.params.soil_input_dir
soil_output_folder = snakemake.params.soil_output_dir
mukey_input_folder = snakemake.params.mukey_input_dir
states = snakemake.params.states
CSB_years = snakemake.params.CSB_years
soil_depth_cm = snakemake.params.soil_depth_cm # specify the maximum soil depth (in cm) to be used in variable calculations (weighting)


def process_gSSURGO_tabular_data(state, soil_depth_cm, soil_input_folder, soil_output_folder):
        
    ## Uploade tabular datasets and keep specific soil variables
    # Path to gSSURGO dataset for a specific state
    gdb_path = f"{soil_input_folder}/gSSURGO_{state}/gSSURGO_{state}.gdb"

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

    mapunit_component_chorizon_corerestictions_muaggatt_legend_sacatalog = mapunit_component_chorizon_corerestictions_muaggatt.merge(legend, on="lkey", how="left").merge(sacatalog, on="areasymbol", how="left")
    
    # Check duplicates
    if (mapunit_component_chorizon_corerestictions_muaggatt_legend_sacatalog.duplicated(subset=['mukey', 'cokey', 'hzdept_r'], keep=False)).sum() != 0:
        print(f"WARNING DUPLICATING ROWS: {mapunit_component_chorizon_corerestictions_muaggatt_legend_sacatalog[mapunit_component_chorizon_corerestictions_muaggatt_legend_sacatalog.duplicated(subset=['mukey', 'cokey', 'hzdept_r'], keep=False)]}") 
    
    # Save the integrated soil dataset
    #integrated_soil_output_path = os.path.join(soil_output_folder, f"{state}_integrated_soil_variables.parquet")
    #mapunit_component_chorizon_corerestictions_muaggatt_legend_sacatalog.to_parquet(integrated_soil_output_path, compression="zstd")
    
    
    ## Creating soil variables (the list of variables is provided by the soil team)
    top = 0
    bottom = soil_depth_cm
    
    # Set a mukey dataset for spatial join
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
    mukey_soil_variables = mukey_soil_variables.merge(dominant_component.filter(regex="mukey|dominant"), on="mukey", how="left")
    
    # Add variables that need to be weighted by soil layer composition: clay, sand, and pH
    mapunit_component_chorizon_copy = mapunit_component_chorizon.copy()
    
    # Compute horizon overlap with the set soil layer depth
    mapunit_component_chorizon_copy["overlap"] = np.maximum(
        0,
        np.minimum(mapunit_component_chorizon_copy["hzdepb_r"], bottom) - np.maximum(mapunit_component_chorizon_copy["hzdept_r"], top))

    # Function to compute horizon-weighted mean per component (depth-weighted mean for the first 30 cm using horizon thickness overlap)
    def horizon_weighted_target_depth(x, col):
        vals = pd.to_numeric(x[col], errors="coerce")
        weights = x["overlap"]

        mask = weights > 0
        vals = vals[mask]
        weights = weights[mask]

        if weights.sum() == 0:
            return np.nan

        return np.sum(vals * weights) / np.sum(weights)
    
    # Compute component-level values
    component_weighted_target_depth = (
        mapunit_component_chorizon_copy.groupby(["mukey", "cokey", "comppct_r"])
        .apply(lambda x: pd.Series({
            f"claytotal_r_{soil_depth_cm}cm": horizon_weighted_target_depth(x, "claytotal_r"),
            f"sandtotal_r_{soil_depth_cm}cm": horizon_weighted_target_depth(x, "sandtotal_r"),
            f"ph1to1h2o_r_{soil_depth_cm}cm": horizon_weighted_target_depth(x, "ph1to1h2o_r"),
            f"om_r_{soil_depth_cm}cm": horizon_weighted_target_depth(x, "om_r")
        }), include_groups=False)
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
        component_weighted_target_depth.groupby("mukey")
        .apply(lambda x: pd.Series({
            f"claytotal_r_{soil_depth_cm}cm_weighted": component_weighted_mean(x, f"claytotal_r_{soil_depth_cm}cm"),
            f"sandtotal_r_{soil_depth_cm}cm_weighted": component_weighted_mean(x, f"sandtotal_r_{soil_depth_cm}cm"),
            f"ph1to1h2o_r_{soil_depth_cm}cm_weighted": component_weighted_mean(x, f"ph1to1h2o_r_{soil_depth_cm}cm"),
            f"om_r_{soil_depth_cm}cm_weighted": component_weighted_mean(x, f"om_r_{soil_depth_cm}cm")
        }), include_groups=False)
        .reset_index()
    )

    # Add weighted variables variables to the integrated dataset
    mukey_soil_variables = mukey_soil_variables.merge(mukey_weighted, on="mukey", how="left")

    # Add slope data to the integrated dataset
    mukey_soil_variables = mukey_soil_variables.merge(muaggatt, on="mukey", how="left")

    return mukey_soil_variables



def merge_CSB_mukey_soilvars(state, year, mukey_soil_variables, CSB_input_folder, rasterized_CSB_input_folder, CSB_checks_folder, mukey_input_folder, CSB_output_folder):
            
    CSB_table_input_path = os.path.join(CSB_input_folder, f"{state}_CSB{year}_table.parquet")
    CSB_raster_input_path = os.path.join(rasterized_CSB_input_folder, f"{state}_CSB{year}_raster_to_gSSURGO_grid.tif")
    gssurgo_mukey_input_path = os.path.join(mukey_input_folder, "gSSURGO Mukey Grid", f"{state}_MURASTER_30m.tif")
    overlapping_fields_input_path = os.path.join(CSB_checks_folder, f"{state}_CSB{year}_overlapping_fields.parquet")
    CSBID_pid_input_path = os.path.join(rasterized_CSB_input_folder, f"{state}_CSB{year}_CSBID_pid_correspondence.parquet")
    output_path_table = os.path.join(CSB_output_folder, f"{state}_CSB{year}_supplement_8_table.parquet")
    
    
    # Load CSB_dises joined datasets
    CSB_table = pd.read_parquet(CSB_table_input_path)
    CSB_table = CSB_table[['CSBID']]
    
    # Create soil dataset
    CSB_soil = CSB_table.copy()

    
    
    ### Process mukey raster file to create pairs of CSBID CSBID: mukeys (pixel values) ###
    # Upload the dataset with a mapping from CSB CSBID to unique field integers
    id_map = pd.read_parquet(CSBID_pid_input_path)
    id_map = id_map.set_index("pid")
    
    # Read both CSB and mukey rasters at once
    with rasterio.open(gssurgo_mukey_input_path) as src_mukey, \
        rasterio.open(CSB_raster_input_path) as src_CSB:
        mukey_raster  = src_mukey.read(1)
        CSB_raster = src_CSB.read(1)

    # Extract pixel values for each field
    # Indices of all pixels that belong to some field
    mask = CSB_raster >= 0
    CSB_pids = CSB_raster[mask]
    mukeys = mukey_raster[mask].astype(str)
    
    # Remove raster files from the memory
    del mukey_raster, CSB_raster
    gc.collect()

    # Build correnspondence between CSBID and pixel values
    CSBID_mukeys = defaultdict(list)
    CSB_ids = id_map.loc[CSB_pids, "CSBID"].values

    for CSB_id, mukey in zip(CSB_ids, mukeys):
        CSBID_mukeys[CSB_id].append(mukey)

    del CSB_pids, mukeys
    gc.collect()
    
    

    ### Compute soil variables for each CSBID field ###
    # Split numerical and categorical columns
    num_cols = mukey_soil_variables.select_dtypes(include='number').columns.drop(['mukey', 'lkey'], errors='ignore')
    cat_cols = mukey_soil_variables.select_dtypes(exclude='number').columns.drop(['mukey', 'lkey'], errors='ignore')

    # Index soil dataset to speed up searching
    soil_indexed = mukey_soil_variables.set_index("mukey")

    # Compute soil variables for each CSBID field by taking average (mode) across all mukeys within a given field
    rows = []
    for CSBID, mukeys in CSBID_mukeys.items():
        if not mukeys:
            continue
        
        sub = soil_indexed.reindex(mukeys).dropna(how="all")
        if sub.empty:
            continue

        num_mean = sub[num_cols].mean()
        cat_mode = sub[cat_cols].apply(lambda x: x.dropna().value_counts().idxmax() if not x.dropna().empty else None)

        row = pd.concat([num_mean, cat_mode])
        row["CSBID"] = CSBID

        rows.append(row)
    
    soil_data_for_CSB = pd.DataFrame(rows)
    CSB_soil = CSB_soil.merge(soil_data_for_CSB, on='CSBID', how="left")

    # Map each smaller_CSBID to its corresponding larger_CSBID and then copy the values from the larger fields back into the target dataset
    CSB_overlaps = pd.read_parquet(overlapping_fields_input_path)
    CSB_soil = CSB_soil.set_index("CSBID")
    CSB_soil.loc[CSB_overlaps["CSBID_smaller"]] = CSB_soil.loc[CSB_overlaps["CSBID_larger"]].values
    CSB_soil = CSB_soil.reset_index()
    
    ### Save output files ###
    # Convert all float64 to float32
    float64_cols = CSB_soil.select_dtypes(include="float64").columns
    CSB_soil[float64_cols] = CSB_soil[float64_cols].astype("float32")

    CSB_soil.to_parquet(output_path_table, compression="zstd")
    print(f"Creating and saving soil dataset for {state} and {year} is complete")
        

# Main code
for year in CSB_years:
    for state in states:

        mukey_soil_variables = process_gSSURGO_tabular_data(state, soil_depth_cm, soil_input_folder, soil_output_folder)
        
        merge_CSB_mukey_soilvars(state, year, mukey_soil_variables, CSB_input_folder, rasterized_CSB_input_folder, CSB_checks_folder, mukey_input_folder, CSB_output_folder)