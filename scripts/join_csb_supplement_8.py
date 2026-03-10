import numpy as np
import pandas as pd
import geopandas as gpd
import rasterio
from rasterio.features import rasterize
import os
from pathlib import Path


# Input and output folders for CSB
input_folder_CSB = "data/edited/CSB/"
output_folder_CSB = "data/edited/CSB/"
output_folder_soil = "data/edited/Soil/"

# Pull list of states for running the code
states = snakemake.params.states
# Specify the maximum soil depth (in cm) to be used in variable calculations (weighting)
soil_depth_cm = snakemake.params.soil_depth_cm

for state in states:
    
    input_path_CSB = os.path.join(input_folder_CSB, f"{state}_CSB1724_CSBID_geometry.parquet")
    output_path_spatial = os.path.join(output_folder_CSB, f"{state}_CSB1724_supplement_8_spatial.parquet")
    output_path_table = os.path.join(output_folder_CSB, f"{state}_CSB1724_supplement_8_table.parquet")
    output_path_integrated_soil = os.path.join(output_folder_soil, f"{state}_integrated_soil_variables.parquet")
    
    # Load CSB_dises joined datasets
    CSB_geometry = gpd.read_parquet(input_path_CSB)
    
    # Reproject geometry to the same CRS (NAD83/CONUS Albers)
    CSB_geometry = CSB_geometry.to_crs(epsg=5070)
    
    # Create soil dataset
    CSB_soil = CSB_geometry[['CSBID']]
    
    
    
    ### Clean CSB: CSB geometries contain overlaps which harms rasterization, we need to remove overlaps ###
    gdf = CSB_geometry.copy()
    gdf["geometry"] = gdf["geometry"].buffer(0)
    gdf["area"] = gdf.geometry.area

    # Overlay CSBID on itself. This will produce all overlapping polygons (including boundaries)
    CSB_overlaps = gpd.overlay(gdf, gdf, how="intersection")

    # Remove self-overlaps
    CSB_overlaps = CSB_overlaps[CSB_overlaps["CSBID_1"] != CSB_overlaps["CSBID_2"]]

    # Compute overlap ratio
    # Determine the smaller area in each pair
    CSB_overlaps["min_area"] = CSB_overlaps[["area_1", "area_2"]].min(axis=1)

    # Compute overlap fraction relative to smaller polygon
    CSB_overlaps["overlap_fraction"] = CSB_overlaps.geometry.area / CSB_overlaps["min_area"]

    # Keep only those with overlap >= 50%
    CSB_overlaps = CSB_overlaps[CSB_overlaps["overlap_fraction"] >= 0.5]

    # Keep only unique pairs (A,B) same as (B,A)
    # Make a tuple of (smaller_area_parcel, larger_area_parcel)
    CSB_overlaps["pair"] = CSB_overlaps.apply(
        lambda row: (row["CSBID_1"], row["CSBID_2"]) if row["area_1"] <= row["area_2"] else (row["CSBID_2"], row["CSBID_1"]), axis=1)
    # Drop duplicates
    CSB_overlaps = CSB_overlaps.drop_duplicates(subset="pair").reset_index(drop=True)
    # Assign smaller field IDs in each pair to CSBID_1 column and larger field IDs to CSBID_2
    CSB_overlaps.loc[:, 'CSBID_1'] = CSB_overlaps["pair"].apply(lambda x: x[0] if pd.notna(x) else None)
    CSB_overlaps.loc[:, 'CSBID_2'] = CSB_overlaps["pair"].apply(lambda x: x[1] if pd.notna(x) else None)
    # Rename CSBID_1 and CSBID_2 accordingly
    CSB_overlaps.rename(columns={"CSBID_1": "CSBID_smaller", "CSBID_2": "CSBID_larger"}, inplace = True)
    # Keep only necessary columns
    CSB_overlaps = CSB_overlaps[['CSBID_smaller', 'CSBID_larger', 'overlap_fraction']]
    
    # Clean CSB fields to remove overlapping fields
    # Extract all smaller parcel IDs from the pairs
    overlap_smaller_CSBIDs = CSB_overlaps["CSBID_smaller"].unique()
    # Drop rows whose CSBID is in smaller_ids and reset index
    CSB_geometry = (CSB_geometry[~CSB_geometry["CSBID"].isin(overlap_smaller_CSBIDs)]).reset_index(drop=True)
        
    
    
    ### Process mukey raster file to create pairs of CSBID CSBID: mukeys (pixel values) ###
    # CSBID CSBID → unique field integers
    id_map = {s: i for i, s in enumerate(CSB_geometry["CSBID"].unique())}
    # Unique field integers → CSBID CSBID
    reverse_id_map = {v: k for k, v in id_map.items()} 
    CSB_geometry["pid"] = CSB_geometry["CSBID"].map(id_map)

    # Open raster
    with rasterio.open(f"data/edited/Soil/gSSURGO Mukey Grid/{state}_MURASTER_30m.tif") as src:
        # Match CRS between vector CSBID and gSSURGO raster
        CSB_geometry = CSB_geometry.to_crs(src.crs)

        # Read gSSURGO raster file
        band = src.read(1)
        transform = src.transform
        height = src.height
        width = src.width

        # Rasterize fields: assign unique field integers to pixels whose centers lie inside a given polygon
        shapes = ((geom, pid) for geom, pid in zip(CSB_geometry.geometry, CSB_geometry.pid))
        
        parcel_raster = rasterize(
            shapes=shapes,
            out_shape=(height, width),
            transform=transform,
            fill=-1,            # pixels not belonging to any parcel
            all_touched=False,  # only pixels whose center is inside polygon
            dtype="int32",
            skip_invalid = False
        )

    # Extract pixel values for each field
    # Indices of all pixels that belong to some field
    rows, cols = np.where(parcel_raster >= 0)
    parcel_ids = parcel_raster[rows, cols]
    values = band[rows, cols].astype(str)

    # Build dictionary of CSBID → pixel values
    CSBID_mukeys = {CSBID: [] for CSBID in CSB_geometry.CSBID}
    for pid, val in zip(parcel_ids, values):
        CSBID_mukeys[reverse_id_map[pid]].append(val)
        
    # Find and print fields with no pixels assigned
    print("CSBID fields with no pixels assigned \n")
    for k, v in CSBID_mukeys.items():
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
    chorizon = chorizon[['cokey', 'hzdept_r', 'hzdepb_r', 'sandtotal_r', 'claytotal_r', 'ph1to1h2o_r']]
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
    #mapunit_component_chorizon_corerestictions_muaggatt_legend_sacatalog.to_parquet(output_path_integrated_soil, compression="zstd")
    
    
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
            "ph1to1h2o_r_30cm": horizon_weighted_30(x, "ph1to1h2o_r")
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
            "ph1to1h2o_r_30cm_weighted": component_weighted_mean(x, "ph1to1h2o_r_30cm")
        }))
        .reset_index()
    )

    # Add weighted variables variables to the integrated dataset
    mukey_soil_variables = mukey_soil_variables.merge(mukey_weighted, on="mukey", how="left")

    # Add slope data to the integrated dataset
    mukey_soil_variables = mukey_soil_variables.merge(muaggatt, on="mukey", how="left")



    ### Compute soil variables for each CSBID field ###
    # Split numerical and categorical columns
    num_cols = mukey_soil_variables.select_dtypes(include='number').columns.drop(['mukey', 'lkey'], errors='ignore')
    cat_cols = mukey_soil_variables.select_dtypes(exclude='number').columns.drop(['mukey', 'lkey'], errors='ignore')

    # Index soil and target datasets to speed up searching
    soil_indexed = mukey_soil_variables.set_index("mukey")
    CSB_soil = CSB_soil.set_index("CSBID")

    # Compute soil variables for each CSBID field by taking average (mode) across all mukeys within a given field
    for fieldID, mukeys in CSBID_mukeys.items():
        if mukeys:
            sub = soil_indexed.reindex(mukeys).dropna(how="all")
            if sub.empty:
                continue

            num_mean = sub[num_cols].mean()
            cat_mode = sub[cat_cols].apply(lambda x: x.dropna().value_counts().idxmax() if not x.dropna().empty else None)

            row = pd.concat([num_mean, cat_mode])

            CSB_soil.loc[fieldID, row.index] = row.values

    # Map each smaller_CSBID to its corresponding larger_CSBID and then copy the values from the larger fields back into the target dataset
    CSB_soil.loc[CSB_overlaps["CSBID_smaller"]] = CSB_soil.loc[CSB_overlaps["CSBID_larger"]].values
    CSB_soil = CSB_soil.reset_index()
    
    ### Save output files ###
    # Convert all float64 to float32
    float64_cols = CSB_soil.select_dtypes(include="float64").columns
    CSB_soil[float64_cols] = CSB_soil[float64_cols].astype("float32")

    CSB_soil.to_parquet(output_path_table, compression="zstd")
    print(f"Creating and saving soil dataset for {state} is complete")