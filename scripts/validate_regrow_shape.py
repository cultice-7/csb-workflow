import geopandas as gp
import pandas as pd
import rasterio
import os
from rasterstats import zonal_stats

# Load the shapefile with polygons
gdf = gp.read_file("../data/edited/Regrow/regrow_clean.geojson")

# Folder containing the CDL rasters
raster_folder = "../data/CDL/"

# Initialize summary list
summary_data = []

# Loop through each year (2014 to 2024) and perform zonal statistics
for year in range(2014, 2025):
    raster_path = os.path.join(raster_folder, f"{year}_30m_cdls.tif")
    
    if os.path.exists(raster_path):
        print(f"Processing {year}...")
        
        # Compute zonal statistics
        stats = zonal_stats(
            vectors=gdf,
            raster=raster_path,
            stats="majority",
            geojson_out=False,
            nodata=0
        )
        
        # Add majority land cover value
        cdl_col = f"cdl_{year}"
        gdf[cdl_col] = [s.get("majority", None) for s in stats]

        # Dynamically check which crop columns exist
        crop_cols = [col for col in [f"crop{year}_1", f"crop{year}_2", f"crop{year}_3"] if col in gdf.columns]
        
        # Add validation column
        valid_col = f"valid_{year}"
        gdf[valid_col] = gdf.apply(
            lambda row: 1 if row[cdl_col] in [row.get(c) for c in crop_cols] else 0,
            axis=1
        )

        # Add to summary
        valid_count = int(gdf[valid_col].sum())
        total_count = int(gdf[valid_col].count())
        invalid_count = total_count - valid_count
        percent_valid = (valid_count / total_count) * 100 if total_count > 0 else 0

        summary_data.append({
            "year": year,
            "total_count": total_count,
            "valid_count": valid_count,
            "invalid_count": invalid_count,
            "percent_valid": percent_valid
        })
    else:
        print(f"Raster not found for year {year}: {raster_path}")

# Save the updated GeoDataFrame
output_geojson = "../data/edited/Regrow/regrow_with_cdl_validation.geojson"
gdf.to_file(output_geojson, driver="GeoJSON")
print(f"Saved updated GeoJSON to: {output_geojson}")

# Save the summary table
summary_df = pd.DataFrame(summary_data)
summary_csv = "../data/edited/Regrow/regrow_validity_summary_by_year.csv"
summary_df.to_csv(summary_csv, index=False)
print(f"Saved summary table to: {summary_csv}")
