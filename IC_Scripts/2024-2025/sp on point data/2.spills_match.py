import os
import pandas as pd
import geopandas as gpd
import numpy as np
from shapely.geometry import MultiLineString, LineString, Point
from shapely.ops import nearest_points
from multiprocessing import Pool, cpu_count
from sklearn.neighbors import BallTree

import contextily as ctx
from matplotlib.colors import ListedColormap, Normalize
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt

# ------------------------------------------------------------------------------
# Set working directory
# ------------------------------------------------------------------------------
os.chdir('/Users/ichittumuri/Desktop/MINES/COGCC-Risk-Analysis/Data')

# ------------------------------------------------------------------------------
# 1. Load your augmented crude-oil lines and your raw spills
# ------------------------------------------------------------------------------
crudeoil_gdf = gpd.read_file("full_length_flowlines.geojson")
spills_gdf   = gpd.read_file("spills.geojson")

# ------------------------------------------------------------------------------
# 2. Reproject both to a projected CRS for accurate meter distances (auto-UTM)
# ------------------------------------------------------------------------------
projected_crs = spills_gdf.estimate_utm_crs()
crudeoil_gdf  = crudeoil_gdf.to_crs(projected_crs)
spills_gdf    = spills_gdf.to_crs(projected_crs)

# ------------------------------------------------------------------------------
# 3. Prepare counters and storage for matching
# ------------------------------------------------------------------------------
matched_spills      = []
matched_spill_count = 0
missing_geom_spill  = 0
no_op_spill         = 0
no_match_spill      = 0

# ------------------------------------------------------------------------------
# 4. Loop over every spill point and find nearest line of same operator
# ------------------------------------------------------------------------------
for idx, spill in spills_gdf.iterrows():
    pt = spill.geometry
    if pt is None or pt.is_empty:
        print(f"[Spill {idx}] Missing geometry – skipping.")
        missing_geom_spill += 1
        continue

    op = spill.get("Operator Name", "")
    if pd.isnull(op) or not op.strip():
        print(f"[Spill {idx}] No Operator Name – skipping.")
        no_op_spill += 1
        continue

    op_name = op.strip().lower()
    candidates = crudeoil_gdf[
        crudeoil_gdf["Operator"]
            .str.strip()
            .str.lower()
            .eq(op_name)
    ].copy()
    candidates = candidates[candidates.geometry.notnull()]

    if candidates.empty:
        print(f"[Spill {idx}] No crude-oil lines for operator “{op}”.")
        no_match_spill += 1
        continue

    dists = candidates.geometry.distance(pt).dropna()
    if dists.empty:
        print(f"[Spill {idx}] All distances NaN – skipping.")
        no_match_spill += 1
        continue

    nearest_idx = dists.idxmin()
    min_dist    = dists.min()
    nearest_line = crudeoil_gdf.loc[nearest_idx]

    # exact nearest point on the line
    _, nearest_pt = nearest_points(pt, nearest_line.geometry)

    matched_spill_count += 1
    print(f"[Match {matched_spill_count}] spill {idx} → crude-oil {nearest_idx} at {min_dist:.2f} m")

    # build a new row: copy all original spill attributes but add match info
    new_spill = spill.copy()
    new_spill["match_point"]         = nearest_pt
    new_spill["match_distance_m"]    = min_dist
    new_spill["matched_crudeoil_idx"] = nearest_idx

    matched_spills.append(new_spill)

# ------------------------------------------------------------------------------
# 5. Summary of matching
# ------------------------------------------------------------------------------
print("\n=== Spills Matching Complete ===")
print(f"Total spills:              {len(spills_gdf)}")
print(f"Matched spills:            {matched_spill_count}")
print(f"Missing geometry spills:   {missing_geom_spill}")
print(f"No Operator Name spills:   {no_op_spill}")
print(f"No line match spills:      {no_match_spill}")

# ------------------------------------------------------------------------------
# 6. Create GeoDataFrame of matched spills
# ------------------------------------------------------------------------------
matched_spills_gdf = gpd.GeoDataFrame(matched_spills, crs=spills_gdf.crs)

# ------------------------------------------------------------------------------
# 7. Swap in match_point as the new geometry
# ------------------------------------------------------------------------------
matched_spills_gdf = matched_spills_gdf.rename(columns={"geometry": "orig_geometry"})
matched_spills_gdf["geometry"] = matched_spills_gdf["match_point"]
matched_spills_gdf = matched_spills_gdf.set_geometry("geometry")

# 7b. Re-declare CRS so GeoPandas knows geometry is still in projected UTM
matched_spills_gdf = matched_spills_gdf.set_crs(projected_crs, allow_override=True)

# ------------------------------------------------------------------------------
# 8. Reproject matched spills back to WGS84 lon/lat
# ------------------------------------------------------------------------------
matched_spills_gdf = matched_spills_gdf.to_crs(epsg=4326)

# ------------------------------------------------------------------------------
# 9. Drop temporary columns if desired
# ------------------------------------------------------------------------------
matched_spills_gdf = matched_spills_gdf.drop(columns=["orig_geometry", "match_point"])

# ------------------------------------------------------------------------------
# 10. Write out updated_spills.geojson (in EPSG:4326)
# ------------------------------------------------------------------------------
matched_spills_gdf.to_file("updated_spills.geojson", driver="GeoJSON")
print("Wrote updated_spills.geojson in EPSG:4326.")

# ------------------------------------------------------------------------------
# 11. Load full_length_flowlines (in WGS84) and updated spills
# ------------------------------------------------------------------------------
crudeoil_wgs = gpd.read_file("full_length_flowlines.geojson").to_crs(epsg=4326)
spills_wgs   = gpd.read_file("updated_spills.geojson")  # already EPSG:4326

# ------------------------------------------------------------------------------
# 12. Build joined_rows by combining spill geometry + flowline attributes
# ------------------------------------------------------------------------------
joined_rows = []
for _, spill in spills_wgs.iterrows():
    line_idx = spill["matched_crudeoil_idx"]
    line     = crudeoil_wgs.loc[line_idx]

    # Start with all crudeoil attributes (excluding its geometry)
    new_row = line.drop(labels=["geometry"]).copy()

    # Add all spill attributes (excluding geometry, matched_crudeoil_idx, match_distance_m)
    for col in spill.index:
        if col in ["geometry", "matched_crudeoil_idx", "match_distance_m"]:
            continue
        new_row[f"spill_{col}"] = spill[col]

    # Explicitly add match distance
    new_row["spill_match_distance_m"] = spill["match_distance_m"]

    # Use spill geometry instead of line geometry
    new_row["geometry"] = spill.geometry

    joined_rows.append(new_row)

# ------------------------------------------------------------------------------
# 13. Create GeoDataFrame of spills with appended flowline attributes
# ------------------------------------------------------------------------------
joined_gdf = gpd.GeoDataFrame(
    pd.DataFrame(joined_rows),
    geometry="geometry",
    crs=spills_wgs.crs
)

# ------------------------------------------------------------------------------
# 14. Add risk column and set all values to 1
# ------------------------------------------------------------------------------
joined_gdf["risk"] = 1

# ------------------------------------------------------------------------------
# 15. Extract numeric X, Y, and coords from spill geometry
# ------------------------------------------------------------------------------
joined_gdf["X"]      = joined_gdf.geometry.x
joined_gdf["Y"]      = joined_gdf.geometry.y
joined_gdf["coords"] = list(zip(joined_gdf["X"], joined_gdf["Y"]))

# ------------------------------------------------------------------------------
# 16. Save the final spills GeoDataFrame (with flowline attributes, X, Y, coords, and risk)
# ------------------------------------------------------------------------------
joined_gdf.to_file(
    "spills_w_flowline_attributes.geojson",
    driver="GeoJSON"
)
print(f"Created {len(joined_gdf)} spill-points with appended flowline attributes, X/Y/coords, and risk=1.")

# ------------------------------------------------------------------------------
# 17. (Optional) Visualize distribution of match distances
# ------------------------------------------------------------------------------
plt.figure(figsize=(10, 6))
plt.hist(matched_spills_gdf["match_distance_m"].dropna(), bins=30, edgecolor="black")
plt.title("Distance Between Original and Matched Spill Points")
plt.xlabel("Distance (meters)")
plt.ylabel("Frequency")
plt.grid(True)
plt.tight_layout()
plt.show()

plt.figure(figsize=(10, 6))
plt.hist(matched_spills_gdf["match_distance_m"].dropna(), bins=30, edgecolor="black")
plt.title("Distance Between Original and Matched Spill Points (Zoomed In)")
plt.xlabel("Distance (meters)")
plt.ylabel("Frequency")
plt.xlim(0, 5000)
plt.grid(True)
plt.tight_layout()
plt.show()
