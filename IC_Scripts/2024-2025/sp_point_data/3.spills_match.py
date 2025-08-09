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

os.chdir('/Users/ichittumuri/Desktop/MINES/COGCC-Risk-Analysis/Data')

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

# Log scale histogram:

# Drop NA and filter positive distances
distances = matched_spills_gdf["match_distance_m"].dropna()
distances = distances[distances > 0]

# Compute mean and median
mean_dist = distances.mean()
median_dist = distances.median()

# Log-spaced bins
min_val = distances.min()
max_val = distances.max()
bins = np.logspace(np.log10(min_val), np.log10(max_val), 50)

plt.figure(figsize=(10, 6))
plt.hist(distances, bins=bins, edgecolor="black", color="skyblue")

# Add vertical lines for mean and median
plt.axvline(mean_dist, color='red', linestyle='--', label=f'Mean: {mean_dist:.1f} m')
plt.axvline(median_dist, color='green', linestyle='--', label=f'Median: {median_dist:.1f} m')

plt.xscale("log")
plt.title("Distance Between Original and Matched Spill Points (Log X-Axis)")
plt.xlabel("Distance (meters, log scale)")
plt.ylabel("Frequency")
plt.grid(axis='y')  # Only horizontal lines
plt.legend()
plt.tight_layout()
plt.show()
