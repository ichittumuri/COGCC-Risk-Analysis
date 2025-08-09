import os
import pandas as pd
import numpy as np
import geopandas as gpd
from shapely.geometry import MultiLineString, LineString, Point
from shapely.ops import nearest_points, transform
import matplotlib.pyplot as plt
import seaborn as sns
import pyproj


# Set working directory and display options
os.chdir('/Users/ichittumuri/Desktop/MINES/COGCC-Risk-Analysis/Data')
pd.options.display.max_columns = None

# ----------------------------------------------------------------------
# 1. Load flowlines, add 'risk' column, and compute global min/max lengths
# ----------------------------------------------------------------------

matched_flowlines_gdf = gpd.read_file('full_length_flowlines.geojson')
matched_flowlines_gdf["risk"] = 0  # start all flowlines with risk = 0

# Ensure geometries are in geographic CRS (WGS84) for geodesic length calculations
matched_flowlines_gdf = matched_flowlines_gdf.to_crs(epsg=4326)

# Create a Geod object for geodesic calculations
geod = pyproj.Geod(ellps="WGS84")

global_min_length = float('inf')
global_max_length = float('-inf')

for geom in matched_flowlines_gdf['geometry']:
    if isinstance(geom, MultiLineString):
        for line in geom.geoms:
            length_m = geod.geometry_length(line)
            global_min_length = min(global_min_length, length_m)
            global_max_length = max(global_max_length, length_m)
    elif isinstance(geom, LineString):
        length_m = geod.geometry_length(geom)
        global_min_length = min(global_min_length, length_m)
        global_max_length = max(global_max_length, length_m)

print("Global min length (m):", global_min_length)
print("Global max length (m):", global_max_length)

# ----------------------------------------
# 2. Plot distribution of LineString lengths
# ----------------------------------------

def geodesic_length(line):
    return geod.geometry_length(line)

line_lengths = []
for geom in matched_flowlines_gdf['geometry']:
    if isinstance(geom, MultiLineString):
        for line in geom.geoms:
            line_lengths.append(geodesic_length(line))
    elif isinstance(geom, LineString):
        line_lengths.append(geodesic_length(geom))

plt.figure(figsize=(10, 6))
sns.histplot(line_lengths, bins=50, kde=True)
plt.title('Distribution of LineString Geodesic Lengths')
plt.xlabel('Length (meters)')
plt.ylabel('Frequency')
plt.grid(True)
plt.show()

max_length_to_plot = 500
filtered_lengths = [length for length in line_lengths if length <= max_length_to_plot]

plt.figure(figsize=(10, 6))
sns.histplot(filtered_lengths, bins=50, kde=True)
plt.title(f'Distribution of Geodesic LineString Lengths (0 to {max_length_to_plot} m)')
plt.xlabel('Length (meters)')
plt.ylabel('Frequency')
plt.grid(True)
plt.show()

# Logscale plot of LineString lengths

def geodesic_length(line):
    return geod.geometry_length(line)

line_lengths = []
for geom in matched_flowlines_gdf['geometry']:
    if isinstance(geom, MultiLineString):
        for line in geom.geoms:
            line_lengths.append(geodesic_length(line))
    elif isinstance(geom, LineString):
        line_lengths.append(geodesic_length(geom))

# Convert to NumPy array for stats
line_lengths = np.array(line_lengths)
line_lengths = line_lengths[line_lengths > 0]  # avoid log(0) issues

mean_len = np.mean(line_lengths)
median_len = np.median(line_lengths)

plt.figure(figsize=(10, 6))
sns.histplot(line_lengths, bins=50, kde=False, log_scale=(True, False), color='skyblue')

# Add vertical lines for mean and median
plt.axvline(mean_len, color='red', linestyle='--', label=f'Mean: {mean_len:.1f} m')
plt.axvline(median_len, color='green', linestyle='--', label=f'Median: {median_len:.1f} m')

plt.title('Distribution of LineString Geodesic Lengths (Log X-Axis)')
plt.xlabel('Length (meters, log scale)')
plt.ylabel('Frequency')
plt.grid(True, which='both', axis='x')
plt.legend()
plt.show()

# ----------------------------------------
# 3. Generate evenly spaced points along each flowline,
#    ensuring no duplicate points (with small offsets if needed)
# ----------------------------------------

def points_along_line(line, spacing):
    """
    Given a LineString 'line' and a spacing in the same units as line.length,
    return a list of Points at every 'spacing' distance along the line.
    """
    num_points = int(line.length // spacing) + 1
    return [line.interpolate(i * spacing) for i in range(num_points)]

def points_along_multiline(multiline, spacing):
    """
    Given a MultiLineString 'multiline' and a spacing, return all Points
    along each constituent LineString.
    """
    all_points = []
    for line in multiline.geoms:
        all_points.extend(points_along_line(line, spacing))
    return all_points

# Desired spacing in the same units as line.length (for geographic CRS, this is in degrees; adjust as needed)
point_spacing = 50  # e.g., 50 "units" along the geometry (if geometries are projected, units are meters)

seen_coords = {}  # to track how many times each exact (x, y) has been used
point_features = []  # will hold each point's new row

for idx, row in matched_flowlines_gdf.iterrows():
    geom = row.geometry

    if isinstance(geom, MultiLineString):
        raw_points = points_along_multiline(geom, point_spacing)
    elif isinstance(geom, LineString):
        raw_points = points_along_line(geom, point_spacing)
    else:
        continue  # skip non-line geometries

    for pt in raw_points:
        x, y = pt.x, pt.y
        coord_key = (x, y)

        if coord_key not in seen_coords:
            seen_coords[coord_key] = 0
            unique_pt = pt
        else:
            # Duplicate detected—shift by a tiny offset (in degrees) so it won’t be exactly equal
            seen_coords[coord_key] += 1
            offset_count = seen_coords[coord_key]
            offset_amount = 0.00001 * offset_count
            unique_pt = Point(x + offset_amount, y + offset_amount)
            coord_key = (unique_pt.x, unique_pt.y)
            seen_coords[coord_key] = 0

        # Create a new row, copying all attributes (including 'risk') but replacing geometry
        new_row = row.copy()
        new_row["geometry"] = unique_pt
        point_features.append(new_row)

points_gdf = gpd.GeoDataFrame(point_features, columns=matched_flowlines_gdf.columns, crs=matched_flowlines_gdf.crs)

# ----------------------------------------------------------------------
# 4. Extract X/Y into columns and build a 'coords' tuple column
# ----------------------------------------------------------------------

points_gdf['X'] = points_gdf.geometry.x
points_gdf['Y'] = points_gdf.geometry.y
points_gdf['coords'] = list(zip(points_gdf['X'], points_gdf['Y']))

# ----------------------------------------------------------------------
# 5. Drop duplicates, but ONLY remove rows where risk == 0
# ----------------------------------------------------------------------

# Split into two subsets:
#  - keep_always := all rows with risk != 0
#  - candidates := all rows with risk == 0
keep_always = points_gdf[points_gdf['risk'] != 0].copy()
candidates = points_gdf[points_gdf['risk'] == 0].copy()

# 5a) If any (X,Y) appear in keep_always, drop those same (X,Y) from candidates
coords_in_keep_always = set(zip(keep_always['X'], keep_always['Y']))
mask_to_drop_from_candidates = candidates[['X','Y']].apply(lambda r: (r['X'], r['Y']) in coords_in_keep_always, axis=1)
dropped_due_to_nonzero_risk_conflict = mask_to_drop_from_candidates.sum()
candidates = candidates.loc[~mask_to_drop_from_candidates].copy()

print(f"Dropped {dropped_due_to_nonzero_risk_conflict} rows where risk=0 but (X,Y) matched a risk≠0 point.")

# 5b) Now, within the remaining 'candidates' (all have risk=0), drop any true duplicates of (X, Y),
#      leaving exactly one row per (X, Y)
dupe_mask_candidates = candidates.duplicated(subset=['X','Y'])
num_dupes = dupe_mask_candidates.sum()
candidates = candidates.loc[~dupe_mask_candidates].copy()

print(f"Dropped {num_dupes} duplicate rows (all risk=0) within candidates so each (X,Y) appears just once.")

# 5c) Re‐assemble the final GeoDataFrame:
final_gdf = pd.concat([keep_always, candidates], ignore_index=True)

# ----------------------------------------------------------------------
# 6. Final sanity check
# ----------------------------------------------------------------------
all_coords = list(zip(final_gdf['X'], final_gdf['Y']))
unique_coords = set(all_coords)

print(f"Total rows in final_gdf: {len(final_gdf)}")
print(f"Unique (X, Y) pairs:       {len(unique_coords)}")
if len(final_gdf) == len(unique_coords):
    print("✔ Final (X,Y) pairs are all unique.")
else:
    print("⚠ Still duplicates present somewhere—please investigate!")

# ----------------------------------------------------------------------
# 7. Save to GeoJSON
# ----------------------------------------------------------------------
final_gdf.to_file("final_point_data.geojson", driver="GeoJSON")
