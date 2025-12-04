# =============================================================================
# 1. Setup & Imports
# =============================================================================
import os
import warnings
import pandas as pd
import geopandas as gpd
from shapely.ops import nearest_points
import matplotlib.pyplot as plt
import numpy as np

# =============================================================================
# 2. Load Data
# =============================================================================
os.chdir('/Users/ichittumuri/Desktop/MINES/COGCC-Risk-Analysis/Data')
warnings.filterwarnings("ignore", category=RuntimeWarning)

flowlines_gdf = gpd.read_file("interpolated_clean_flowlines.geojson")
spills_gdf = gpd.read_file("clean_spills.geojson")

# =============================================================================
# 3. Match spills to flowlines  (NO THRESHOLD)
# =============================================================================
def match_spills_to_flowlines(spills_gdf, flowlines_gdf, output_file):

    matched_spills      = []
    matched_spill_count = 0
    missing_geom_spill  = 0
    no_op_spill         = 0
    no_match_spill      = 0

    projected_crs = spills_gdf.estimate_utm_crs()

    for idx, spill in spills_gdf.iterrows():
        pt = spill.geometry
        if pt is None or pt.is_empty:
            missing_geom_spill += 1
            continue

        op = spill.get("operator_name", "")
        if pd.isnull(op) or not op.strip():
            no_op_spill += 1
            continue

        op_name = op.strip().lower()
        candidates = flowlines_gdf[
            flowlines_gdf["operator_name"]
                .astype(str)
                .str.strip()
                .str.lower()
                .eq(op_name)
        ].copy()
        candidates = candidates[candidates.geometry.notnull()]

        if candidates.empty:
            no_match_spill += 1
            continue

        dists = candidates.geometry.distance(pt).dropna()
        if dists.empty:
            no_match_spill += 1
            continue

        nearest_idx  = dists.idxmin()
        min_dist     = dists.min()
        nearest_line = flowlines_gdf.loc[nearest_idx]

        _, nearest_pt = nearest_points(pt, nearest_line.geometry)

        matched_spill_count += 1
        new_spill = spill.copy()
        new_spill["match_point"]          = nearest_pt
        new_spill["match_distance_m"]     = min_dist
        new_spill["matched_flowline_idx"] = nearest_idx
        matched_spills.append(new_spill)

    print("\n=== Spills Matching Complete ===")
    print(f"Total spills:              {len(spills_gdf)}")
    print(f"Matched spills:            {matched_spill_count}")
    print(f"Missing geometry spills:   {missing_geom_spill}")
    print(f"No Operator Name spills:   {no_op_spill}")
    print(f"No line match spills:      {no_match_spill}")

    matched_spills_gdf = gpd.GeoDataFrame(matched_spills)
    matched_spills_gdf = matched_spills_gdf.rename(columns={"geometry": "orig_geometry"})
    matched_spills_gdf["geometry"] = matched_spills_gdf["match_point"]
    matched_spills_gdf = matched_spills_gdf.set_geometry("geometry")

    matched_spills_gdf = matched_spills_gdf.set_crs(projected_crs, allow_override=True)
    matched_spills_gdf = matched_spills_gdf.to_crs(epsg=4326)
    matched_spills_gdf = matched_spills_gdf.drop(columns=["orig_geometry", "match_point"])

    matched_spills_gdf.to_file(output_file, driver="GeoJSON")
    print(f"Wrote {output_file} in EPSG:4326.")

    return matched_spills_gdf


matched_spills_gdf = match_spills_to_flowlines(
    spills_gdf,
    flowlines_gdf,
    "updated_spills.geojson"
)

# =============================================================================
# 4. Plot simple histograms (OPTIONAL)
# =============================================================================
def plot_match_distance_hist(matched_spills_gdf, bins=30, zoom_max=None, figsize=(10, 6), title=None): 
    distances = matched_spills_gdf["match_distance_m"].dropna()

    plt.figure(figsize=figsize)
    plt.hist(distances, bins=bins, edgecolor="black")
    if title:
        plt.title(title)
    plt.xlabel("Distance (meters)")
    plt.ylabel("Frequency")
    if zoom_max is not None:
        plt.xlim(0, zoom_max)
    plt.grid(True)
    plt.tight_layout()
    plt.show()

plot_match_distance_hist(matched_spills_gdf, bins=40)
plot_match_distance_hist(matched_spills_gdf, bins=6000, zoom_max=1000)

# =============================================================================
# 5. Join ALL matched spills with flowline attributes
# =============================================================================
joined_rows = []

match_idx_col = "matched_flowline_idx"

for i, spill in matched_spills_gdf.iterrows():
    line_idx = spill.get(match_idx_col)
    if pd.isna(line_idx) or line_idx not in flowlines_gdf.index:
        continue

    flowline_row = flowlines_gdf.loc[line_idx]
    flow_dict = flowline_row.drop(labels=["geometry"]).to_dict()

    spill_keep_cols = [c for c in spill.index if c not in {"geometry", match_idx_col}]
    spill_dict = spill[spill_keep_cols].to_dict()

    merged = {**flow_dict, **spill_dict}
    merged["geometry"] = spill.geometry

    joined_rows.append(merged)

joined_gdf = gpd.GeoDataFrame(joined_rows, geometry="geometry", crs=matched_spills_gdf.crs)
joined_gdf = joined_gdf.to_crs(epsg=4326)

# Drop unnecessary duplicates
joined_gdf["X"] = joined_gdf.geometry.x
joined_gdf["Y"] = joined_gdf.geometry.y
cols_to_drop = ["X", "Y"]
joined_gdf = joined_gdf.drop(columns=[c for c in cols_to_drop if c in joined_gdf.columns])

# Column order
desired_order = [
    "unique_id", "operator_name", "operator_number",
    "flowline_id", "location_id", "status",
    "flowline_action", "location_type", "fluid",
    "material", "diameter_in", "length_ft",
    "max_operating_pressure", "line_age_yr", "construct_date",
    "match_distance_m", "incident_date", "root_cause", "risk",
    "geometry"
]
current_cols = joined_gdf.columns.tolist()
ordered_cols = [col for col in desired_order if col in current_cols] + \
               [col for col in current_cols if col not in desired_order]

joined_gdf = joined_gdf[ordered_cols]

joined_gdf.to_file("spills_w_flowline_attributes.geojson", driver="GeoJSON")
print(f"Saved {len(joined_gdf)} spill-points with flowline attributes.")
