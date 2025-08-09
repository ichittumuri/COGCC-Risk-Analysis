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
from scipy.signal import find_peaks
from scipy.ndimage import gaussian_filter1d
import matplotlib.colors as mcolors
import contextily as ctx
from shapely.geometry import box
from matplotlib.colors import ListedColormap
import matplotlib

# =============================================================================
# 2. Load Data
# =============================================================================
os.chdir('/Users/ichittumuri/Desktop/MINES/COGCC-Risk-Analysis/Data')
warnings.filterwarnings("ignore", category=RuntimeWarning)

flowlines_gdf = gpd.read_file("interpolated_clean_flowlines.geojson")
spills_gdf = gpd.read_file("clean_spills.geojson")

# =============================================================================
# 3. Match spills to flowlines
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
            print(f"[Spill {idx}] Missing geometry – skipping.")
            missing_geom_spill += 1
            continue

        op = spill.get("operator_name", "")
        if pd.isnull(op) or not op.strip():
            print(f"[Spill {idx}] No operator_name – skipping.")
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
            print(f"[Spill {idx}] No flowlines for operator “{op}”.")
            no_match_spill += 1
            continue

        dists = candidates.geometry.distance(pt).dropna()
        if dists.empty:
            print(f"[Spill {idx}] All distances NaN – skipping.")
            no_match_spill += 1
            continue

        nearest_idx  = dists.idxmin()
        min_dist     = dists.min()
        nearest_line = flowlines_gdf.loc[nearest_idx]

        _, nearest_pt = nearest_points(pt, nearest_line.geometry)

        matched_spill_count += 1
        print(f"[Match {matched_spill_count}] spill {idx} → flowline {nearest_idx} at {min_dist:.2f} m")

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
# 4. Plot simple histograms
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

plot_match_distance_hist(matched_spills_gdf, bins=30, title="Distance Between Original and Matched Spill Points")
plot_match_distance_hist(matched_spills_gdf, bins=3000, zoom_max=1000,
                         title="Distance Between Original and Matched Spill Points (Zoomed In)")

# =============================================================================
# 5. Find threshold distance between modes
# =============================================================================
def find_threshold_distance(
    gdf,
    column="match_distance_m",
    bins=50,
    sigma=1.5,
    prominence=5,
    count_cutoff=30
):
    distances = gdf[column].dropna()
    distances = distances[distances > 0]
    if distances.empty:
        return None, {"reason": "no_positive_distances"}

    mean_dist = float(distances.mean())
    median_dist = float(distances.median())

    min_val = float(distances.min())
    max_val = float(distances.max())
    if min_val <= 0 or max_val <= 0 or min_val == max_val:
        bins_array = bins
    else:
        bins_array = np.logspace(np.log10(min_val), np.log10(max_val), bins)

    counts, bin_edges = np.histogram(distances, bins=bins_array)
    counts_smooth = gaussian_filter1d(counts.astype(float), sigma=sigma)

    peak_idx, _ = find_peaks(counts_smooth, prominence=prominence)
    peak_idx = np.sort(peak_idx)

    split_index = None
    last_idx = None
    threshold_value = None

    if len(peak_idx) >= 2:
        first_mode, second_mode = int(peak_idx[0]), int(peak_idx[1])
        valley_slice = slice(first_mode, second_mode)
        valley_rel = int(np.argmin(counts_smooth[valley_slice]))
        split_index = first_mode + valley_rel

        for i in range(split_index - 1, -1, -1):
            if counts[i] > count_cutoff:
                last_idx = i
                threshold_value = float(bin_edges[last_idx + 1])
                break

    details = {
        "counts": counts,
        "bin_edges": bin_edges,
        "counts_smooth": counts_smooth,
        "peak_idx": peak_idx,
        "mean": mean_dist,
        "median": median_dist,
        "min_val": min_val,
        "max_val": max_val,
        "split_index": split_index,
        "last_idx": last_idx
    }
    return threshold_value, details

def plot_distance_histogram(
    gdf,
    details,
    threshold_value=None,
    column="match_distance_m",
    figsize=(10, 6),
    color="skyblue",
    count_cutoff=None
):
    distances = gdf[column].dropna()
    distances = distances[distances > 0]

    plt.figure(figsize=figsize)
    plt.hist(distances, bins=details["bin_edges"], edgecolor="black", color=color)

    plt.axvline(details["mean"], color='red', linestyle='--', label=f'Mean: {details["mean"]:.1f} m')
    plt.axvline(details["median"], color='green', linestyle='--', label=f'Median: {details["median"]:.1f} m')

    if threshold_value is not None:
        plt.axvline(threshold_value, color='blue', linestyle='--', linewidth=2,
                    label=f'Threshold: {threshold_value:.2f} m')

    if details["min_val"] > 0 and details["min_val"] != details["max_val"]:
        plt.xscale("log")

    plt.xlabel("Distance (meters, log scale)")
    plt.ylabel("Frequency")
    plt.grid(axis='y')
    plt.legend()

    title_parts = ["Log Distance Between Original and Matched Spill Points"]
    plt.title(" ".join(title_parts))
    plt.tight_layout()
    plt.show()

# Threshold
threshold, details = find_threshold_distance(
    matched_spills_gdf,
    bins=50,
    sigma=1.5,
    prominence=5,
    count_cutoff=30
)

if threshold is not None:
    print(f"Threshold distance: {threshold:.2f} m")
else:
    print("No threshold found before the second mode.")

# Plot log histogram
plot_distance_histogram(
    matched_spills_gdf,
    details,
    threshold_value=threshold,
    count_cutoff=30
)

# =============================================================================
# 6. Save spills under threshold
# =============================================================================
def save_spills_under_threshold(gdf, threshold, output_file):
    filtered_gdf = gdf[gdf["match_distance_m"] < threshold].copy()
    filtered_gdf = filtered_gdf.to_crs(epsg=4326)
    filtered_gdf.to_file(output_file, driver="GeoJSON")
    print(f"Wrote {output_file} in EPSG:4326.")
    return filtered_gdf

spills_under_85m = save_spills_under_threshold(
    matched_spills_gdf,
    threshold=85.07,
    output_file="updated_spills_under_85m.geojson"
)

# =============================================================================
# 7. Plot map of spills & flowlines
# =============================================================================
spills_under_85m = spills_under_85m[spills_under_85m.geometry.type.isin(['Point', 'MultiPoint'])].copy()

def plot_flowlines_and_spills(flowlines_gdf, spills_gdf, distance_col='match_distance_m', zoom_bounds=None):
    flowlines = flowlines_gdf.to_crs(epsg=3857)
    spills = spills_gdf.to_crs(epsg=3857)

    fig, ax = plt.subplots(figsize=(12, 12))

    flowlines.plot(
        ax=ax,
        color="#0B50B8",
        linewidth=2,
        alpha=0.6,
        zorder=2
    )

    distances = spills[distance_col]
    norm = mcolors.Normalize(vmin=distances.min(), vmax=distances.max())
    cmap = matplotlib.colormaps["YlOrRd"]

    spills.plot(
        ax=ax,
        column=distance_col,
        cmap=cmap,
        norm=norm,
        markersize=20,
        alpha=0.85,
        legend=True,
        zorder=3,
        legend_kwds={
            'label': "Distance to Flowline (m)",
            'shrink': 0.5
        }
    )

    ctx.add_basemap(ax, source=ctx.providers.CartoDB.Positron, zorder=1)

    if zoom_bounds:
        x_min, x_max, y_min, y_max = zoom_bounds
        ax.set_xlim(x_min, x_max)
        ax.set_ylim(y_min, y_max)

    ax.set_title("Spills Relocated to Nearest Flowline, Colored by Match Distance (≤85.07m)", fontsize=16)
    ax.set_axis_off()
    plt.tight_layout()
    plt.show()

# Call the function
plot_flowlines_and_spills(
    flowlines_gdf=flowlines_gdf,
    spills_gdf=spills_under_85m,
    distance_col='match_distance_m'
)

# =============================================================================
# 8. Plot cropped map (Boulder–Greeley area)
# =============================================================================
def plot_cropped_spills(
    flowlines_gdf,
    spills_gdf,
    bbox_wgs84,  # (minx, miny, maxx, maxy) in EPSG:4326
    distance_col="match_distance_m",
    title="Spills ≤ 85.07 m — Cropped View",
    out_png=None
):
    # bbox geodataframe in Web Mercator
    bbox_poly = box(*bbox_wgs84)
    bbox_gdf = gpd.GeoDataFrame(geometry=[bbox_poly], crs="EPSG:4326").to_crs(epsg=3857)
    x_min, y_min, x_max, y_max = bbox_gdf.total_bounds
    bbox_geom = bbox_gdf.geometry.iloc[0]

    # Ensure spills are point/multipoint only
    spills_pts = spills_gdf[spills_gdf.geometry.type.isin(["Point", "MultiPoint"])].copy()
    if spills_pts.empty:
        print("No point spills to plot.")
        return

    # Project + clip
    flowlines_m = flowlines_gdf.to_crs(epsg=3857).clip(bbox_geom)
    spills_m    = spills_pts.to_crs(epsg=3857).clip(bbox_geom)

    if flowlines_m.empty and spills_m.empty:
        print("No features within the provided bbox after clipping.")
        return

    # Truncated YlOrRd colormap
    base = cmap = matplotlib.colormaps["YlOrRd"]
    truncated_cmap = ListedColormap(base(np.linspace(0.25, 1.0, 256)))

    # Plot
    fig, ax = plt.subplots(figsize=(10, 10))

    if not flowlines_m.empty:
        flowlines_m.plot(ax=ax, color="#0B50B8", linewidth=2, alpha=0.6, zorder=2)

    if not spills_m.empty:
        d = spills_m[distance_col]
        norm = mcolors.Normalize(vmin=float(d.min()), vmax=float(d.max()))
        spills_m.plot(
            ax=ax,
            column=distance_col,
            cmap=truncated_cmap,
            norm=norm,
            markersize=8,
            alpha=0.85,
            legend=True,
            zorder=3,
            legend_kwds={"label": "Distance to Flowline (m)", "shrink": 0.5}
        )

    ctx.add_basemap(ax, source=ctx.providers.CartoDB.Positron, zorder=1)

    ax.set_xlim(x_min, x_max)
    ax.set_ylim(y_min, y_max)
    ax.set_title(title, fontsize=14, color="#444444")
    ax.set_axis_off()
    plt.tight_layout()
    plt.show()

bbox = (-106.0, 39.3, -103.9, 40.9)  # WGS84
plot_cropped_spills(
    flowlines_gdf=flowlines_gdf,
    spills_gdf=spills_under_85m,
    bbox_wgs84=bbox,
    title="Spills Relocated to Nearest Flowline (≤85.07m) — Boulder–Greeley Area")

# =============================================================================
# 9. Dataset summary
# =============================================================================
def summarize_spill_datasets(spills_file, matched_file, matched_under_file):
    spills_raw = gpd.read_file(spills_file)
    matched_spills = gpd.read_file(matched_file)
    matched_under = gpd.read_file(matched_under_file)

    n_raw = len(spills_raw)
    n_matched = len(matched_spills)
    n_under = len(matched_under)

    print(f"Original spill points ({spills_file}):               {n_raw}")
    print(f"Matched spills ({matched_file}):                     {n_matched}")
    print(f"Matched spills under threshold ({matched_under_file}): {n_under}")
    print(f"Retention rate under threshold:                      {n_under / n_raw:.2%}")

    return {
        "n_raw": n_raw,
        "n_matched": n_matched,
        "n_under": n_under,
        "retention_rate": n_under / n_raw
    }

summary = summarize_spill_datasets(
    "spills.geojson",
    "updated_spills.geojson",
    "updated_spills_under_85m.geojson"
)

# =============================================================================
# 10. Join spills with flowline attributes
# =============================================================================
joined_rows = []

match_idx_col = "matched_flowline_idx"
if match_idx_col not in spills_under_85m.columns:
    match_idx_col = "matched_crudeoil_idx"

for i, spill in spills_under_85m.iterrows():
    line_idx = spill.get(match_idx_col, None)
    if pd.isna(line_idx) or line_idx not in flowlines_gdf.index:
        print(f"[Skip spill {i}] No valid match index: {line_idx}")
        continue

    line = flowlines_gdf.loc[line_idx]
    new_row = line.drop(labels=["geometry"]).copy()

    for col in spill.index:
        if col in ["geometry", match_idx_col, "match_distance_m"]:
            continue
        new_row[f"spill_{col}"] = spill[col]

    new_row["spill_match_distance_m"] = spill.get("match_distance_m", None)
    new_row["geometry"] = spill.geometry
    joined_rows.append(new_row)

joined_gdf = gpd.GeoDataFrame(pd.DataFrame(joined_rows), geometry="geometry", crs=spills_under_85m.crs)

if joined_gdf.crs is None or joined_gdf.crs.to_epsg() != 4326:
    joined_gdf = joined_gdf.to_crs(epsg=4326)

joined_gdf["X"] = joined_gdf.geometry.x
joined_gdf["Y"] = joined_gdf.geometry.y
joined_gdf["coords"] = list(zip(joined_gdf["X"], joined_gdf["Y"]))

joined_gdf.to_file("spills_w_flowline_attributes.geojson", driver="GeoJSON")
print(f"Saved {len(joined_gdf)} spill-points with flowline attributes.")