# plot_spills_flowlines.py
import os
import warnings
import geopandas as gpd
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import contextily as ctx
from shapely.geometry import box

# -----------------------------
# 1) Paths + basic setup
# -----------------------------
os.chdir('/Users/ichittumuri/Desktop/MINES/COGCC-Risk-Analysis/Data')
warnings.filterwarnings("ignore", category=UserWarning)

FLOWLINES_FP = "combined_flowlines_final.geojson"
SPILLS_FP    = "updated_spills_under_85m.geojson"
DIST_COL     = "match_distance_m"

# -----------------------------
# 2) Load
# -----------------------------
flowlines = gpd.read_file(FLOWLINES_FP)
spills = gpd.read_file(SPILLS_FP)
spills = spills[spills.geometry.type.isin(["Point", "MultiPoint"])].copy()

# -----------------------------
# 3) Full map
# -----------------------------
f_3857 = flowlines.to_crs(3857)
s_3857 = spills.to_crs(3857)

fig, ax = plt.subplots(figsize=(10, 10))
f_3857.plot(ax=ax, color="#0B50B8", linewidth=1, alpha=0.9, zorder=2)

d = s_3857[DIST_COL]
norm = mcolors.Normalize(vmin=float(d.min()), vmax=float(d.max()))
s_3857.plot(
    ax=ax,
    column=DIST_COL,
    cmap="YlOrRd",
    norm=norm,
    markersize=10,
    edgecolor="black",
    linewidth=0.2,
    legend=True,
    zorder=3,
    legend_kwds={"label": "Distance to Flowline (m)", "shrink": 0.6},
)

ctx.add_basemap(ax, source=ctx.providers.CartoDB.Positron, zorder=1)
# ax.set_title("Spills ≤ 85.07 m — Full Extent", fontsize=12)
ax.set_axis_off()
plt.tight_layout()
plt.savefig("spills_flowlines_full.png", dpi=200)
plt.show()

# -----------------------------
# 4) Cropped map (Boulder–Greeley)
# -----------------------------
bbox_wgs84 = (-106.0, 39.3, -103.9, 40.9)
bbox_gdf = gpd.GeoDataFrame(geometry=[box(*bbox_wgs84)], crs=4326).to_crs(3857)
xmin, ymin, xmax, ymax = bbox_gdf.total_bounds
clip_geom = bbox_gdf.geometry.iloc[0]

f_clip = f_3857.clip(clip_geom)
s_clip = s_3857.clip(clip_geom)

fig, ax = plt.subplots(figsize=(10, 10))
if not f_clip.empty:
    f_clip.plot(ax=ax, color="#0B50B8", linewidth=1, alpha=0.95, zorder=2)
if not s_clip.empty:
    d2 = s_clip[DIST_COL]
    norm2 = mcolors.Normalize(vmin=float(d2.min()), vmax=float(d2.max()))
    s_clip.plot(
        ax=ax,
        column=DIST_COL,
        cmap="YlOrRd",
        norm=norm2,
        markersize=10,
        edgecolor="black",
        linewidth=0.2,
        legend=True,
        zorder=3,
        legend_kwds={"label": "Distance to Flowline (m)", "shrink": 0.6},
    )

ctx.add_basemap(ax, source=ctx.providers.CartoDB.Positron, zorder=1)
ax.set_xlim(xmin, xmax)
ax.set_ylim(ymin, ymax)
# ax.set_title("Spills ≤ 85.07 m — Cropped (Boulder–Greeley)", fontsize=12)
ax.set_axis_off()
plt.tight_layout()
plt.savefig("spills_flowlines_cropped.png", dpi=200)
plt.show()
