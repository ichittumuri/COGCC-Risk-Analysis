import os
import pandas as pd
import geopandas as gpd
import numpy as np
from shapely.geometry import MultiLineString, LineString, Point
from shapely.ops import nearest_points
from multiprocessing import Pool, cpu_count
from sklearn.neighbors import BallTree


os.chdir('/Users/ichittumuri/Desktop/MINES/COGCC-Risk-Analysis/Data')
pd.options.display.max_columns = None

# Load Data
flowlines_gdf = gpd.read_file('split_flowlines.geojson')
spills_gdf = gpd.read_file('spills.geojson')

# Check and count missing geometries in spills
missing_geometry_count = spills_gdf[spills_gdf.geometry.isna()].shape[0]
# Print the count of missing geometries
print(f"Number of spills with missing geometry: {missing_geometry_count}")

# Drop rows with missing geometries
spills_gdf = spills_gdf.dropna(subset=['geometry'])

