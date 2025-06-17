# --- Start of Script 1: 1.excel_to_gbd.ipynb ---
# This code transforms the original data formats into GeoJSON files, removing unnecessary attributes along the way.
# The gbd files are in CRS : EPSG-26913

# ----------------------------------------
# 0. Import Libraries
# ----------------------------------------
import os
import pandas as pd
import openpyxl
import geopandas as gpd
from shapely.geometry import Point, LineString

os.chdir('/Users/ichittumuri/Desktop/MINES/COGCC-Risk-Analysis/Data')
pd.options.display.max_columns = None

import os
import pandas as pd
import geopandas as gpd
import shapely
from shapely.geometry import Point, Polygon, MultiLineString, LineString,MultiPolygon, MultiPoint
from shapely.ops import nearest_points
os.chdir('/Users/ichittumuri/Desktop/MINES/COGCC-Risk-Analysis/Data')
pd.options.display.max_columns = None

# ----------------------------------------
# 1. excel_to_gbd
# ----------------------------------------

# Load Data
gdf0 = gpd.GeoDataFrame(gpd.read_file("ECMC_Flowline_Data_Access/COGCC_Form44_Off_Location_Flowlines_Approved_CONFIDENTIAL.gdb"))
gdf1 = gpd.GeoDataFrame(gpd.read_file("ECMC_Flowline_Data_Access/COGCC_Form44_Crude_Oil_Produced_Water_Transfer_Flowlines_Approved_CONFIDENTIAL.gdb"))
flowlines = pd.read_excel('FlowlineSpreadsheet_Mines.xlsx')

# Check if the CRS of gdf0 is different from that of gdf1
if gdf0.crs != gdf1.crs:
    # If true, convert the CRS of gdf1 to match that of gdf0, modifying gdf1 in place
    gdf1.to_crs(gdf.crs, inplace=True)

# Concatenate gdf0 and gdf1 into a single GeoDataFrame, ignoring the original indices to create a new continuous index
gdf = pd.concat([gdf0,gdf1],ignore_index=True)
gdf.drop(columns=['Doc_Num'], inplace=True)

# List all columns in the 'flowlines' DataFrame
columns_list = flowlines.columns.tolist()
print("Columns in spills DataFrame:", columns_list)

# Keep the following columns
selected_columns = ['LOCATION_ID', 'FLOWLINEID', 'STARTLOCATIONID', 'FLOWLINEACTION', 'ENTIRELINEREMOVED', 'ACTIONDESCRIPTION', 'RECEIVE_DATE', 'OPERATOR_NUM', 'COMPANY_NAME', 'LOCATIONTYPE', 'ENDLAT', 'ENDLONG', 'STARTLAT', 'STARTLONG',
                    'PIPEMATERIAL', 'BEDDINGMATERIAL', 'TYPEOFFLUIDTRANS', 'MAXOPPRESSURE', 'CONSTRUCTDATE']
flowlines_selected = flowlines[selected_columns]

# Check if the 'geometry' column does not exist in the flowlines DataFrame and initialize it as an empty string if true
if 'geometry' not in flowlines_selected.columns:
    flowlines_selected['geometry'] = ''

# Iterate over each row in the flowlines DataFrame
for index, row in flowlines_selected.iterrows():
    # Create a LineString geometry from the start and end coordinates of each flowline
    geom = LineString([(row['STARTLONG'],row['STARTLAT']),(row['ENDLONG'],row['ENDLAT'])])
    # Assign the created LineString geometry to the 'geometry' column at the current index
    flowlines_selected.at[index,'geometry'] = geom

# Convert the flowlines DataFrame into a GeoDataFrame, explicitly setting the 'geometry' column and the Coordinate Reference System (CRS) to 'EPSG:4326' (WGS 84)
fl_gdf = gpd.GeoDataFrame(flowlines_selected, geometry='geometry', crs='EPSG:4326')

# Check if the CRS of fl_gdf is different from that of gdf
if fl_gdf.crs != gdf.crs:
    # If true, convert the CRS of fl_gdf to match that of gdf, modifying fl_gdf in place
    fl_gdf.to_crs(gdf.crs, inplace=True)
    print('Change fl crs to gdf crs')

# List all columns in the 'spills' DataFrame
spills = pd.read_excel('Flowline-Related Spills (through 2024).xlsx')
columns_list = spills.columns.tolist()
print("Columns in spills DataFrame:", columns_list)

# Keep the following columns
selected_columns = ['trkg_num', 'Operator Name', 'facility_type', 'Spill_Desc', 'Spill Type', 'Root Cause', 'Preventative Measure',
                    'Root Cause Type', 'Detailed Root Cause Type', 'Long', 'Lat','facility_status', 'Gathering?', 'Metallic?', 'incident_date']
spills_selected = spills[selected_columns]

# Number of spills
num_rows = len(spills_selected)
print("Number of rows:", num_rows)

# Remove rows where the 'Gathering?' column has information (non-NaN values)
spills_cleaned = spills_selected[spills_selected['Gathering?'].isna()]

# Check if all gathering intances are removed
only_nan_columns = spills_cleaned.isna().all()
columns_with_all_nan = only_nan_columns[only_nan_columns].index.tolist()
print("Columns with all NaN values:", columns_with_all_nan)

# Remove rows where 'Operator Name' contains the word "GATHERING"
spills_filtered = spills_cleaned[~spills_cleaned['Operator Name'].str.contains("GATHERING", case=False, na=False)]

# Remove 'Gathering?' and 'Operator Name' columns from the DataFrame
final_spills = spills_filtered.drop(columns=['Gathering?'])

# Number of spills after cleaning
num_rows = len(final_spills)
print("New number of rows:", num_rows)

# Convert 'spills' DataFrame to GeoDataFrame 'spl_gdf' with Point geometries from 'Long' and 'Lat', setting CRS to 'EPSG:4326'.
spl_gdf = gpd.GeoDataFrame(final_spills, geometry=gpd.points_from_xy(final_spills.Long,final_spills.Lat), crs='EPSG:4326')

# Check if the CRS of fl_gdf is different from that of gdf
if spl_gdf.crs != gdf.crs:
    # If true, convert the CRS of fl_gdf to match that of gdf, modifying fl_gdf in place
    spl_gdf.to_crs(gdf.crs, inplace=True)
    print('Change spl crs to gdf crs')

gdf.to_file('crudeoil_offlocation.geojson', driver='GeoJSON')
fl_gdf.to_file('flowlines.geojson', driver='GeoJSON')
spl_gdf.to_file('spills.geojson', driver='GeoJSON')

# ----------------------------------------
# 2. combine_all_flowlines
# ----------------------------------------

# Load Data
flowlines_gdf = gpd.read_file('flowlines.geojson')
crudeoil_gdf = gpd.read_file('crudeoil_offlocation.geojson')

# Check if CRS is the same for both files
if flowlines_gdf.crs != crudeoil_gdf.crs:
    flowlines_gdf = flowlines_gdf.to_crs(crudeoil_gdf.crs)

clean_crudeoil_gdf = crudeoil_gdf.dropna()

import geopandas as gpd
import pandas as pd
from shapely.ops import nearest_points

# 1. Load your GeoJSONs
flowlines_gdf = gpd.read_file("flowlines.geojson")
crudeoil_gdf  = gpd.read_file("crudeoil_offlocation.geojson")

# 2. Reproject if needed
if crudeoil_gdf.crs != flowlines_gdf.crs:
    crudeoil_gdf = crudeoil_gdf.to_crs(flowlines_gdf.crs)

# 3. Prepare counters and storage
matched_rows   = []
matched_count  = 0
skipped_count  = 0

# 4. Loop over the entire crude-oil dataset
for idx, crude in crudeoil_gdf.iterrows():
    geom = crude.geometry
    if geom is None:
        print(f"[{idx}] Missing geometry – skipping.")
        skipped_count += 1
        continue

    op = crude.get("Operator", "")
    if pd.isnull(op) or not op.strip():
        print(f"[{idx}] No Operator – skipping.")
        skipped_count += 1
        continue

    op_name = op.strip().lower()
    # Filter flowlines by matching operator name and drop null geometries
    candidates = flowlines_gdf[
        flowlines_gdf["COMPANY_NAME"]
            .str.strip()
            .str.lower()
            .eq(op_name)
    ].copy()
    candidates = candidates[candidates.geometry.notnull()]

    if candidates.empty:
        print(f"[{idx}] No flowlines for operator “{op}”.")
        skipped_count += 1
        continue

    # Compute true geometry distances
    dists = candidates.geometry.distance(geom).dropna()
    if dists.empty:
        print(f"[{idx}] All distances NaN – skipping.")
        skipped_count += 1
        continue

    nearest_idx = dists.idxmin()
    min_dist    = dists.min()

    matched_count += 1
    print(f"[Match {matched_count}] crudeoil index {idx} → flowline index {nearest_idx} at {min_dist:.2f} m")

    # Build a new row: keep crude-oil geometry, inject flowline attributes (except geometry)
    new_row = crude.copy()
    for col in flowlines_gdf.columns:
        if col == flowlines_gdf.geometry.name:
            continue
        new_row[col] = flowlines_gdf.loc[nearest_idx, col]
    new_row["flowline_match_distance_m"] = min_dist

    matched_rows.append(new_row)

# 5. Summary
total = len(crudeoil_gdf)
print("\n=== Complete ===")
print(f" Total crude-oil features:         {total}")
print(f" Successfully matched:             {matched_count}")
print(f" Skipped (no match or missing data): {skipped_count}")

# 6. Create a GeoDataFrame of the matched crude-oil features
matched_crudeoil_gdf = gpd.GeoDataFrame(matched_rows, crs=crudeoil_gdf.crs)

# 8. Inspect the first few rows
print(matched_crudeoil_gdf.head())

# 5. Summary
total = len(crudeoil_gdf)
print("\n=== Complete ===")
print(f" Total crude-oil features:         {total}")
print(f" Successfully matched:             {matched_count}")
print(f" Skipped (no match or missing data): {skipped_count}")

# Create a simple 1,2,3… unique identifier
matched_crudeoil_gdf["unique_id"] = range(1, len(matched_crudeoil_gdf) + 1)

# 2. Save back out to GeoJSON
matched_crudeoil_gdf.to_file(
    "full_length_flowlines.geojson",
    driver="GeoJSON"
)