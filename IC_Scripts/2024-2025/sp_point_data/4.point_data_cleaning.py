import os
import numpy as np
import pandas as pd
import geopandas as gpd
from sklearn.impute import KNNImputer
from sklearn.preprocessing import LabelEncoder

pd.options.display.max_columns = None

# set working directory
os.chdir('/Users/ichittumuri/Desktop/MINES/COGCC-Risk-Analysis/Data')

# 1. Load both GeoJSONs and reproject to WGS84
flowline_pts = gpd.read_file("final_point_data.geojson").to_crs(epsg=4326)
spills_pts   = gpd.read_file("spills_w_flowline_attributes.geojson").to_crs(epsg=4326)

# 2. Add explicit X, Y columns
flowline_pts["X"] = flowline_pts.geometry.x
flowline_pts["Y"] = flowline_pts.geometry.y
spills_pts["X"]   = spills_pts.geometry.x
spills_pts["Y"]   = spills_pts.geometry.y

# 3. Split each dataset into risk=1 vs. risk=0
flow_keep1   = flowline_pts[flowline_pts["risk"] == 1].copy()
flow_cand    = flowline_pts[flowline_pts["risk"] == 0].copy()
spills_keep1 = spills_pts[spills_pts["risk"] == 1].copy()
spills_cand  = spills_pts[spills_pts["risk"] == 0].copy()

# 4. Drop any flow_cand rows whose (X,Y) appear in flow_keep1
coords_in_flow_keep1 = set(zip(flow_keep1["X"].round(6), flow_keep1["Y"].round(6)))
mask_flow_drop       = flow_cand[["X", "Y"]].round(6).apply(
    lambda r: (r["X"], r["Y"]) in coords_in_flow_keep1,
    axis=1
)
flow_cand = flow_cand.loc[~mask_flow_drop].copy()

# 5. Drop duplicates within flow_cand
flow_cand = flow_cand.drop_duplicates(subset=["X", "Y"])

# 6. Drop any spills_cand rows whose (X,Y) appear in spills_keep1
coords_in_spills_keep1 = set(zip(spills_keep1["X"].round(6), spills_keep1["Y"].round(6)))
mask_spills_drop       = spills_cand[["X", "Y"]].round(6).apply(
    lambda r: (r["X"], r["Y"]) in coords_in_spills_keep1,
    axis=1
)
spills_cand = spills_cand.loc[~mask_spills_drop].copy()

# 7. Drop duplicates within spills_cand
spills_cand = spills_cand.drop_duplicates(subset=["X", "Y"])

# 8. Re‐assemble final cleaned GeoDataFrame
combined_gdf = pd.concat([flow_keep1, flow_cand, spills_keep1, spills_cand], ignore_index=True)
combined_gdf = gpd.GeoDataFrame(combined_gdf, geometry="geometry", crs="EPSG:4326")

# 9. ROUND coordinates and drop any remaining duplicates across all sources
combined_gdf["Xr"] = combined_gdf["X"].round(6)
combined_gdf["Yr"] = combined_gdf["Y"].round(6)
combined_gdf = combined_gdf.drop_duplicates(subset=["Xr", "Yr"])

# 10. Verify no duplicates remain
dupes = combined_gdf.duplicated(subset=["Xr", "Yr"])
print("Any duplicates left? ", dupes.any())

# 11. Convert CONSTRUCTDATE to datetime and compute line_age_yr
combined_gdf['CONSTRUCTDATE'] = pd.to_datetime(combined_gdf['CONSTRUCTDATE'], errors='coerce')
today = pd.Timestamp.now()
combined_gdf['line_age_yr'] = (today - combined_gdf['CONSTRUCTDATE']).dt.days / 365.25
print(combined_gdf[['CONSTRUCTDATE', 'line_age_yr']].head())

# 12. Normalize operator names and build mapping table
mapping = {
    'KINDER MORGAN CO2 CO LP': 'KINDER MORGAN CO2 CO LLC',
    'BEEMAN OIL & GAS INC': 'BEEMAN OIL & GAS LLC',
}
combined_gdf['Operator'] = combined_gdf['Operator'].replace(mapping)
combined_gdf.rename(columns={'OPERATOR_NUM': 'operator_number', 'Operator': 'operator_name'}, inplace=True)

combined_gdf_operator_mapping = combined_gdf[['operator_number', 'operator_name']].drop_duplicates().reset_index(drop=True)
print(combined_gdf_operator_mapping)

# 13. Rename spill_* columns (except spill_match_distance_m)
rename_dict = {
    col: col[len("spill_"):]
    for col in combined_gdf.columns
    if col.startswith("spill_") and col != "spill_match_distance_m"
}
print("Will rename these columns:")
for old, new in rename_dict.items():
    print(f"  {old} → {new}")
combined_gdf = combined_gdf.rename(columns=rename_dict)
print("\nNew column list:")
print(combined_gdf.columns.tolist())

# 14. Drop unwanted columns (if they exist)
columns_to_remove = [
    'trkg_num', 'Operator Name', 'facility_type', 'Spill_Desc', 'Spill Type',
    'Root Cause', 'Preventative Measure', 'Detailed Root Cause Type',
    'Long', 'Lat', 'facility_status', 'Metallic?', 'ACTIONDESCRIPTION',
    'BEDDINGMATERIAL', 'COMPANY_NAME', 'ENDLAT', 'ENDLONG',
    'ENTIRELINEREMOVED', 'PIPEMATERIAL', 'RECEIVE_DATE', 'STARTLAT',
    'STARTLOCATIONID', 'STARTLONG', 'TYPEOFFLUIDTRANS', 'operator_name',
    'SHAPE_Length', 'matched_crudeoil_idx', 'flowline_match_distance_m',
    'spill_match_distance_m'
]
cols_to_drop = [c for c in columns_to_remove if c in combined_gdf.columns]
print(f"Dropping these existing columns: {cols_to_drop}")
combined_gdf = combined_gdf.drop(columns=cols_to_drop)

# 15. Reorder DataFrame columns (verify existence first)
new_order = [
    'unique_id', 'operator_number', 'FLOWLINEID', 'LOCATION_ID',
    'Status', 'FLOWLINEACTION', 'LOCATIONTYPE', 'Fluid',
    'Material', 'Diam_in', 'Length_ft', 'MAXOPPRESSURE', 'line_age_yr',
    'CONSTRUCTDATE', 'incident_date', 'Root Cause Type', 'coords',
    'geometry', 'risk'
]
missing_cols = [c for c in new_order if c not in combined_gdf.columns]
if missing_cols:
    raise KeyError(f"Missing columns before reordering: {missing_cols}")
combined_gdf = combined_gdf[new_order]
print("After reordering, columns are:", combined_gdf.columns.tolist())

# 16. Print unique values for certain columns
columns_to_check = ['Status', 'FLOWLINEACTION', 'LOCATIONTYPE', 'Fluid', 'Material']
for column in columns_to_check:
    if column in combined_gdf.columns:
        print(f"Unique values in {column}: {combined_gdf[column].unique().tolist()}")
    else:
        print(f"{column}: Column not found in DataFrame.")

status_mapping = {
    'Active': 'Active', 'ACTIVE': 'Active', 'Actove': 'Active', 'Avtive': 'Active',
    'Actve': 'Active', 'active': 'Active', 'Actuve': 'Active', 'Flowline': 'Active',
    'Out of Service': 'Out of Service', 'OOS': 'Out of Service', 'OutofService': 'Out of Service',
    'Out-of-Service': 'Out of Service', 'Out Of Service': 'Out of Service',
    'Out of service': 'Out of Service', 'OUT OF SERVICE': 'Out of Service',
    'OOSLAT': 'Out of Service',
    'Abandoned': 'Abandoned', 'abandoned': 'Abandoned', 'Abandoned in Place': 'Abandoned',
    'ABANDONED': 'Abandoned', 'Abandon': 'Abandoned', 'Abadnon': 'Abandoned',
    'TA': 'Abandoned', 'Abandonded': 'Abandoned', 'Plugged': 'Abandoned', 'AB': 'Abandoned',
    'Inactive': 'Inactive', 'InActive': 'Inactive', 'INACTIVE': 'Inactive', 'IDLE': 'Inactive',
    'Pending Analysis': 'Pending Analysis', 'PA': 'Pending Analysis',
    'ABiP': 'Pending Analysis', 'ABIP': 'Pending Analysis', 'AIP': 'Pending Analysis',
    'Shut In': 'Shut In', 'shut in': 'Shut In', 'Shut in': 'Shut In', 'SI': 'Shut In',
    'Pre-Abandonment': 'Pre-Abandonment', 'Pre Abandonment': 'Pre-Abandonment',
    'PreAbandonment': 'Pre-Abandonment', 'Pre-Abandoned': 'Pre-Abandonment',
    'Removed': 'Removed', 'REMOVED': 'Removed', 'REMOVED DUPLICATE': 'Removed',
    'Partial Removed see comment': 'Removed',
    'New Construction': 'New Construction',
    'PreCommission': 'Pre-Commissioned', 'Pre-commissioned': 'Pre-Commissioned',
    'Pre-Commissioned': 'Pre-Commissioned',
    'Flushed and capped': 'Decommissioned',
    'Flow filled and capped': 'Decommissioned',
    'Discontinued': 'Decommissioned',
    'Injection': 'Injection', 'Injecting': 'Injection',
    'Future': 'Planned', 'Planned': 'Planned',
    'Nitrogen': 'DROP_ROW',
    'Air Line': 'DROP_ROW',
    '3rd Party Line': 'DROP_ROW',
    'Turned into Air Line': 'DROP_ROW',
    'Unknown': None, 'Status': None, 'unk': None, '0': None
}
combined_gdf['Status'] = combined_gdf['Status'].replace(status_mapping)
combined_gdf = combined_gdf[combined_gdf['Status'] != 'DROP_ROW']
print(combined_gdf['Status'].unique())

flowlineaction_mapping = {
    'Out of Service': 'Out of Service', 'Removed From Service': 'Out of Service',
    'Pre-Abandonment Notice': 'Pre-Abandonment Notice',
    'Abandonment Verification': 'Abandonment',
    'Realignment': 'Realignment', 'Registration': 'Registration',
    'Abandonment': 'Abandonment'
}
combined_gdf['FLOWLINEACTION'] = combined_gdf['FLOWLINEACTION'].replace(flowlineaction_mapping)
print(combined_gdf['FLOWLINEACTION'].unique())

locationtype_mapping = {
    'Production Facilities': 'Production Facilities', 'Well Site': 'Well Site',
    'Manifold': 'Manifold', 'Compressor Station': 'Compressor Station',
    'Gathering Line': 'Gathering Line', 'Crude Oil Transfer Line': 'Crude Oil Transfer Line',
    'Produced Water Transfer System': 'Produced Water Transfer System'
}
combined_gdf['LOCATIONTYPE'] = combined_gdf['LOCATIONTYPE'].replace(locationtype_mapping)
print(combined_gdf['LOCATIONTYPE'].unique())

# 18. Normalize "Fluid"
fluid_mapping = {
    'Natual Gas': 'Natural Gas', 'Natural Gas Production': 'Natural Gas',
    'Natural Gas Lift': 'Natural Gas', 'Natuarl Gas': 'Natural Gas',
    'Natural Gas High Pressure': 'Natural Gas', 'Natural Gas Supply': 'Natural Gas',
    'Natural Gas Sales': 'Natural Gas', 'Supply Gas': 'Natural Gas',
    'Fuel Gas': 'Natural Gas',
    'Oil': 'Crude Oil', 'Crude Oil': 'Crude Oil',
    'Crude Oil Unprocessed': 'Crude Oil',
    'Crude Oil Emulsion': 'Crude Oil Emulsion',
    'Crude Oil Emmulsion, Water And Oil': 'Crude Oil Emulsion',
    'Crude Oill Emulsion': 'Crude Oil Emulsion',
    'Crude Oil And Water Emulsion': 'Crude Oil Emulsion',
    'Oil Water Emulsion': 'Crude Oil Emulsion',
    'Oil/Water': 'Crude Oil Emulsion', 'Oil Water': 'Crude Oil Emulsion',
    'Oil And Water': 'Crude Oil Emulsion', 'Emulsion': 'Crude Oil Emulsion',
    'Crude Oil / Produced Water': 'Crude Oil Emulsion',
    'Produced Water': 'Produced Water', 'Water': 'Produced Water',
    'Saltwater': 'Produced Water', 'Injection Produced Water': 'Produced Water',
    'Produced  Water': 'Produced Water', 'Produced/Waste Water': 'Produced Water',
    'Waste Water/Produced Water/Formation Water': 'Produced Water',
    'Fresh Water': 'Produced Water', '(Other) Treated Produced Water': 'Produced Water',
    'Co2': 'Co2/Produced Water', 'C02/Prod Water': 'Co2/Produced Water',
    'Co2/Prod Water': 'Co2/Produced Water', 'Co2Produced Water': 'Co2/Produced Water',
    'Co2/Produced Wtaer': 'Co2/Produced Water', 'Co2 Produced Water': 'Co2/Produced Water',
    'Co2/Produced': 'Co2/Produced Water',
    '3 Phase': 'Multiphase', 'Multiphase': 'Multiphase', 'Multi-Phase': 'Multiphase',
    'Mulitphase': 'Multiphase', 'Multi Phase': 'Multiphase', 'Mulit Phase': 'Multiphase',
    'Mutliphase': 'Multiphase', 'Oil, Water, Gas': 'Multiphase', 'Oil-Gas-Water': 'Multiphase',
    'Gas, Oil And Water': 'Full Well Stream', 'Oil /Water/Gas': 'Full Well Stream',
    'Oil/Gas/Water': 'Full Well Stream', 'Oil, Gas, Water': 'Full Well Stream',
    'Gas,  Oil And Water': 'Full Well Stream', 'Natural Gas/Condensate/Produced Water': 'Full Well Stream',
    'Condensate': 'Condensate',
    'Poly': 'Polymer Fluids', 'Polymer fluids': 'Polymer Fluids',
    'Liquid': 'Other', 'Liquids (Wtr/Cond)': 'Other',
    'Unprocessed Production Fluids': 'Other', 'Production Fluids': 'Other',
    'Produced Fluids': 'Other', 'Other': 'Other',
    'Other (Natural Gas Liquid)': 'Other', 'Other- Nitrogen/ Natural Gas': 'Other',
    'Service Line': 'Other', 'Air': 'Other', 'Cox V2': 'Other',
    'Steel': 'Other', 'Fluid': 'Other', 'Vent': 'Other',
    'Unk': None, 'Unknown': None
}
combined_gdf['Fluid'] = combined_gdf['Fluid'].str.strip().str.title().replace(fluid_mapping)
print(sorted(combined_gdf['Fluid'].dropna().unique()))

# 19. Normalize "Material"
material_mapping = {
    'Carbon Steel': 'Carbon Steel', 'Carbonsteel': 'Carbon Steel', 'Carbon  Steel': 'Carbon Steel',
    'Carbon Steel Sch 80': 'Carbon Steel', 'Carbon Steel - Hdpe': 'Carbon Steel/HDPE',
    'Carbon Steel-Hdpe': 'Carbon Steel/HDPE', 'Carbon Steel And Hdpe': 'Carbon Steel/HDPE',
    'Carbon Steel/Hdpe': 'Carbon Steel/HDPE', 'Carbon Steel/HDPE': 'Carbon Steel/HDPE',
    'Hdpe/Steel': 'Carbon Steel/HDPE', 'HDPE/Steel': 'Carbon Steel/HDPE',
    'Hdpe Lined Steel': 'Carbon Steel/HDPE', 'Hdpe/Steel, Flexsteel': 'Carbon Steel/HDPE/Flexsteel',
    'Carbon Steel, Hdpe,Stainless Steel': 'Carbon Steel/HDPE/Stainless Steel',
    'Carbon Steel, Hdpe, Stainless Steel': 'Carbon Steel/HDPE/Stainless Steel',
    'Carbon Steel/Stainless Steel/Hdpe': 'Carbon Steel/HDPE/Stainless Steel',
    'Carbon Steel/Hdpe/Stainless': 'Carbon Steel/HDPE/Stainless Steel',
    'Carbon Steel/Hdpe/Stainless Steel': 'Carbon Steel/HDPE/Stainless Steel',
    'Carbon Steel/Stainless/Hdpe': 'Carbon Steel/HDPE/Stainless Steel',
    'Stainless/Carbon Steel/Hdpe': 'Carbon Steel/HDPE/Stainless Steel',
    'Stainless/ Carbon Steel/Hdpe': 'Carbon Steel/HDPE/Stainless Steel',
    'Stainless Steel/Carbon Steel/Hdpe': 'Carbon Steel/HDPE/Stainless Steel',
    'Stainless/Carbonsteel/Hdpe': 'Carbon Steel/HDPE/Stainless Steel',
    'Satinless/Carbon Steel/Hdpe': 'Carbon Steel/HDPE/Stainless Steel',
    'Stainless/Carbon Steel/Hspe': 'Carbon Steel/HDPE/Stainless Steel',
    'Stainless/Carbon Steel/ Hdpe': 'Carbon Steel/HDPE/Stainless Steel',
    'Carbon Steel/Fiberglass': 'Fiberglass/Carbon Steel',
    'Carbon Steel Mixed With Fiberglass': 'Fiberglass/Carbon Steel',
    'Fiberglass': 'Fiberglass', 'Fibergalss': 'Fiberglass', 'Fiberspar': 'Fiberglass',
    'Fiber Glass': 'Fiberglass', 'Flexspar': 'Fiberglass', 'Fiberoptic': 'Fiberglass',
    'Fiberglass & Fiberspar': 'Fiberglass', 'Fiberglass And Carbon Steel': 'Fiberglass/Carbon Steel',
    'Fiberglass Sleaved W/ Hdpe': 'Fiberglass/HDPE', 'Fiberglass And Hdpe': 'Fiberglass/HDPE',
    'Hdpe/Fiberglass': 'Fiberglass/HDPE',
    'Steel': 'Steel', 'Stainless': 'Steel', 'Stainless Steel': 'Steel',
    'Coated Steel': 'Steel', 'Lined Steel': 'Steel', 'Flexpipe': 'Steel',
    'Flex Steel': 'Steel', 'Flexsteel': 'Steel',
    'Other (Flexsteel)': 'Steel', 'Other (Flex Steel)': 'Steel', 'Other (Flex Pipe)': 'Steel',
    'Hdpe': 'HDPE', 'Hdpe Poly': 'HDPE', 'Hdpi Poly': 'HDPE',
    'Hdpe Line Sdr 7': 'HDPE', 'Hdpe Sdr7': 'HDPE', 'Hdpe Sdr11': 'HDPE',
    'Composite Hdpe': 'HDPE', 'High-Density Polyethylene (Hdpe)': 'HDPE',
    'Poly': 'Polycarbonate', 'Polyline': 'Polycarbonate', 'Polycarbonate': 'Polycarbonate',
    'Other (Poly)': 'Polycarbonate', 'Sdr 7 Poly': 'Polycarbonate',
    'Poly & Steel': 'Polycarbonate/Steel', 'Steel/Poly': 'Polycarbonate/Steel',
    'Poly/Steel': 'Polycarbonate/Steel', 'Polycarbonate/Steel': 'Polycarbonate/Steel',
    'Poly Sleeved Steel': 'Polycarbonate/Steel', 'Hdpe Poly Sdr 11': 'Polycarbonate/HDPE',
    'Sdr7 Polyethelyne': 'Polyethylene', 'Sdr 11 Poly Pipe': 'Polyethylene',
    'Sdr 11 Poly': 'Polyethylene', 'Poly Pipe': 'Polyethylene',
    'Sdr_Poly': 'Polyethylene', 'Sdr-11': 'Polyethylene', 'Sdr-11 Poly': 'Polyethylene',
    'Poly Sdr 7': 'Polypropylene', 'Poly Sdr-7': 'Polypropylene',
    'Pvc': 'PVC',
    'Duplex': 'Duplex',
    'Zaplock': 'Other',
    'Plastic': 'Other', 'Polypipe': 'Other', 'Core Linepipe': 'Other', 'Flex Pipe': 'Other',
    'Shawcor Fp150': 'Other', 'Zapock': 'Other', 'Other (Hdpe And Tubing)': 'Other',
    'Other (Please Specify)': 'Other', 'Other (Stainless Steel)': 'Other',
    'Fplp': 'Other', 'Oil': 'Other', 'Gas': 'Other', 'Co2/Produced Water': 'Other',
    'Other': 'Other',
    'Unknown': None, 'Other (Unknown)': None, 'Unk': None, '0': None,
    'Material': None, 'Sdr': None, 'Flowline': None
}
combined_gdf['Material'] = combined_gdf['Material'].str.strip().str.title().replace(material_mapping)
print(combined_gdf['Material'].unique())

# 20. Normalize "Root Cause Type"
root_cause_mapping = {
    'Corrosion': 'Corrosion', 'Unknown': 'Unknown', 'Incorrect Operation': 'Incorrect Operation',
    'Equipment Failure': 'Equipment Failure', 'Equipment failure': 'Equipment Failure',
    'Other Outside Force Damage': 'Other Outside Force Damage', 'Natural Force Damage': 'Natural Force Damage',
    'Pipe, Weld, or Joint Failure': 'Pipe, Weld, or Joint Failure',
    'Pipe, Weld Joint Failure': 'Pipe, Weld, or Joint Failure', 'Excavation Damage': 'Excavation Damage',
    'Other Outside Force': 'Other Outside Force Damage',  # unify name
    'Pipe, Weld, Joint Failure': 'Pipe, Weld, or Joint Failure'
}
combined_gdf['Root Cause Type'] = combined_gdf['Root Cause Type'].replace(root_cause_mapping)
print(combined_gdf['Root Cause Type'].unique())

# 21. Rename columns to lowercase and snake_case using a mapping rather than full assignment
old_to_new = {
    'FLOWLINEID': 'flowline_id',
    'LOCATION_ID': 'location_id',
    'Status': 'status',
    'FLOWLINEACTION': 'flowline_action',
    'LOCATIONTYPE': 'location_type',
    'Fluid': 'fluid',
    'Material': 'material',
    'Diam_in': 'diameter_in',
    'Length_ft': 'length_ft',
    'MAXOPPRESSURE': 'max_operating_pressure',
    'CONSTRUCTDATE': 'construct_date',
    'incident_date': 'spill_date',
    'Root Cause Type': 'root_cause',
    'coords': 'coords',  # keep same name
    'Xr': 'Xr',
    'Yr': 'Yr',
    'line_age_yr': 'line_age_yr',
    'geometry': 'geometry',
    'risk': 'risk',
    'unique_id': 'unique_id',
    'operator_number': 'operator_number',
}
# Only rename existing columns
rename_existing = {old: new for old, new in old_to_new.items() if old in combined_gdf.columns}
combined_gdf = combined_gdf.rename(columns=rename_existing)
print("Columns after renaming:", combined_gdf.columns.tolist())

# 22. Drop rows with missing critical values (max_operating_pressure AND risk=0)
if 'max_operating_pressure' in combined_gdf.columns:
    combined_gdf = combined_gdf[~((combined_gdf['max_operating_pressure'].isna()) & (combined_gdf['risk'] == 0))]

# 23. KNN Imputation for numeric columns
cols_for_impute = [c for c in ['max_operating_pressure', 'diameter_in', 'length_ft', 'line_age_yr', 'material', 'fluid'] if c in combined_gdf.columns]
df_for_imputation = combined_gdf[cols_for_impute].copy()

le = LabelEncoder()
if 'material' in df_for_imputation.columns:
    df_for_imputation['material_encoded'] = le.fit_transform(df_for_imputation['material'].fillna(''))
if 'fluid' in df_for_imputation.columns:
    df_for_imputation['fluid_encoded'] = le.fit_transform(df_for_imputation['fluid'].fillna(''))
df_for_imputation = df_for_imputation.drop(columns=[c for c in ['material','fluid'] if c in df_for_imputation.columns])

imputer = KNNImputer(n_neighbors=5)
df_imputed = imputer.fit_transform(df_for_imputation)

idx = 0
for c in ['max_operating_pressure', 'diameter_in', 'length_ft', 'line_age_yr']:
    if c in cols_for_impute:
        combined_gdf[c] = df_imputed[:, idx]
        idx += 1

# 24. Count NaNs and fill for risk=1 rows
def count_only_nan(series):
    return series.apply(lambda x: 1 if isinstance(x, float) and np.isnan(x) else 0).sum()

na_columns = combined_gdf.apply(count_only_nan)
columns_with_only_nan = na_columns[na_columns > 0]
print("Columns with NaN values (excluding None) and their counts:")
print(columns_with_only_nan)

columns_with_na = combined_gdf.columns[combined_gdf.isna().any()]
na_with_risk_1 = {}
for column in columns_with_na:
    na_with_risk_1[column] = combined_gdf[(combined_gdf[column].isna()) & (combined_gdf['risk'] == 1)].shape[0]
print("Number of NaNs with risk of 1 in each column:")
print(na_with_risk_1)

to_fill = [c for c in ['status','flowline_action','fluid','material'] if c in combined_gdf.columns]
for col in to_fill:
    combined_gdf.loc[(combined_gdf['risk'] == 1) & combined_gdf[col].isna(), col] = 'None'

# Drop any remaining rows missing flowline_id, location_id, or diameter_in
drop_na_cols = [c for c in ['flowline_id','location_id','diameter_in'] if c in combined_gdf.columns]
combined_gdf = combined_gdf.dropna(subset=drop_na_cols)

# 25. Convert numeric columns to integer where appropriate
for col in ['diameter_in', 'length_ft', 'max_operating_pressure', 'line_age_yr']:
    if col in combined_gdf.columns:
        combined_gdf[col] = combined_gdf[col].astype(int)

# Print class counts
total_rows    = combined_gdf.shape[0]
risk_1_count  = combined_gdf[combined_gdf['risk'] == 1].shape[0]
risk_0_count  = combined_gdf[combined_gdf['risk'] == 0].shape[0]
print(f"Total number of rows: {total_rows}")
print(f"Total number of rows with risk = 1: {risk_1_count}")
print(f"Total number of rows with risk = 0: {risk_0_count}")

# 26. Write cleaned GeoJSON and operator mapping CSV
combined_gdf.to_file("cleaned_gdf.geojson", driver='GeoJSON')
combined_gdf_operator_mapping.to_csv('operator_mapping.csv', index=False)

# 27. Load additional datasets and join population density and elevation
pop_density = gpd.read_file('Population_Density_(Census_Tracts)').to_crs(combined_gdf.crs)
flowlines   = gpd.read_file('cleaned_gdf.geojson')

print('Summary of Census Tract Data:')
print(pop_density.info())
print('\nFirst few rows of the data:')
print(pop_density.head())

joined = gpd.sjoin(flowlines, pop_density, how='left', predicate='within')
flowlines['avg_population'] = joined['Populati_1']

import rasterio
from rasterio import features
from pyproj import CRS

dem = rasterio.open('output_USGS30m.tif')
dem_crs = CRS(dem.crs)
if flowlines.crs != dem_crs:
    flowlines = flowlines.to_crs(dem_crs)

def get_elevation(point, dem):
    try:
        val = list(dem.sample([(point.x, point.y)]))[0][0]
        if dem.nodata is not None and val == dem.nodata:
            return None
        return val
    except:
        return None

flowlines['avg_elevation'] = flowlines.geometry.apply(lambda pt: get_elevation(pt, dem))
flowlines = flowlines.drop(columns=['index_right'], errors='ignore')

# 28. Reorder final dataset columns
final_order = [
    'unique_id', 'operator_number', 'flowline_id', 'location_id', 'status',
    'flowline_action', 'location_type', 'fluid', 'material', 'diameter_in',
    'length_ft', 'max_operating_pressure', 'avg_population', 'avg_elevation',
    'line_age_yr', 'construct_date', 'spill_date', 'root_cause', 'risk',
    'coords', 'geometry'
]
final_order = [c for c in final_order if c in flowlines.columns]
flowlines = flowlines.reindex(columns=final_order)

# 29. Save final point-based GeoJSON
flowlines.to_file("final_dataset.geojson", driver='GeoJSON')
print("✅ Saved final_dataset.geojson with population and elevation for point-based flowlines.")
