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
line_data = gpd.read_file("full_length_flowlines.geojson").to_crs(epsg=4326)
line_data["risk"] = 0  # start all flowlines with risk = 0

# 11. Convert CONSTRUCTDATE to datetime and compute line_age_yr
line_data['CONSTRUCTDATE'] = pd.to_datetime(line_data['CONSTRUCTDATE'], errors='coerce')
today = pd.Timestamp.now()
line_data['line_age_yr'] = (today - line_data['CONSTRUCTDATE']).dt.days / 365.25
print(line_data[['CONSTRUCTDATE', 'line_age_yr']].head())

# 12. Normalize operator names and build mapping table
mapping = {
    'KINDER MORGAN CO2 CO LP': 'KINDER MORGAN CO2 CO LLC',
    'BEEMAN OIL & GAS INC': 'BEEMAN OIL & GAS LLC',
}
line_data['Operator'] = line_data['Operator'].replace(mapping)
line_data.rename(columns={'OPERATOR_NUM': 'operator_number', 'Operator': 'operator_name'}, inplace=True)

# 14. Drop unwanted columns (if they exist)
columns_to_remove = [
    'trkg_num', 'Operator Name', 'facility_type', 'Spill_Desc', 'Spill Type',
    'Root Cause', 'Preventative Measure', 'Detailed Root Cause Type',
    'Long', 'Lat', 'facility_status', 'Metallic?', 'ACTIONDESCRIPTION',
    'BEDDINGMATERIAL', 'COMPANY_NAME', 'ENDLAT', 'ENDLONG',
    'ENTIRELINEREMOVED', 'PIPEMATERIAL', 'RECEIVE_DATE', 'STARTLAT',
    'STARTLOCATIONID', 'STARTLONG', 'TYPEOFFLUIDTRANS', 'operator_name',
    'SHAPE_Length', 'matched_crudeoil_idx'
]
cols_to_drop = [c for c in columns_to_remove if c in line_data.columns]
print(f"Dropping these existing columns: {cols_to_drop}")
line_data = line_data.drop(columns=cols_to_drop)

# 15. Reorder DataFrame columns (verify existence first)
new_order = [
    'unique_id', 'operator_number', 'FLOWLINEID', 'LOCATION_ID',
    'Status', 'FLOWLINEACTION', 'LOCATIONTYPE', 'Fluid',
    'Material', 'Diam_in', 'Length_ft', 'MAXOPPRESSURE', 'line_age_yr',
    'CONSTRUCTDATE', 'geometry', 'risk'
]
missing_cols = [c for c in new_order if c not in line_data.columns]
if missing_cols:
    raise KeyError(f"Missing columns before reordering: {missing_cols}")
line_data = line_data[new_order]
print("After reordering, columns are:", line_data.columns.tolist())

# 16. Print unique values for certain columns
columns_to_check = ['Status', 'FLOWLINEACTION', 'LOCATIONTYPE', 'Fluid', 'Material']
for column in columns_to_check:
    if column in line_data.columns:
        print(f"Unique values in {column}: {line_data[column].unique().tolist()}")
    else:
        print(f"{column}: Column not found in DataFrame.")

# 17. Normalize categorical values
status_mapping = {
    'Active': 'Active', 'ACTIVE': 'Active', 'Actove': 'Active', 'Avtive': 'Active',
    'Actve': 'Active', 'active': 'Active', 'Out of Service': 'Out of Service',
    'OOS': 'Out of Service', 'OutofService': 'Out of Service',
    'Out-of-Service': 'Out of Service', 'Out Of Service': 'Out of Service',
    'Out of service': 'Out of Service', 'Abandoned': 'Abandoned', 'abandoned': 'Abandoned',
    'Abandoned in Place': 'Abandoned', 'ABANDONED': 'Abandoned', 'Abandon': 'Abandoned',
    'Abadnon': 'Abandoned', 'TA': 'Abandoned', 'Inactive': 'Inactive',
    'InActive': 'Inactive', 'INACTIVE': 'Inactive', 'PA': 'Pending Analysis',
    'ABiP': 'Pending Analysis', 'Shut in': 'Shut In', 'shut in': 'Shut In',
    'SI': 'Shut In', 'Status': 'Unknown', 'Future': 'Future',
    'REMOVED': 'Removed', 'Pre Abandonment': 'Pre-Abandonment',
    'PreAbandonment': 'Pre-Abandonment'
}
line_data['Status'] = line_data['Status'].replace(status_mapping)
print(line_data['Status'].unique())

flowlineaction_mapping = {
    'Out of Service': 'Out of Service', 'Removed From Service': 'Out of Service',
    'Pre-Abandonment Notice': 'Pre-Abandonment Notice',
    'Abandonment Verification': 'Abandonment',
    'Realignment': 'Realignment', 'Registration': 'Registration',
    'Abandonment': 'Abandonment'
}
line_data['FLOWLINEACTION'] = line_data['FLOWLINEACTION'].replace(flowlineaction_mapping)
print(line_data['FLOWLINEACTION'].unique())

locationtype_mapping = {
    'Production Facilities': 'Production Facilities', 'Well Site': 'Well Site',
    'Manifold': 'Manifold', 'Compressor Station': 'Compressor Station',
    'Gathering Line': 'Gathering Line', 'Crude Oil Transfer Line': 'Crude Oil Transfer Line',
    'Produced Water Transfer System': 'Produced Water Transfer System'
}
line_data['LOCATIONTYPE'] = line_data['LOCATIONTYPE'].replace(locationtype_mapping)
print(line_data['LOCATIONTYPE'].unique())

# 18. Normalize "Fluid"
line_data['Fluid'] = line_data['Fluid'].str.strip().str.title().replace({
    'Natual Gas': 'Natural Gas', 'Natural Gas Production': 'Natural Gas',
    'Co2': 'Co2/Produced Water', 'C02/Prod Water': 'Co2/Produced Water',
    'Co2/Prod Water': 'Co2/Produced Water', 'Co2Produced Water': 'Co2/Produced Water',
    'Co2/Produced Wtaer': 'Co2/Produced Water', 'Gas': 'Natural Gas',
    'Gas, Oil And Water': 'Full Well Stream', 'Oil': 'Crude Oil',
    'Crude Oil': 'Crude Oil', 'Crude Oil Emulsion': 'Crude Oil Emulsion',
    'Emulsion': 'Crude Oil Emulsion', 'Crude Oil Emmulsion, Water And Oil': 'Crude Oil Emulsion',
    'Crude Oil And Water Emulsion': 'Crude Oil Emulsion', 'Oil Water Emulsion': 'Crude Oil Emulsion',
    'Oil/Water': 'Crude Oil Emulsion', 'Oil Water': 'Crude Oil Emulsion',
    'Oil And Water': 'Crude Oil Emulsion', 'Oil /Water/Gas': 'Full Well Stream',
    'Oil/Gas/Water': 'Full Well Stream', 'Oil, Gas, Water': 'Full Well Stream',
    '3 Phase': 'Multiphase', 'Multiphase': 'Multiphase', 'Multi-Phase': 'Multiphase',
    'Mulitphase': 'Multiphase', 'Multi Phase': 'Multiphase', 'Mulit Phase': 'Multiphase',
    'Multi-Phase\xa0': 'Multiphase', 'Injection Produced Water': 'Produced Water',
    'Produced Water': 'Produced Water', 'Water': 'Produced Water', 'Saltwater': 'Produced Water',
    'Condensate': 'Condensate', 'Liquid': 'Other', 'Liquids (Wtr/Cond)': 'Other',
    'Unprocessed Production Fluids': 'Other', 'Production Fluids': 'Other',
    'Produced Fluids': 'Other', 'Full Well Stream': 'Full Well Stream',
    'Other': 'Other', 'Gas,  Oil And Water': 'Full Well Stream',
    'Natural Gas Lift': 'Natural Gas', 'Natuarl Gas': 'Natural Gas',
    'Natural Gas High Pressure': 'Natural Gas', 'Natural Gas Supply': 'Natural Gas',
    'Crude Oill Emulsion': 'Crude Oil Emulsion', 'Unk': 'Unknown', 'Poly': 'Polymer fluids'
})
print(line_data['Fluid'].unique())

# 19. Normalize "Material"
line_data['Material'] = line_data['Material'].str.strip().str.title().replace({
    'Fiberglass': 'Fiberglass', 'Fibergalss': 'Fiberglass', 'Fiberspar': 'Fiberglass',
    'Fiber Glass': 'Fiberglass', 'Carbon Steel': 'Carbon Steel', 'Carbonsteel': 'Carbon Steel',
    'Carbon Steel Sch 80': 'Carbon Steel', 'Carbon Steel - Hdpe': 'Carbon Steel/HDPE',
    'Carbon Steel, Hdpe,Stainless Steel': 'Carbon Steel/HDPE/Stainless Steel',
    'Carbon Steel, Hdpe, Stainless Steel': 'Carbon Steel/HDPE/Stainless Steel',
    'Carbon Steel/Stainless Steel/Hdpe': 'Carbon Steel/HDPE/Stainless Steel',
    'Carbon Steel/Hdpe/Stainless': 'Carbon Steel/HDPE/Stainless Steel',
    'Carbon Steel/Hdpe': 'Carbon Steel/HDPE', 'Satinless/Carbon Steel/Hdpe': 'Carbon Steel/HDPE/Stainless Steel',
    'Carbon Steel/Stainless/Hdpe': 'Carbon Steel/HDPE/Stainless Steel', 'Steel': 'Steel',
    'Lined Steel': 'Steel', 'Coated Steel': 'Steel', 'Flexsteel': 'Steel',
    'Flexpipe': 'Steel', 'Fiber Glass And Carbon Steel': 'Fiberglass/Carbon Steel',
    'Fiberglass And Hdpe': 'Fiberglass/HDPE', 'Hdpe': 'HDPE', 'Hdpe Poly': 'HDPE',
    'Composite Hdpe': 'HDPE', 'Hdpe/Steel': 'HDPE/Steel', 'Hdpe Lined Steel': 'HDPE/Steel',
    'Hdpe/Steel, Flexsteel': 'HDPE/Steel', 'Poly': 'Polycarbonate', 'Polyline': 'Polycarbonate',
    'Poly & Steel': 'Polycarbonate/Steel', 'Steel/Poly': 'Polycarbonate/Steel',
    'Poly/Steel': 'Polycarbonate/Steel', 'Polycarbonate': 'Polycarbonate',
    'Polycarbonate/Steel': 'Polycarbonate/Steel', 'Pvc': 'PVC', 'Flexspar': 'Fiberglass',
    'Stainless': 'Steel', 'Stainless/Carbon Steel/Hdpe': 'Carbon Steel/HDPE/Stainless Steel',
    'Carbon Steel/Hdpe/Stainless Steel': 'Carbon Steel/HDPE/Stainless Steel',
    'Unknown': 'Unknown', 'Other': 'Other', 'Other (Poly)': 'Polycarbonate',
    'Sdr7 Polyethelyne': 'Polyethylene', 'Sdr 11 Poly Pipe': 'Polyethylene',
    'Sdr 11 Poly': 'Polyethylene', 'Poly Pipe': 'Polyethylene', 'Sdr_Poly': 'Polyethylene',
    'Poly Sdr 7': 'Polypropylene', 'Poly Sdr-7': 'Polypropylene', 'Duplex': 'Duplex',
    'Fplp': 'Other', 'Flowline': 'Other', 'Flex Steel': 'Steel',
    'Other (Flex Steel)': 'Steel', 'Fiberglass And Carbon Steel': 'Carbon Steel/Fiberglass',
    'Stainless Steel': 'Steel', 'HDPE Lined Steel': 'HDPE/Steel', 'Fiberglass/Hdpe': 'Fiberglass/HDPE',
    'Unk': 'Unknown', 'Other (Unknown)': 'Unknown', 'Other': 'Unknown',
})
print(line_data['Material'].unique())

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
    'Xr': 'Xr',
    'Yr': 'Yr',
    'line_age_yr': 'line_age_yr',
    'geometry': 'geometry',
    'risk': 'risk',
    'unique_id': 'unique_id',
    'operator_number': 'operator_number',
}
# Only rename existing columns
rename_existing = {old: new for old, new in old_to_new.items() if old in line_data.columns}
line_data = line_data.rename(columns=rename_existing)
print("Columns after renaming:", line_data.columns.tolist())

# 22. Drop rows with missing critical values (max_operating_pressure AND risk=0)
if 'max_operating_pressure' in line_data.columns:
    line_data = line_data[~((line_data['max_operating_pressure'].isna()) & (line_data['risk'] == 0))]

# 23. KNN Imputation for numeric columns
cols_for_impute = [c for c in ['max_operating_pressure', 'diameter_in', 'length_ft', 'line_age_yr', 'material', 'fluid'] if c in line_data.columns]
df_for_imputation = line_data[cols_for_impute].copy()

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
        line_data[c] = df_imputed[:, idx]
        idx += 1

# 24. Count NaNs and fill for risk=1 rows
def count_only_nan(series):
    return series.apply(lambda x: 1 if isinstance(x, float) and np.isnan(x) else 0).sum()

na_columns = line_data.apply(count_only_nan)
columns_with_only_nan = na_columns[na_columns > 0]
print("Columns with NaN values (excluding None) and their counts:")
print(columns_with_only_nan)

columns_with_na = line_data.columns[line_data.isna().any()]
na_with_risk_1 = {}
for column in columns_with_na:
    na_with_risk_1[column] = line_data[(line_data[column].isna()) & (line_data['risk'] == 1)].shape[0]
print("Number of NaNs with risk of 1 in each column:")
print(na_with_risk_1)

to_fill = [c for c in ['status','flowline_action','fluid','material'] if c in line_data.columns]
for col in to_fill:
    line_data.loc[(line_data['risk'] == 1) & line_data[col].isna(), col] = 'None'

# Drop any remaining rows missing flowline_id, location_id, or diameter_in
drop_na_cols = [c for c in ['flowline_id','location_id','diameter_in'] if c in line_data.columns]
line_data = line_data.dropna(subset=drop_na_cols)

# 25. Convert numeric columns to integer where appropriate
for col in ['diameter_in', 'length_ft', 'max_operating_pressure', 'line_age_yr']:
    if col in line_data.columns:
        line_data[col] = line_data[col].astype(int)

# Print class counts
total_rows    = line_data.shape[0]
risk_1_count  = line_data[line_data['risk'] == 1].shape[0]
risk_0_count  = line_data[line_data['risk'] == 0].shape[0]
print(f"Total number of rows: {total_rows}")
print(f"Total number of rows with risk = 1: {risk_1_count}")
print(f"Total number of rows with risk = 0: {risk_0_count}")

line_data.to_file("cleaned_line_data.geojson", driver='GeoJSON')

# 27. Load additional datasets and join population density and elevation
pop_density = gpd.read_file('Population_Density_(Census_Tracts)').to_crs(line_data.crs)
flowlines   = gpd.read_file('cleaned_line_data.geojson')

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
    'line_age_yr', 'construct_date', 'risk', 'geometry'
]
final_order = [c for c in final_order if c in flowlines.columns]
flowlines = flowlines.reindex(columns=final_order)

# … after you’ve built & re-ordered flowlines …
print("Columns in flowlines:", flowlines.columns.tolist())
print(flowlines.head())

# 1) Show all dtypes so we can spot any hidden geometry column
print(flowlines.dtypes)

# 2) List all columns whose dtype is geometry
geom_cols = [col for col, dt in flowlines.dtypes.items() if dt.name == 'geometry']
print("Geometry-dtype columns:", geom_cols)

# 3) Show which one is currently the active geometry
print("Active geometry column:", flowlines.geometry.name)

# 4) If you see more than one in geom_cols, drop the extra(s)
for col in geom_cols:
    if col != 'geometry':
        flowlines = flowlines.drop(columns=[col])
        print(f"Dropped extra geometry column: {col!r}")

# 5) Re-set the one you want as the active geometry
flowlines = flowlines.set_geometry('geometry') 
print(flowlines.head())

# 6) Now write it out
flowlines.to_file("final_line_data.geojson", driver="GeoJSON")
print("✅ Saved final_line_data.geojson")