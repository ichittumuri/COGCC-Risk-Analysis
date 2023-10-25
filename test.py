import os
import numpy as np
import pandas as pd
import geopandas as gpd

data_path = "ECMC_CSM_Flowline_Data_Access/COGCC_Form44_Off_Location_Flowlines_Approved_CONFIDENTIAL.gdb"

print(f"{os.path.exists(data_path) = }")

df = gpd.read_file(data_path)

print(df)
print(df.head())



