# # import geopandas as gpd;
# from pyproj import Proj, transform

# gdf = gpd.GeoDataFrame(gpd.read_file("COGCC_Form44_Crude_Oil_Produced_Water_Transfer_Flowlines_Approved_CONFIDENTIAL.gdb"))

# gdf.to_csv("crudeOilProduced.csv")

# #print("Done")
import pandas as pd
from shapely.wkt import loads
from pyproj import Transformer
import csv

# Read the CSV file into a DataFrame
df = pd.read_csv("crudeOilProduced_initial.csv")

# Define the source and destination coordinate systems using Transformer from PyProj
transformer = Transformer.from_crs("epsg:26913", "epsg:4326", always_xy=True)

# Open a CSV file for writing the transformed coordinates
with open('transformed_coordinates.csv', 'w', newline='') as csvfile:
    csvwriter = csv.writer(csvfile)
    # Write header
    csvwriter.writerow(['Line', 'Longitude', 'Latitude'])

    # Iterate over each row in the DataFrame
    for index, row in df.iterrows():
        # Parse the WKT to a Shapely geometry
        multiline = loads(row['geometry'])

        # Check if the geometry is indeed a MultiLineString
        if multiline.geom_type == 'MultiLineString':
            for line in multiline.geoms:
                for point in line.coords:
                    # Perform the coordinate transformation
                    lon, lat = transformer.transform(point[0], point[1])
                    # Write each point's coordinates to the CSV
                    csvwriter.writerow([index + 1, lon, lat])
        else:
            print(f"Geometry at index {index} is not a MultiLineString")

print("Coordinate transformation and CSV export done")

