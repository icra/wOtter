import geopandas
import numpy
import os

# Load data
directory = os.path.join(os.path.dirname(os.getcwd()), 'data')

# input
AGG_WWTP_location = rivers_shapefile = os.path.join(directory, "AGG_WWTP.shp")
moved_discharge_location = os.path.join(directory, "Raw data/move_discharge.csv")

AGG_WWTP_df = geopandas.read_file(AGG_WWTP_location)


# output
AGG_WWTP_adapted = os.path.join(directory, "AGG_WWTP_adapted.shp")


# location changes:
moved_pixels = geopandas.read_file(moved_discharge_location)
destinations = moved_pixels.iloc[:, 0]
key_values = moved_pixels.iloc[:, 1]

for i in range(len(destinations)):
    destination = destinations[i].split(',')
    location = AGG_WWTP_df['dcpLongitu'] == float(key_values[i])
    AGG_WWTP_df['dcpLatitud'][location] = destination[0]
    AGG_WWTP_df['dcpLongitu'][location] = destination[1]

AGG_WWTP_df.to_file(AGG_WWTP_adapted)



