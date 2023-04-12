import pickle
from osgeo import gdal
import numpy
import pandas
import os

from src.library import graph_functions
from src.library import shapefile_raster_functions



directory = os.path.join(os.getcwd(), 'data')

# input locations
data_location = os.path.join(directory, "Wilkinson_z_normalized_adapted.csv")
reference_raster_location = os.path.join(directory, "reference_raster.tif")
graph_location = os.path.join(directory, "river_graph.pkl")

# temp locations
river_raster_location = "data/Rrivers_from_graph.tif"
temp_discharge_location = "disch.tif"

# output location
output_name = os.path.join(directory, "pollution_observed_adapted.csv")

# split coordinates into two columns
data = pandas.read_csv(data_location)
data['latitude'] = 0
data['longitude'] = 0
for i in range(len(data)):
    coordinates = data['Coord'].iloc[i]
    coordinates = coordinates.split(',')
    data['latitude'].iloc[i] = float(coordinates[0])
    data['longitude'].iloc[i] = float(coordinates[1])

# put the observation points into an indicator raster according to their location
reference_raster = gdal.Open(reference_raster_location)
discharge_array = numpy.zeros([reference_raster.RasterYSize, reference_raster.RasterXSize])
rows, columns = numpy.shape(discharge_array)
df_entries = len(data)
pixel_loc = numpy.zeros([df_entries]) - 10

for j in range(df_entries):
    i = df_entries - j - 1  # allows for dropping elements without skipping over numbers
    long, lat = data['longitude'].loc[i], data['latitude'].loc[i]
    lat_nr, long_nr = shapefile_raster_functions.give_pixel([lat, long], reference_raster)
    found = False
    if data['pollutant'][i] > -1:  # excludes nans
        if 0 <= lat_nr <= rows:
            if 0 <= long_nr <= columns:
                found = True
                discharge_array[lat_nr, long_nr] = 1
                pixel_loc[i] = lat_nr*columns + long_nr
    if not found:
        data = data.drop(data.index[i])
pixel_loc_list = []
for i in range(df_entries):
    if pixel_loc[i] > -1:
        pixel_loc_list.append(pixel_loc[i])

# pixel_loc_list links the order of appearance of observations in the raster to that of the dataframe. This allows us
# to copy information from the old locations to the new locations
data['index'] = pixel_loc_list

temp_discharge_raster = gdal.GetDriverByName('GTiff').Create(temp_discharge_location, reference_raster.RasterXSize,
                                                             reference_raster.RasterYSize, 1, gdal.GDT_Byte)
temp_discharge_raster.GetRasterBand(1).WriteArray(discharge_array)
temp_discharge_raster.SetProjection(reference_raster.GetProjectionRef())
temp_discharge_raster.SetGeoTransform(reference_raster.GetGeoTransform())
temp_discharge_raster = None

# find closest river points
graph_functions.print_graph(graph_location, [], reference_raster_location, river_raster_location)
river_discharge = shapefile_raster_functions.find_discharge(temp_discharge_location, river_raster_location, maximum_radius=3, precision=0.5)
river_discharge = river_discharge.transpose()

# move observation points to closest river
rows, columns = numpy.shape(river_discharge)
river_discharge = river_discharge.flatten()
discharge_array = discharge_array.flatten()
pixel_locations = []
pollution_data = []
latitudes = []
longitudes = []
index = -1

open_graph = open(graph_location, "rb")
river_graph = pickle.load(open_graph)
open_graph.close()

for i in range(len(discharge_array)):
    is_lake = True
    if discharge_array[i] == 1:  # if there is an observation point
        index += 1
        coordinates = river_discharge[i]
        if coordinates != ['none found']:  # if there is a river close to the observation point
            pixel_number = coordinates[1] * columns + coordinates[0]  # write the observation point
            if pixel_number in river_graph:
                while is_lake:  # move the observation point to the end of a lake, if it is in one.
                    is_lake = False
                    if river_graph.nodes[pixel_number]["lakes"] > 0:
                        is_lake = True
                        child = list(river_graph.successors(pixel_number))
                        try:
                            pixel_number = child[0]
                        except IndexError:
                            break

                pixel_locations.append(pixel_number)
                location = data['index'] == i
                pollution = data['pollutant'][location]
                pollution_data.append(pollution.iat[0])
                longitude = data['longitude'][location]
                longitudes.append(longitude.iat[0])
                latitude = data['latitude'][location]
                latitudes.append(latitude.iat[0])

# write output
df = pandas.DataFrame()
df['locations'] = pixel_locations
df['contaminant'] = pollution_data
df['longitude'] = longitudes
df['latitude'] = latitudes

df.to_csv(output_name)
