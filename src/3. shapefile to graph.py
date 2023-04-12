from osgeo import gdal
import numpy
from src.library import shapefile_raster_functions
import networkx
import pickle
import os
import geopandas
from src.library import graph_functions as gf
from time import time

directory = os.path.join(os.getcwd(), 'data')

# input files
rivers_shapefile = os.path.join(directory, "hydro_rivers_adapted.shp")

# gives the dimensions/projections of the output
reference_raster_location = os.path.join(directory, "reference_raster.tif")

# directions corresponding to the hydrorivers shapefile
direction_raster_location = os.path.join(directory, "15s_directions.tif")
slopes_raster_location = os.path.join(directory, "15s_slopes_10km.tif")  # calculated in calculate slope.py

# shapefile of lakes
lakes_shapefile_location = os.path.join(directory, "Raw data/HydroSHEDS", "HydroLakes_polys_v10.shp")

# temporary output
id_location = "id.tif"  # will contain the river ids
discharge_location = "discharge.tif"  # will contain the discharges of the rivers
indicator_location = "indicator.tif"  # will contain an indicator if a pixel is a river.
shapefile_id_location = "shapefile_id.shp"  # will contain the rivers shapefile with an additional unique id.
lakes_location = "lakes.tif"  # will contain an indicator if a pixel is a lake.
indicator_location = os.path.join(directory, "3s_rivers.tif")

# output files
graph_location = os.path.join(directory, "river_graph.pkl")
topological_sort_location = os.path.join(directory, "sorted_river_list.pkl")
rivers_15s_location = os.path.join(directory, "15s_rivers.tif")
rivers_from_graph_location = os.path.join(directory, "rivers_from_graph.tif")

# parameters of input
slopes_factor = 1/10000  # needs to concord with 'calculate slope.py'
scale_factor = 5  # increases the precision of rasterizing the shapefile
minimum_discharge = 0.01  # ensures that residence time does not become infinite.

# Change the river shapefile and give it an additional id. The HYRIV_ID is too large to be rasterized without precision
# error.

frame_shapefile = geopandas.read_file(rivers_shapefile)
frame_shapefile["ID"] = [i for i in range(len(frame_shapefile))]  # add an additional id for each shapefile object
frame_shapefile.to_file(shapefile_id_location)  # Save the new shapefile

# rasterize shapefiles to a smaller resolution. The smaller resolution ensures that the shapefile and the raster
# correspond precisely. if a larger resolution is chosen, the rasterized pixel might be adjacent to where the shapefile
# is
reference_raster = gdal.Open(reference_raster_location)



# put the unique river ids into a raster
river_id = shapefile_raster_functions.shapefile_to_raster(shapefile_id_location, reference_raster_location, id_location,
                                  attribute_name_list=["ID"], option=0, scale=scale_factor,
                                  include_ind=False, data_type=gdal.GDT_Byte)
# convert the raster to 15 seconds
rows, columns = numpy.shape(river_id)
rows = int(rows/scale_factor)
columns = int(columns/scale_factor)
river_id_matrix = river_id.reshape([rows, scale_factor, columns, scale_factor])

# put the discharges into a raster
river_discharge = shapefile_raster_functions.shapefile_to_raster(rivers_shapefile, reference_raster_location, discharge_location,
                                         attribute_name_list=["DIS_AV_CMS"], option=0, scale=scale_factor,
                                         include_ind=False)
river_discharge = river_discharge.reshape([rows, scale_factor, columns, scale_factor])

# put the river indicators into a raster
river_indicator = shapefile_raster_functions.shapefile_to_raster(rivers_shapefile, reference_raster_location, indicator_location, option=0,
                                         scale=scale_factor, data_type=gdal.GDT_Byte)
reduced_river_matrix = river_indicator.reshape([rows, scale_factor, columns, scale_factor]).sum(3).sum(1)

# put the lakes into a raster with their id and their total volume.
shapefile_raster_functions.shapefile_to_raster(lakes_shapefile_location, reference_raster_location, lakes_location,
                       attribute_name_list=['Hylak_id', 'Vol_total'], data_type=gdal.GDT_Int16)


# If after downscaling, a pixel would contain more than 2 river pixels in the higher resolution raster, we consider it a
# river in the downscaled raster as well
reduced_river_matrix = reduced_river_matrix > 2

# prepare inputs for the graph
rows, columns = numpy.shape(reduced_river_matrix)
river_graph = networkx.DiGraph()  # the graph object that will be the output
direction_raster = gdal.Open(direction_raster_location)  # determines the directed edges of the graph
slopes_raster = gdal.Open(slopes_raster_location)
lakes_raster = gdal.Open(lakes_location)

direction_matrix = direction_raster.GetRasterBand(1).ReadAsArray()
slopes_matrix = slopes_factor * slopes_raster.GetRasterBand(1).ReadAsArray()
lakes_id_matrix = lakes_raster.GetRasterBand(2).ReadAsArray()
lakes_vol_matrix = lakes_raster.GetRasterBand(3).ReadAsArray()

horizontal_distance, vertical_distance, diagonal_distance = shapefile_raster_functions.cell_dimensions(reference_raster_location)
# create graph
for i in range(rows):
    for j in range(columns):  # for all element in the river matrix
        if reduced_river_matrix[i, j] == 1 and direction_matrix[i, j] > 0:  # if the pixel is a river

            # collect information on largest discharge and associated river id.
            current_cell_number = i * columns + j  # this becomes the main identifier of the node.

            # When downscaling, we take the maximum discharge
            discharge_candidates = river_discharge[i, :, j, :]
            if numpy.max(discharge_candidates) == 0:
                loc = numpy.argmax(river_id_matrix[i, :, j, :])
                row_loc = int(loc / scale_factor)
                column_loc = loc - scale_factor * row_loc
            else:
                discharge_loc = numpy.argmax(river_discharge[i, :, j, :])
                row_loc = int(discharge_loc / scale_factor)
                column_loc = discharge_loc - scale_factor * row_loc

            # assign the river id corresponding to that maximum discharge
            river_id_val = river_id_matrix[i, row_loc, j, column_loc]

            # swaps the id for the corresponding HYRIV_ID
            river_id_val = frame_shapefile["HYRIV_ID"].iloc[int(river_id_val)] - 1

            discharge_val = river_discharge[i, row_loc, j, column_loc]
            if discharge_val < minimum_discharge:
                discharge_val = minimum_discharge

            # find cell distance
            direction = direction_matrix[i, j]
            if direction in (1, 16):
                distance = horizontal_distance[i]
            elif direction in (4, 64):
                distance = vertical_distance[i]
            else:
                distance = diagonal_distance[i]

            # residence time for rivers
            slope_val = slopes_matrix[i, j]
            RT = gf.calculate_residence_time(discharge_val, slope_val, distance)
            discharge_hour = discharge_val * 3600
            long, lat = shapefile_raster_functions.give_pixel(current_cell_number, reference_raster, reverse=1)
            # add node with all the data
            river_graph.add_node(current_cell_number, x=i, y=j, longitude=long, latitude=lat, slope=slope_val,
                                 lakes=lakes_id_matrix[i, j], volume=lakes_vol_matrix[i, j], flow_HR=discharge_hour,
                                 HYRIV_ID=river_id_val, cell_distance=float(distance), RT_HR=float(RT))

            # create edge if a next cell exists
            next_cell = shapefile_raster_functions.next_cell([i, j], direction, rows, columns)  # gives next cell according to the direction
            next_cell_number = columns * next_cell[0] + next_cell[1]

            # add edge to next river cell, if that river cell exists
            if reduced_river_matrix[next_cell[0], next_cell[1]] == 1 and current_cell_number != next_cell_number:
                river_graph.add_edge(current_cell_number, next_cell_number)

# the next code removes river nodes that do not fit into the reference raster, but were included by the previous code as
# the result of an edge
remove_list = []
for node in river_graph:
    if not river_graph[node]:
        remove_list.append(node)
river_graph.remove_nodes_from(remove_list)

sorted_graph = list(networkx.topological_sort(river_graph))  # gives a list where nodes first in the list are preceding
gf.add_RT_lakes(river_graph, sorted_graph, "RT_HR")  # this adds the residence times of lakes

# saving the output
open_file = open(graph_location, "wb")
pickle.dump(river_graph, open_file)
open_file.close()

open_file = open(topological_sort_location, "wb")
pickle.dump(sorted_graph, open_file)
open_file.close()
# deleting temporary output
river_id = None
lakes_raster = None
os.remove(id_location)
os.remove(discharge_location)
shapefile_id_location = shapefile_id_location.split('.')
shapefile_id_location = shapefile_id_location[0]
for extension in ['.shp', '.shx', '.dbf', '.cpg', '.prj']:
    os.remove(shapefile_id_location + extension)
os.remove(lakes_location)

# This creates additional otuput

# rivers from the graph
gf.print_graph(graph_location, [], reference_raster_location, rivers_from_graph_location)

# create 15s river raster
gtiff_driver = gdal.GetDriverByName('GTiff')
out_ds = gtiff_driver.Create(rivers_15s_location, reference_raster.RasterXSize, reference_raster.RasterYSize,
                             1, gdal.GDT_Byte)  # this creates a raster document with dimensions, bands, datatype

out_ds.SetProjection(reference_raster.GetProjection())  # copy direction projection to output raster
out_ds.SetGeoTransform(reference_raster.GetGeoTransform())  # copy direction resolution/location to output raster
out_ds.GetRasterBand(1).WriteArray(reduced_river_matrix)
out_ds = None