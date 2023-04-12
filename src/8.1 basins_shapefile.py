from osgeo import gdal, ogr, osr
import numpy
import library.shapefile_raster_functions as sh
import networkx
import pickle
from time import time
import os
import geopandas

from src.library import graph_functions


directory = os.path.join(os.getcwd(), 'data')

# input files
direction_raster_location = os.path.join(directory, "15s_directions.tif")  # directions hydroRIVERS shapefile
graph_location = os.path.join(directory, 'river_graph.pkl')

# output files
basins_raster_location = os.path.join(directory, "connected_basins.tif")
river_basins_location = os.path.join(directory, 'river_basins.tif')
river_basins_shape_location = os.path.join(directory, 'river_basins.shp')

open_graph = open(graph_location, "rb")
river_graph = pickle.load(open_graph)
open_graph.close()

undirected_rivers = river_graph.to_undirected()
basins = list(networkx.connected_components(undirected_rivers))
networkx.set_node_attributes(river_graph, 0, name='basin')
for index in range(len(basins)):
    for node in basins[index]:
        river_graph.nodes[node]['basin'] = index

graph_functions.print_graph(river_graph, ['basin'], direction_raster_location, river_basins_location, gdal.GDT_Int16)
null_basin = len(basins) + 1
# Load data
direction_raster = gdal.Open(direction_raster_location)  # determines the directed edges of the graph
direction_matrix = direction_raster.GetRasterBand(1).ReadAsArray()
rows, columns = numpy.shape(direction_matrix)
basin_raster = gdal.Open(river_basins_location)
basin_matrix = basin_raster.GetRasterBand(1).ReadAsArray()
for i in range(rows):
    for j in range(columns):
        if basin_matrix[i, j] < 0 and direction_matrix[i, j] > 0:
            current_network = [[i, j]]
            basin_matrix[i, j] = null_basin
            search_for_next_cell = True
            current_i = i
            current_j = j
            while search_for_next_cell:
                direction = direction_matrix[current_i, current_j]
                next_cell = sh.next_cell([current_i, current_j], direction, rows, columns)
                if next_cell != [current_i, current_j]:
                    current_i, current_j = next_cell
                    if basin_matrix[current_i, current_j] < 0:
                        basin_matrix[current_i, current_j] = null_basin
                        current_network.append([current_i, current_j])
                    else:
                        search_for_next_cell = False
                        value = basin_matrix[current_i, current_j]
                        for cell in current_network:
                            basin_matrix[cell[0], cell[1]] = value
                else:
                    search_for_next_cell = False

basin_matrix[basin_matrix == null_basin] = -1
gtiff_driver = gdal.GetDriverByName('GTiff')
out_ds = gtiff_driver.Create(basins_raster_location, direction_raster.RasterXSize, direction_raster.RasterYSize,
                             1, gdal.GDT_Int32)  # this creates a raster document with dimensions, bands, datatype

out_ds.SetProjection(direction_raster.GetProjection())  # copy direction projection to output raster
out_ds.SetGeoTransform(direction_raster.GetGeoTransform())  # copy direction resolution/location to output raster
out_ds.GetRasterBand(1).WriteArray(basin_matrix)
out_ds.GetRasterBand(1).SetNoDataValue(-1)
out_ds = None

open_file = open(graph_location, "wb")
pickle.dump(river_graph, open_file)
open_file.close()

#  get raster datasource
src_ds = gdal.Open(basins_raster_location)
#
srcband = src_ds.GetRasterBand(1)
dst_layername = 'basin'
drv = ogr.GetDriverByName("ESRI Shapefile")
dst_ds = drv.CreateDataSource(river_basins_shape_location)

sp_ref = osr.SpatialReference()
sp_ref.SetFromUserInput('EPSG:4326')

dst_layer = dst_ds.CreateLayer(dst_layername, srs=sp_ref)

basin_number = ogr.FieldDefn("basin_numb", ogr.OFTInteger)
dst_layer.CreateField(basin_number)
dst_field = dst_layer.GetLayerDefn().GetFieldIndex("basin_numb")

gdal.Polygonize(srcband, srcband, dst_layer, dst_field, [], callback=None)

del src_ds
del dst_ds

