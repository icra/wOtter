import math
import os
import string
import networkx
import pandas
from osgeo import gdal
from osgeo import ogr
from typing import Union
import numpy
import random
import networkx as nx
import geopandas
import pickle
import time

def point_shapefile_sum_to_raster(shapefile_location: str, reference_raster_location: str,
                                  attributes: list = [], option: int = 1, output_name: str = '',
                                  print_info: bool = 0) -> Union[gdal.Band, numpy.ndarray, None]:
    """
    point_shapefile_sum_to_raster converts a shapefile with points to a raster of the size,
    resolution and projection of reference_raster. It sums the content of the shapefile.
     The raster can be returned as a matrix or as a .tif file
    :rtype: either a gdal.band object (a raster) or a numpy array
    :shapefile_location: string: the location of the shapefile on the computer
    :reference_raster_location: string: the location of a raster with the correct dimension on the computer
    :output_name: string: determines the name of the .tif file if a raster is returned
    :layer_number: int: the layer of the shapefile that contains the required information. Default is 0
    :option: bool: change to 0 for a matrix raster. Default (=1) gives a .tif file. Option 2 only saves
    :print_info: bool: change to 1 to print information
    :return: this function returns an indicator raster of locations of the shapefile. The return format
    can either be a matrix or a raster
    """

    # Loading the shapefile
    shapefile = ogr.Open(shapefile_location)
    shapefile_layer = shapefile.GetLayer(0)
    # Loading the directions raster for correct projection and resolution/mapping
    direction = gdal.Open(reference_raster_location)
    transformation = direction.GetGeoTransform()  # contains parameters for pixel to coordinate transformations

    # Collecting pixel locations of treatment centers
    pixel_location = []
    attributes_by_pixel = []
    empty = 0  # count empty geometries
    for feature in shapefile_layer:
        pt = feature.geometry()  # returns a geometry object
        if pt is not None:  # do not consider empty geometry objects
            pt_x = pt.GetX()
            pt_y = pt.GetY()
            inverse_transform = gdal.InvGeoTransform(transformation)  # coordinate to pixel instructions
            x_location, y_location = map(int, gdal.ApplyGeoTransform(inverse_transform, pt_x, pt_y))  # may contain
            # rounding error
            current_location = numpy.reshape([x_location, y_location], [1, 2])
            pixel_location.append([x_location, y_location])

            attributes_by_pixel.append([feature.GetField(attribute) for attribute in attributes])
        else:
            empty += 1  # counts the empty geometries

    if print_info == 1:  # if set to print
        print("There were " + str(len(pixel_location)) + " valid locations. " + str(empty) +
              " locations were dropped as they lacked coordinates.")

    # writing the pixels
    raster_data = numpy.zeros([direction.RasterXSize, direction.RasterYSize, 1 + len(attributes)], dtype=object)
    for i in range(len(pixel_location)):
        x_location = pixel_location[i][0]
        y_location = pixel_location[i][1]
        if 0 <= x_location <= direction.RasterXSize:  # the if statements ensure the pixel_location fits on the map
            if 0 <= y_location <= direction.RasterYSize:
                raster_data[x_location, y_location, 0] = 1  # write location to 1 if it contains a plant
                for j in range(len(attributes)):
                    raster_data[x_location, y_location, j+1] += attributes_by_pixel[i][j]

    if option == 0:
        return raster_data

    if output_name == '':  # returns a random string if no name is specified
        letters = string.ascii_lowercase
        output_name = ''.join(random.choice(letters) for i in range(8))  # filename has 8 letters
        output_name = output_name + ".tif"
        print("the output name was not specified, the randomly generated name is " + output_name)

    if not output_name.endswith('.tif'):  # ensures correct output format
        output_name += ".tif"
        print(".tif was added to the output name")

    # preparing output raster
    gtiff_driver = gdal.GetDriverByName('GTiff')
    out_ds = gtiff_driver.Create(output_name, direction.RasterXSize, direction.RasterYSize, len(attributes) + 1,
                                 gdal.GDT_Int32)
    out_ds.SetProjection(direction.GetProjection())  # copy direction projection to output raster
    out_ds.SetGeoTransform(direction.GetGeoTransform())  # copy direction resolution/location to output raster
    write_array = raster_data[:, :, 0]
    out_ds.GetRasterBand(1).WriteArray(write_array.transpose())
    for i in range(1, len(attributes)+1):
        write_array = raster_data[:, :, i]
        out_ds.GetRasterBand(i+1).WriteArray(write_array.transpose())
    if option == 1:
        return out_ds
    else:
        out_ds = None

    pass


def paint_rivers(raster_band: gdal.Band, direction_raster_location: str, layer_number: int = 0,
                 option: int = 1, output_name: str = '', print_info: bool = 0) -> Union[gdal.Band, numpy.ndarray, None]:
    """
    paint_rivers draws rivers using a shapefile with plant locations and a direction matrix. The size, resolution
    and projection of the river raster is copied from the direction raster. The raster can be returned as a matrix
    or a .tif file.
    :shapefile_location: string: the location of the shapefile on the computer
    :direction_raster_location: string: the location of the direction raster on the computer
    :layer_number: int: the layer of the shapefile that contains the required information. Default is 0
    :option: int: change to 0 for a matrix raster. Default (=1) gives a .tif file. 2 gives only save mode.
    :print_info: bool: change to 1 to print information
    :return: this function returns a river matrix based on initial locations and directions. The return format
    can either be a matrix or a raster. It can also
    """

    # obtain matrix raster of shapefile
    indicator_matrix = raster_band.ReadAsArray()
    indicator_matrix = indicator_matrix.transpose()
    # open direction matrix
    direction = gdal.Open(direction_raster_location)
    direction_matrix = direction.ReadAsArray()  # numpy array
    direction_matrix = direction_matrix.transpose()

    # create appropriate river matrix
    rows, columns = indicator_matrix.shape
    river_matrix = numpy.zeros([rows, columns])

    # for all entries in the river matrix
    for i in range(rows):
        for j in range(columns):
            # if the pixel contains a plant, but not a river
            if indicator_matrix[i, j] == 1 and not river_matrix[i, j] == 1:
                river_matrix[i, j] = 1  # river starting point
                current_x_cell = i
                current_y_cell = j

                river_continues = True
                while river_continues:  # keep painting the river if it continues
                    direction_value = direction_matrix[current_x_cell, current_y_cell]  # value determines direction
                    break_out = True  # breaks out of loop if river ends up at the edge of the raster

                    """" Now follow 9 distinct cases, pertaining to 8 direction and 1 invalid direction. If the
                     direction is invalid, stop (river_continues = false). If the direction is valid, check that
                      the next pixel in that direction is in the raster. If the next pixel is in the raster,
                       continue, otherwise stop (break_out is true -> river_continues is false).
                    """

                    if direction_value == 1:  # east
                        if current_x_cell < rows - 1:
                            current_x_cell += 1
                            break_out = False

                    elif direction_value == 2:  # SE
                        if current_x_cell < rows - 1 and current_y_cell < columns - 1:
                            current_x_cell += 1
                            current_y_cell += 1
                            break_out = False

                    elif direction_value == 4:  # S
                        if current_y_cell < columns - 1:
                            current_y_cell += 1
                            break_out = False

                    elif direction_value == 8:  # SW
                        if current_x_cell != 0 and current_y_cell < columns - 1:
                            current_x_cell += -1
                            current_y_cell += 1
                            break_out = False

                    elif direction_value == 16:  # W
                        if current_x_cell != 0:
                            current_x_cell += -1
                            break_out = False

                    elif direction_value == 32:  # NW
                        if current_x_cell != 0 and current_y_cell != 0:
                            current_x_cell += -1
                            current_y_cell += -1
                            break_out = False

                    elif direction_value == 64:  # N
                        if current_y_cell != 0:
                            current_y_cell += -1
                            break_out = False

                    elif direction_value == 128:  # NE
                        if current_x_cell < rows - 1 and current_y_cell != 0:
                            current_x_cell += 1
                            current_y_cell += -1
                            break_out = False

                    else:  # if there is no direction
                        river_continues = False

                    # if at the end of the raster, or if river loads into another
                    if break_out or river_matrix[current_x_cell, current_y_cell] == 1:
                        river_continues = False

                    river_matrix[current_x_cell, current_y_cell] = 1  # update the new (or unchanged) cell to 1

    if option == 0:  # return river as a matrix
        return river_matrix
    else:  # return river as a raster
        if output_name == '':  # returns a random string if no name is specified
            letters = string.ascii_lowercase
            output_name = ''.join(random.choice(letters) for i in range(8))  # filename has 8 letters
            output_name = output_name + ".tif"
            print("the output name was not specified, the randomly generated name is " + output_name)

        if not output_name.endswith('.tif'):  # ensures correct output format
            output_name += ".tif"
            print(".tif was added to the output name")

        gtiff_driver = gdal.GetDriverByName('GTiff')
        out_ds = gtiff_driver.Create(output_name, direction.RasterXSize, direction.RasterYSize, 1,
                                     gdal.GDT_Byte)
        out_ds.SetProjection(direction.GetProjection())  # copy direction projection to output raster
        out_ds.SetGeoTransform(direction.GetGeoTransform())  # copy direction resolution/location to output raster
        out_ds.GetRasterBand(1).WriteArray(river_matrix.transpose())
        if option == 1:
            return out_ds
        else:
            out_ds = None


def transpose_raster(input_raster_location: str, option: bool = 1, output_name: str = '') -> Union[None, gdal.Band]:
    """
    transpose_raster takes an input raster and transposes it. The raster can be returned or saved
    :input_raster_location: string: the location of the input raster on the computer
    :option: bool: change to 0 for saving but void. Default (=1) gives a gdal.band object
    :return: this function returns a transposed raster. The return format
    can either be a None (but saved) or a raster
    """
    input_raster = gdal.Open(input_raster_location)
    input_raster_array = input_raster.ReadAsArray()
    rows, columns = numpy.shape(input_raster_array)

    if output_name == '':  # returns a random string if no name is specified
        letters = string.ascii_lowercase
        output_name = ''.join(random.choice(letters) for i in range(8))  # filename has 8 letters
        output_name = output_name + ".tif"
        print("the output name was not specified, the randomly generated name is " + output_name)

    if not output_name.endswith('.tif'):  # ensures correct output format
        output_name += ".tif"
        print(".tif was added to the output name")

    gtiff_driver = gdal.GetDriverByName('GTiff')
    out_ds = gtiff_driver.Create(output_name, rows, columns, 1, gdal.GDT_Float32)
    out_ds.SetProjection(input_raster.GetProjection())  # copy direction projection to output raster
    out_ds.SetGeoTransform(input_raster.GetGeoTransform())  # copy direction resolution/location to output raster
    out_ds.GetRasterBand(1).WriteArray(input_raster_array.transpose())

    if option == 0:
        return out_ds
    else:
        out_ds = None


def multiply_rasters(indicator_raster_location: str, value_raster_location: str, option: int = 1,
                     output_name: str = '') -> Union[gdal.Band, numpy.ndarray, None]:
    """
    multiply_rasters multiplies two rasters elementwise. The raster can be returned as a matrix
    or a .tif file or can be saved.
    :indicator_raster_location: string: the location of the first raster on the computer
    :value_raster_location: string: the location of the second raster on the computer
    :option: int: change to 0 for a matrix raster. Default (=1) gives a .tif file. 2 gives only save mode.
    :return: this function returns a river matrix based on initial locations and directions. The return format
    can either be a matrix or a raster. It can also return none.
    """
    # opening rasters
    indicator_raster = gdal.Open(indicator_raster_location)
    indicator_raster_matrix = indicator_raster.ReadAsArray()
    indicator_raster_matrix = indicator_raster_matrix.transpose()
    value_raster = gdal.Open(value_raster_location)
    value_raster_matrix = value_raster.ReadAsArray()
    value_raster_matrix = value_raster_matrix.transpose()

    # multiplying the rasters
    output_matrix = numpy.multiply(indicator_raster_matrix, value_raster_matrix)

    # output options
    if option == 0:  # return output as a matrix
        return output_matrix
    else:  # as a raster
        if output_name == '':  # returns a random string if no name is specified
            letters = string.ascii_lowercase
            output_name = ''.join(random.choice(letters) for i in range(8))  # filename has 8 letters
            output_name = output_name + ".tif"
            print("the output name was not specified, the randomly generated name is " + output_name)

        if not output_name.endswith('.tif'):  # ensures correct output format
            output_name += ".tif"
            print(".tif was added to the output name")

        gtiff_driver = gdal.GetDriverByName('GTiff')
        out_ds = gtiff_driver.Create(output_name, indicator_raster.RasterXSize, indicator_raster.RasterYSize, 1,
                                     gdal.GDT_Byte)
        out_ds.SetProjection(indicator_raster.GetProjection())  # copy direction projection to output raster
        out_ds.SetGeoTransform(
            indicator_raster.GetGeoTransform())  # copy direction resolution/location to output raster
        out_ds.GetRasterBand(1).WriteArray(output_matrix.transpose())
        if option == 1:
            return out_ds
        else:
            out_ds = None


def shapefile_to_raster(input_shapefile: str, input_reference_raster: str, output_name: str = '',
                        attribute_name_list: list = [], option: int = 1, raster_number:int = 1,
                        all_attributes: bool = 0, scale: int = 1,
                        data_type: gdal.gdalconst = gdal.GDT_Float32, include_ind: bool = 1) -> \
                        Union[gdal.Band, numpy.ndarray, None]:
    """
    shapefile_to_raster converts an input_shapefile to a raster of the size,
    resolution and projection of input_reference_raster. The raster can be returned as a matrix
    or as a .tif file
    :rtype: either a gdal.band object (a raster) or a numpy array
    :input_shapefile: string: the location of the shapefile on the computer
    :input_reference_raster: string: the location of a raster with the correct dimension on the computer
    :output_name: string: determines the name of the .tif file if a raster is returned
    :attribute_name_list: string list: names of the attributes that the raster includes
    :option: bool: change to 0 for a matrix raster. Default (=1) gives a .tif file. Option 2 only saves
    :raster_number: int: if option is 0, this is the raster band you wish to have a matrix from.
    :all_attributes: bool: change to 1 to save all shapefile attributes in the raster
    :scale: int: determines the scale of the output raster. The output dimensions are (scale*rows) x (scale*columns)
    :return: this function returns an indicator raster of locations of the shapefile. The return format
    can either be a matrix, a raster or null.
    """

    burn_value = 1  # indicator for position

    ##########################################################
    # Get projection info from reference image
    reference_raster = gdal.Open(input_reference_raster, gdal.GA_ReadOnly)

    # Open Shapefile
    shapefile = ogr.Open(input_shapefile)
    print(input_shapefile)
    shapefile_layer = shapefile.GetLayer()

    # determine name of the output
    if output_name == '':
        letters = string.ascii_lowercase
        output_name = ''.join(random.choice(letters) for i in range(8))  # filename has 8 letters
        output_name = output_name + ".tif"
        print("the output name was not specified, the randomly generated name is " + output_name)

    if not output_name.endswith('.tif'):  # ensures correct output format
        output_name += ".tif"
        print(".tif was added to the output name")

    # all attributes option
    if all_attributes:
        attribute_name_list = []
        shapefile_layer_definition = shapefile_layer.GetLayerDefn()

        for i in range(shapefile_layer_definition.GetFieldCount()):
            attribute_name_list.append(shapefile_layer_definition.GetFieldDefn(i).name)

    if not attribute_name_list: # if only indicator, bytes
        data_type = gdal.GDT_Byte

    # Here the raster output is prepared with the dimensions, projection, name and datatype
    output = gdal.GetDriverByName('GTiff').Create(output_name, reference_raster.RasterXSize*scale,
                                                  reference_raster.RasterYSize*scale, len(attribute_name_list) +
                                                  include_ind, data_type)
    output.SetProjection(reference_raster.GetProjectionRef())
    transform = reference_raster.GetGeoTransform()

    if scale != 1:
        transform = list(transform)
        transform[1] = transform[1] / scale
        transform[5] = transform[5] / scale
        transform = tuple(transform)
        output.SetGeoTransform(transform)
    else:
        output.SetGeoTransform(reference_raster.GetGeoTransform())

    current_layer = 1
    # Next, the shapefile is converted into the newly created raster
    if include_ind:
        gdal.RasterizeLayer(output, [1], shapefile_layer, burn_values=[burn_value]) # position indicator
        current_layer += 1


    # Now fill the remaining layers with the attribute names provided
    for attribute in attribute_name_list:
        attribute_string = "ATTRIBUTE=" + attribute
        gdal.RasterizeLayer(output, [current_layer], shapefile_layer, options=[attribute_string])
        output.GetRasterBand(current_layer).SetDescription(attribute)
        current_layer += 1

    # return according to the option specified
    if option == 0:
        return output.GetRasterBand(raster_number).ReadAsArray()

    if option == 1:
        return output
    else:
        output = None
        pass


def find_discharge(point_raster_location: str, flow_raster_location: str, maximum_radius: float = 20,
                   precision: float = 1, print_info: bool = 0) -> numpy.ndarray:
    """
    find_discharge takes an indicator matrix of points and finds the closest flow in the flow_raster. It then stores
    the location of that closest flow in the indicator matrix pixel.
    :rtype: a numpy array that contains in each cell a 2-dim list of the coordinates
    :point_raster_location: string: the location of the point_raster on the computer
    :flow_raster_location: string: the location of a flow_raster on the computer
    :maximum_radius: float: gives the maximum distance between a point and the closest flow
    :precision: float: determines the maximum error in finding the closest river. E.g. with standard value 1, the
    distance between the actual closest river and the reported closest river can be at most 1 pixel.
    :return: this function returns the closest flow in the cell of the point raster.
    """

    # opening rasters
    point_raster = gdal.Open(point_raster_location)
    point_raster_layer = point_raster.GetRasterBand(1)
    point_raster_matrix = point_raster_layer.ReadAsArray()
    point_raster_matrix = point_raster_matrix.transpose()

    flow_raster = gdal.Open(flow_raster_location)
    flow_rasterband = flow_raster.GetRasterBand(1)
    flow_raster_matrix = flow_rasterband.ReadAsArray()
    flow_raster_matrix = flow_raster_matrix.transpose()

    # eliciting dimensions and creating output
    rows, columns = numpy.shape(point_raster_matrix)
    location_matrix = numpy.empty([rows, columns], object)

    point_count = 0  # counts amount of points in the point_raster
    no_river_found_count = 0  # counts points in the point_raster that have close river
    total_distance = 0 # the total distance from WWTP discharge points to the closest river

    # for all pixels in the point_raster
    for i in range(rows):
        for j in range(columns):
            # if the pixel is 1 initialize search
            if point_raster_matrix[i, j]:
                point_count += 1
                point_found = False
                r = 0
                while not point_found:
                    # take out a square around the point
                    rangeLow = int(r)
                    for x in range(-rangeLow, rangeLow + 1):
                        for y in range(-rangeLow, rangeLow + 1):
                            # if the point falls within raster dimensions and the radius
                            if 0 < i + x < rows and 0 < j + y < columns:
                                if x ** 2 + y ** 2 < r ** 2:

                                    # and the point has a river
                                    if flow_raster_matrix[i + x, j + y] == 1:
                                        point_found = True  # we break out of the loop
                                        location_matrix[i, j] = [i + x, j + y]  # we update the location
                                        distance = x ** 2 + y ** 2
                                        total_distance += math.sqrt(distance)

                    r += precision  # increment radius, relevant if no point is found.

                    if r > maximum_radius:
                        location_matrix[i, j] = ["none found"]
                        no_river_found_count += 1
                        break

    if print_info:
        print("There were " + str(no_river_found_count) + " points for which no flow was found within the maximum "
                                                          "radius out of " + str(point_count) + " total points")
        print("the average distance for connected discharge points to an actual river is " + str(total_distance/point_count))

    return location_matrix


def join_indicator(indicator_raster1: str, indicator_raster2: str, output_name: str = '',
                   option: int = 1) -> Union[numpy.ndarray, gdal.Band, None]:

    """
    join_indicator takes two indicator matrices and returns an indicator matrix that is 1 if at least for one
    input indicator matrix the cell is one.
    :rtype: Either a band (raster), a numpy array or none
    :indicator_raster1: string: the location of the first indicator matrix on the computer
    :indicator_raster2: string: the location of the second indicator matrix on the computer
    :output_name: string: determines the file name of the output
    :option: int: if 0, the function returns a numpy array. If 1, a raster and if 2 it only saves.
    :return: This function returns an indicator matrix/raster/file.
    """

    # opening rasters
    indicator_raster1 = gdal.Open(indicator_raster1)
    indicator_raster1_matrix = indicator_raster1.ReadAsArray()
    indicator_raster1_matrix = indicator_raster1_matrix.transpose()

    indicator_raster2 = gdal.Open(indicator_raster2)
    indicator_raster2_matrix = indicator_raster2.ReadAsArray()
    indicator_raster2_matrix = indicator_raster2_matrix.transpose()

    output_matrix = numpy.add(indicator_raster1_matrix, indicator_raster2_matrix) # adds matrices
    output_matrix = numpy.where(output_matrix <2, 0, 1) # equals 0 if entry is smaller than 0.5

    # output options
    if option == 0:  # return output as a matrix
        return output_matrix
    else:  # as a raster
        if output_name == '':  # returns a random string if no name is specified
            letters = string.ascii_lowercase
            output_name = ''.join(random.choice(letters) for i in range(8))  # filename has 8 letters
            output_name = output_name + ".tif"
            print("the output name was not specified, the randomly generated name is " + output_name)

        if not output_name.endswith('.tif'):  # ensures correct output format
            output_name += ".tif"
            print(".tif was added to the output name")

        gtiff_driver = gdal.GetDriverByName('GTiff')
        out_ds = gtiff_driver.Create(output_name, indicator_raster1.RasterXSize,
                                     indicator_raster1.RasterYSize, 1, gdal.GDT_Byte)
        out_ds.SetProjection(indicator_raster1.GetProjection())  # copy direction projection to output raster
        out_ds.SetGeoTransform(
            indicator_raster1.GetGeoTransform())  # copy direction resolution/location to output raster
        out_ds.GetRasterBand(1).WriteArray(output_matrix.transpose())
        if option == 1:
            return out_ds
        else:
            out_ds = None


def reposition_raster(initial_raster: str, location_matrix: numpy.ndarray, output_name: str = '') -> None:
    """
    :param initial_raster: string: location of the raster-file that needs to be repositioned on the computer
    :param location_matrix: ndarray: array of the correct size that contains redirection coordinates in its matrix cells
    :param output_name: string: the resulting name and location of the raster file
    :return: this function creates a raster file, but returns nothing.
    """
    initial_raster = gdal.Open(initial_raster)

    # creating an output name if forgotten
    if output_name == '':  # returns a random string if no name is specified
        letters = string.ascii_lowercase
        output_name = ''.join(random.choice(letters) for i in range(8))  # filename has 8 letters
        output_name = output_name + ".tif"
        print("the output name was not specified, the randomly generated name is " + output_name)

    if not output_name.endswith('.tif'):  # ensures correct output format
        output_name += ".tif"
        print(".tif was added to the output name")

    # Rasterise
    output = gdal.GetDriverByName('GTiff').Create(output_name, initial_raster.RasterXSize,
                                                  initial_raster.RasterYSize, initial_raster.RasterCount,
                                                  gdal.GDT_Int16, options=['COMPRESS=DEFLATE'])
    output.SetProjection(initial_raster.GetProjectionRef())
    output.SetGeoTransform(initial_raster.GetGeoTransform())
    dropped_WWTP_count = 0

    for i in range(1, initial_raster.RasterCount + 1):
        current_band = initial_raster.GetRasterBand(i)
        input_array = current_band.ReadAsArray()
        input_array = input_array.transpose()
        rows, columns = numpy.shape(input_array)
        output_array = numpy.zeros([rows, columns])
        for j in range(rows):
            for k in range(columns):
                if location_matrix[j, k] is not None:
                    if len(location_matrix[j, k]) > 1:
                        location_list = location_matrix[j, k]
                        x_location = location_list[0]
                        y_location = location_list[1]
                        if output_array[j, k] == 0:
                            output_array[x_location, y_location] = input_array[j, k]
                        else:
                            dropped_WWTP_count += 1 # this may increment if a WTTP would overwrite another one in the
                            # location

        output.GetRasterBand(i).WriteArray(output_array.transpose())
    if dropped_WWTP_count > 0:
        print(dropped_WWTP_count + " were dropped as their location coincided with another WTTP")
    output = None
    pass


def next_cell(starting_cell, direction, nrows, ncolumns):

    y = starting_cell[0]
    x = starting_cell[1]
    if direction == 1:  # east
        if x < ncolumns - 1:
            x += 1

    if direction == 2:  # SE
        if x < ncolumns - 1 and y < nrows - 1:
            x += 1
            y += 1

    if direction == 4:  # S
        if y < nrows - 1:
            y += 1

    if direction == 8:  # SW
        if x != 0 and y < nrows - 1:
            x += -1
            y += 1

    if direction == 16:  # W
        if x != 0:
            x += -1

    if direction == 32:  # NW
        if x != 0 and y != 0:
            x += -1
            y += -1

    if direction == 64:  # N
        if y != 0:
            y += -1

    if direction == 128:  # NE
        if x < ncolumns - 1 and y != 0:
            x += 1
            y += -1
    return [y,x]


def encode_direction(starting_cell, next_cell):
    y = starting_cell[0]
    x = starting_cell[1]
    direction = -10
    if isinstance(next_cell, tuple):
        next_cell = [i for i in next_cell]

    if next_cell == [y, x+1]:  # east
        direction = 1

    if next_cell == [y+1, x+1]:  # SE
        direction = 2

    if next_cell == [y+1, x]:  # S
        direction = 4

    if next_cell == [y+1, x-1]:  # SW
        direction = 8

    if next_cell == [y, x-1]:  # W
        direction = 16

    if next_cell == [y-1, x-1]:  # NW
        direction = 32

    if next_cell == [y-1, x]:  # N
        direction = 64

    if next_cell == [y-1, x+1]:  # NE
        direction = 128
    return direction


def shorter_path_2cells(cell_0, cell_f, matrix_river01):

    G = nx.Graph()
    i_dir = (16,32,64,128)
    rows, columns = numpy.shape(matrix_river01)

    for i in range(rows):
        for j in range(columns):
            if matrix_river01[i, j] == 1:
                G.add_node((i, j))
                for direction in i_dir:
                    adjacent_cell = next_cell([i,j], direction, rows, columns)
                    if matrix_river01[adjacent_cell[0], adjacent_cell[1]] == 1:
                        G.add_edge((i,j), tuple(adjacent_cell))

    path = nx.shortest_path(G, cell_0, cell_f)

    return path


def flow_through(matrix, direction):
    rows, columns = numpy.shape(matrix)
    found_incoming = False
    for i in range(rows):
        for j in range(columns):
            succeeding_cell = next_cell([i,j], direction[i,j], rows, columns)
            if matrix[i,j] and succeeding_cell == [1,1]:
                found_incoming = True
            elif matrix[i,j]:
                outgoing_cell = [i,j]
    if found_incoming:
        return encode_direction([1,1], outgoing_cell)
    else:
        pass


def cell_dimensions(reference_raster_location:str) -> list:
    """
    Length_cell creates 3 numpy vectors that contain information on the horizontal, vertical and diagonal dimensions of
    a cell, depending on the row of a raster with appropriate dimensions
    :return: returns 3 vectors that contain the horizontal, vertical and diagonal distances of a cell.
    """
    reference_raster = gdal.Open(reference_raster_location)
    rows = reference_raster.RasterYSize
    transform = reference_raster.GetGeoTransform()
    latitude_step = abs(transform[5])
    lowest_latitude = transform[3] - latitude_step * rows
    radius = 6371007.2  # terrestrial radius

    # Calculation of the horizontal distance
    horizontal = numpy.zeros([rows, 1])
    for i in range(rows):
        latitude = lowest_latitude + latitude_step * (rows - i - 1)

        horizontal[i] = (numpy.sin((latitude + latitude_step / 2) * numpy.pi / 180) - numpy.sin(
            (latitude - latitude_step / 2) * numpy.pi / 180)) * radius  # cells width

    # Calculation of vertical and diagonal distances
    vertical = numpy.zeros([rows, 1]) + latitude_step * numpy.pi / 180 * radius # this value is constant
    diagonal = numpy.sqrt(numpy.square(vertical) + numpy.square(horizontal)) # Pythagoras' theorem

    distances = [horizontal, vertical, diagonal]
    return distances

def crop_rasters(raster_locations: list, upper_left_x: float = -10**9, upper_left_y: float = 10**9,
                 lower_right_x: float = 10**9, lower_right_y: float = -10**9) -> list:
    """
    Crop rasters takes in raster directory locations, and converts those rasters into new rasters that have the
    smallest common dimension between them. Dimensions may be capped by specifying upper_left_x, y and lower_right_x
    and y.
    :rtype: None
    :raster_locations: list: List of directory locations of the rasters
    :upper_left_x: float: gives the x-coordinate of the upper left part of the raster
    :upper_left_y: float: gives the y-coordinate of the upper left part of the raster
    :lower_right_x: float: gives the x-coordinate of the lower right part of the raster
    :lower_right_y: float: gives the y-coordinate of the lower right part of the raster
    :return: returns names of output. Saves rasters under the name crop_ + name
    """


    rasters = [gdal.Open(raster_locations[i]) for i in range(len(raster_locations))]
    for raster in rasters:  # this loop crops the bounds additionally if one raster is too small for the custom bounds

        transform = raster.GetGeoTransform()


        # get the most stringent bounds
        upper_left_x = max(transform[0], upper_left_x)
        upper_left_y = min(transform[3], upper_left_y)
        lower_right_x = min(transform[0] + transform[1] * raster.RasterXSize, lower_right_x)
        lower_right_y = max(transform[3] + transform[5] * raster.RasterYSize, lower_right_y)

    # get the largest bounds
    window = (upper_left_x, upper_left_y, lower_right_x, lower_right_y)
    output_names = []
    for name in raster_locations:
        stem_name = name.split('.')
        stem_name = stem_name[0]
        output_name = stem_name + "_temp.tif"
        output_names.append(output_name)
        gdal.Translate(output_name, name, projWin=window, options=['COMPRESS=DEFLATE'])

    return output_names

def reset_datatype_raster(input_raster_location: str, datatype=gdal.GDT_Byte) -> None:
    """
    This function takes in a raster and a datatype. The old raster is replaced with a new raster based on the new
    datatype.
    :rtype: None
    :input_raster_location: str: location of the raster
    :datatype: gdal datatype: gdal datatype
    """
    input_raster = gdal.Open(input_raster_location)

    output = gdal.GetDriverByName('GTiff').Create("temp_loc.tif", input_raster.RasterXSize,
                                                  input_raster.RasterYSize, input_raster.RasterCount,
                                                  datatype, options=['COMPRESS=DEFLATE'])
    output.SetProjection(input_raster.GetProjectionRef())
    output.SetGeoTransform(input_raster.GetGeoTransform())

    for band in range(1, input_raster.RasterCount + 1):
        input_band = input_raster.GetRasterBand(band).ReadAsArray()
        output.GetRasterBand(band).WriteArray(input_band)
    output = None
    input_raster = None
    os.remove(input_raster_location)
    os.rename("temp_loc.tif", input_raster_location)
    pass


def give_pixel(coord: list, reference_raster: object, return_scalar: bool = False, reverse: bool = False) -> Union[int,
               list]:
    """
    This function gives the pixel of a coordinate (row/latitude, column/longitude). If reverse is specified, the function returns a coordinate for a
    pixel number.
    :rtype: location of the pixel as a list or as the pixel number. If reverse is specified, gives coordenates as a list
    :coord: list: The coordinations of the point. If reverse is specified, this needs to be the pixel number
    :reference_raster: object: A raster with the desired dimensions
    :return_scalar: bool: if true, the output is returned as a scalar (pixel number)
    :reverse: bool: if true, the function takes in a pixel number and returns a coordinate
    """
    transformation = reference_raster.GetGeoTransform()
    if not reverse:
        inverse_transform = gdal.InvGeoTransform(transformation)  # coordinate to pixel instructions
        long_location, lat_location = map(int, gdal.ApplyGeoTransform(inverse_transform, coord[1], coord[0]))

        if return_scalar:
            pixel_number = lat_location * reference_raster.RasterXSize + long_location
            return pixel_number
        return [lat_location, long_location]
    i = int(coord/reference_raster.RasterXSize)
    j = coord - reference_raster.RasterXSize * i
    long_location, lat_location = gdal.ApplyGeoTransform(transformation, int(j), i)
    return [lat_location, long_location]

def csv_to_shapefile(contaminant_location: str, reference_raster_location: str, output_name: str = '', field:
                             str='locations', options: bool = 0) -> None:
    """
       This function makes a shapefile from an observation polution dataframe.
       :rtype: void, saves a shapefile in the same folder with the same name as the dataframe.
       :contaminant_location: str: name of the contaminant dataframe
       :reference_raster_location: str: location of a reference raster with dimensions that concord with 'pixel_number'
       """
    if not isinstance(contaminant_location, str) and output_name == '':
        contamination = contaminant_location
        output_location = "shapefile.shp"
    elif not isinstance(contaminant_location, str) and output_name != '':
        contamination = contaminant_location
        output_location = output_name
    else:
        contamination = pandas.read_csv(contaminant_location)
        contaminant_location_split = contaminant_location.split('.')
        output_location = contaminant_location_split[0] + '.shp'

    # create the point shapefile
    driver = ogr.GetDriverByName('Esri Shapefile')
    ds = driver.CreateDataSource(output_location)
    layer = ds.CreateLayer('', None, ogr.wkbPoint)

    # create fields
    for name in contamination.columns:
        layer.CreateField(ogr.FieldDefn(name, ogr.OFTReal))
    defn = layer.GetLayerDefn()  # stores the type of features this layer has
    reference_raster = gdal.Open(reference_raster_location)
    # Conversion from dataframe to shapefile
    for i in range(len(contamination)):
        # store location
        if options:
            latitude = contamination.iloc[i, 0]
            longitude = contamination.iloc[i, 1]
        else:
            latitude, longitude = give_pixel(contamination[field].iloc[i], reference_raster, reverse=True)

        feature = ogr.Feature(defn)  # create the feature if a location was found

        # next, link the point with the feature
        point = ogr.Geometry(ogr.wkbPoint)
        point.AddPoint(longitude, latitude)
        feature.SetGeometry(point)

        # Store the information
        for j in range(len(contamination.columns)):
            column_name = contamination.columns[j]
            if column_name != field:
                name = column_name[0:10]
                feature.SetField(name, float(contamination[column_name][i]))

        layer.CreateFeature(feature)  # add the feature to the layer
        feature = point = None  # empty the feature and point for next iteration

    ds = layer = None  # saves the layer and the shapefile
    pass


def split_raster(raster_location: str, split_horizontal: int = 1, split_vertical: int = 1, output_name: str ="temp") \
        -> None:
    """
       This function splits a raster into seperate parts specified by the horizontal and vertical splits.
       :rtype: void, rasters are saved
       :raster_location: str: location of the raster on the pc
       :split_horizontal: str: amount of horizontal cuts
       :split_vertical: str: amount of vertical cuts
       :output_name: str: name of the output (excluding the numbering of the output)
       """

    ds = gdal.Open(raster_location)
    band = ds.GetRasterBand(1)
    xsize = band.XSize
    ysize = band.YSize
    matrix = band.ReadAsArray()
    no_data_val = band.GetNoDataValue()
    split_x = int(xsize/split_vertical)
    split_y = int(ysize/split_horizontal)
    split_x_values =[i * split_x for i in range(split_vertical)]
    split_y_values = [i * split_y for i in range(split_horizontal)]
    split_x_values += [xsize]
    split_y_values += [ysize]

    raster_count = (len(split_y_values) - 1) * (len(split_x_values) - 1)
    for i in range(len(split_x_values) - 1):
        for j in range(len(split_y_values) - 1):
            output_matrix = matrix[split_y_values[j] : split_y_values[j+1] : 1,
                            split_x_values[i] : split_x_values[i+1] : 1]
            RasterXSize = split_x_values[i+1] - split_x_values[i]
            RasterYSize =split_y_values[j+1] - split_y_values[j]
            output = gdal.GetDriverByName('GTiff').Create(output_name + str(i*(len(split_y_values)-1) + j) +'.tif', RasterXSize, RasterYSize, 1,
                                                          gdal.GDT_Float64, options=['COMPRESS=DEFLATE'])
            output.SetProjection(ds.GetProjectionRef())
            lat, long = give_pixel(split_y_values[j] * xsize + split_x_values[i], ds, reverse=True)
            transform = list(ds.GetGeoTransform())
            transform[0] = long
            transform[3] = lat
            output.SetGeoTransform(transform)
            output.GetRasterBand(1).WriteArray(output_matrix)
            output.GetRasterBand(1).SetNoDataValue(no_data_val)

    pass




