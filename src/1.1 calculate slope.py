from osgeo import gdal
import numpy
import os
from time import time



directory = os.path.join(os.getcwd(), 'data')
print(directory)
# input rasters

height_raster_location = os.path.join(directory, "Raw data/3s_height.tif")
direction_raster_upscale_location = os.path.join(directory, "Raw data/15s_directions.tif")

# output rasters
slopes_15_location = os.path.join(directory, "15s_slopes_10km.tif")
reference_raster_location = os.path.join(directory, "reference_raster.tif")

# cut direction_raster height_raster




output_names = shapefile_raster_functions.crop_rasters([height_raster_location, direction_raster_upscale_location],
                            lower_right_x=36, lower_right_y=33)


height_raster_location = os.path.join(directory, "3s_height.tif")
direction_raster_upscale_location = os.path.join(directory, "15s_directions.tif")
os.rename(output_names[0], height_raster_location)
os.rename(output_names[1], direction_raster_upscale_location)

# parameters
scale_factor = 10000  # multiplies the slopes such that they can be stored as integers with sufficient precision
downscale_factor = 5  # determines the ratio between input resolution and output resolution (e.g. from 15 seconds to 3
                      # seconds)
minimum_slope = 5  # this gives a lower bound for the slope that is allowed
                   # (minimum_slope/scale_factor * 100% gives the unit in percentages)

# load height and directions
print('1')
height_raster = gdal.Open(height_raster_location)
print('2')
height_matrix = height_raster.ReadAsArray().astype(numpy.int16)
print('3')
rows, columns = numpy.shape(height_matrix)

print('aacsdfsdfdsfv')

# Reserving memory for the matrices
height_difference = numpy.zeros([rows, columns], dtype=numpy.int16)  # int 16 reduces memory use
direction_matrix = numpy.zeros([rows, columns], dtype=numpy.byte)
minus_one_canvas = numpy.zeros([rows, columns], dtype=numpy.int16) - 1
minus_one_canvas = direction_matrix - 1
# calculating the height differences using the directions. Similar code is used for all 8 directions.
height_matrix[height_matrix< -10] = -10
# east
candidates = minus_one_canvas

candidates[:, 0:columns-1:1] = height_matrix[:, 0:columns-1:1] - height_matrix[:, 1:columns:1]
indicators = candidates >= 0
height_difference[indicators] = candidates[indicators]
direction_matrix[indicators] = 1


print('aaaaaaaaaaaaaaaaaaa')


# south-east
candidates = minus_one_canvas
candidates[0:rows-1:1, 0:columns-1:1] = height_matrix[0:rows-1:1, 0:columns-1:1] - height_matrix[1:rows:1, 1:columns:1]
indicators = candidates > height_difference
height_difference[indicators] = candidates[indicators]
direction_matrix[indicators] = 2

# south
candidates = minus_one_canvas
candidates[0:rows-1:1, :] = height_matrix[0:rows-1:1,:] - height_matrix[1:rows:1, :]
indicators = candidates > height_difference
height_difference[indicators] = candidates[indicators]
direction_matrix[indicators] = 4

# south-west
candidates = minus_one_canvas
candidates[0:rows-1:1, 1:columns:1] = height_matrix[0:rows-1:1, 1:columns:1] - height_matrix[1:rows:1, 0:columns-1:1]
indicators = candidates > height_difference
height_difference[indicators] = candidates[indicators]
direction_matrix[indicators] = 8

# west
candidates = minus_one_canvas
candidates[:, 1:columns:1] = height_matrix[:, 1:columns:1] - height_matrix[:, 0:columns-1:1]
indicators = candidates > height_difference
height_difference[indicators] = candidates[indicators]
direction_matrix[indicators] = 16

# north-west
candidates = minus_one_canvas
candidates[1:rows:1, 1:columns:1] = height_matrix[1:rows:1, 1:columns:1] - height_matrix[0:rows-1:1, 0:columns-1:1]
indicators = candidates > height_difference
height_difference[indicators] = candidates[indicators]
direction_matrix[indicators] = 32

# north
candidates = minus_one_canvas
candidates[1:rows:1, 1:columns:1] = height_matrix[1:rows:1, 1:columns:1] - height_matrix[0:rows-1:1, 0:columns-1:1]
indicators = candidates > height_difference
height_difference[indicators] = candidates[indicators]
direction_matrix[indicators] = 64

# north-east
candidates = minus_one_canvas
candidates[1:rows:1, 0:columns-1:1] = height_matrix[1:rows:1, 0:columns-1:1] - height_matrix[0:rows-1:1, 1:columns:1]
indicators = candidates > height_difference
height_difference[indicators] = candidates[indicators]
direction_matrix[indicators] = 128

height_matrix = None  # no longer required
# calculate the slopes

# distances of the cell with directions
distances = cell_dimensions(height_raster_location)  # function that gives cell dimensions using a reference raster

print('bbbbbbbbbbbbbbbbbbbbbbbbbb')

# to speed up runtime, calculate the inverse of the distance directly. This allows for multiplication instead of
# division when calculating the slopes.

# scale_factor multiplies the slopes such that they can be stored as integers with sufficient precision
horizontal_distance_inv = scale_factor / distances[0]  # distances[0] contains horizontal distances
vertical_distance_inv = scale_factor / distances[1]
diagonal_distance_inv = scale_factor / distances[2]

# convert the vectors to integers as multiplication with integers is much faster than with floats.
horizontal_distance_inv = horizontal_distance_inv.astype(numpy.int16)
vertical_distance_inv = vertical_distance_inv.astype(numpy.int16)
diagonal_distance_inv = diagonal_distance_inv.astype(numpy.int16)

inv_distance_matrix = numpy.zeros([rows, columns], dtype=numpy.int16)  # reserve memory for the distance matrix

# calculates the distances of those cells whose direction is horizontal. All other cells keep values 0.
inv_distance_matrix = ((direction_matrix == 1) | (direction_matrix == 16)) * horizontal_distance_inv

# adds the distances of all those cells whose direction is vertical. Only cells that have diagonal direction remain 0.
inv_distance_matrix += ((direction_matrix == 4) | (direction_matrix == 64)) * vertical_distance_inv
diagonal_directions = ((direction_matrix == 2) | (direction_matrix == 8) | (direction_matrix == 32) |
                       (direction_matrix == 128))
inv_distance_matrix += diagonal_directions * diagonal_distance_inv
slopes = height_difference * inv_distance_matrix  # finally, calculate slopes as height difference over distance.

# Downscaling the raster (e.g. from 15 seconds to 3 seconds)
# the downscale parameter determines the ratio between the output resolution and the input resolution.
new_rows = int(rows/downscale_factor)
new_columns = int(columns/downscale_factor)
# the following function packs the matrix into a smaller matrix whose elements have dimensions downscale_factor x
# downscale_factor. E.g. slopes_15_s[i, :, j, :] has dimensions downscale by downscale.
#
slopes = slopes[0:downscale_factor*new_rows, 0:downscale_factor*new_columns]
slopes_15_s = slopes.reshape([new_rows, downscale_factor, new_columns, downscale_factor])
slopes = None
slopes_15_s = slopes_15_s.mean(3, dtype=numpy.int16).mean(1, dtype=numpy.int16)  # downscale the resolution by taking the mean of the elements
slopes_15_s[slopes_15_s < minimum_slope] = minimum_slope


print('ccccccccccccccccccccc')


# Rasterize
output = gdal.GetDriverByName('GTiff').Create(slopes_15_location, new_columns, new_rows, 1, gdal.GDT_Int16,
                                              options=['COMPRESS=DEFLATE'])
output.SetProjection(height_raster.GetProjectionRef())  # gives the raster the same projection as the height raster
transform = list(height_raster.GetGeoTransform())
transform[1] = transform[1] * downscale_factor  # adjusts the pixel sizes to the new resolution
transform[5] = transform[5] * downscale_factor  # adjusts the pixel sizes to the new resolution
output.SetGeoTransform(transform)  # sets the pixel sizes and the location of the raster on a map.
output.GetRasterBand(1).WriteArray(slopes_15_s)  # writes the array data into the new raster.

# create reference raster
direction_raster_upscale = gdal.Open(direction_raster_upscale_location)
output = gdal.GetDriverByName('GTiff').Create(reference_raster_location, new_columns, new_rows, 1, gdal.GDT_Byte,
                                              options=['COMPRESS=DEFLATE'])
output.SetProjection(direction_raster_upscale.GetProjectionRef())  # gives the raster the same projection as the height raster
output.SetGeoTransform(direction_raster_upscale.GetGeoTransform())  # sets the pixel sizes and the location of the raster on a map.
height_raster = None

"""
#os.remove(direction_raster_location)
#os.remove(height_raster_location)
"""
