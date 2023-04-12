import numpy
from osgeo import gdal
from scipy import ndimage
import os
from src.library import shapefile_raster_functions


directory = os.path.join(os.getcwd(), 'data')
# input locations
input_shapefile_location = os.path.join(directory, "WWTP.shp")
reference_raster_location = os.path.join(directory, "Raw data/country_reference_raster.tif")
countries = os.path.join(directory, "Raw data/countries original shp/CNTR_RG_01M_2020_4326.shp")

# temp location
output_WWTP_location = os.path.join(directory, "WWTP.tif")

# output locations
output_land_location = os.path.join(directory, "blurred_land.tif")
output_countries_raster = os.path.join(directory, "certain_land.tif")
output_location = os.path.join(directory, "Population_treated.tif")

# rasterize WWTP shapefile using the treated pe attribute
shapefile_raster_functions.point_shapefile_sum_to_raster(input_shapefile_location, reference_raster_location,
                                 output_name=output_WWTP_location, attributes=["treated PE"])

# rasterize the country shapefile, where 1 is land and 0 is ocean.
land_matrix = shapefile_raster_functions.shapefile_to_raster(countries, reference_raster_location, output_countries_raster, option=0)

# create a land matrix that has slightly less land
land_matrix = 1000 * numpy.where(land_matrix == 1, 1, 0)  # land tiles are 1000
land_matrix = 0.001 * ndimage.uniform_filter(land_matrix, size=3, mode='constant')  # blurs the image
land_matrix = numpy.where(land_matrix == 1, 1, 0)  # only includes non-blurred pixels, hence certain land

# write land matrix to raster
reference_raster = gdal.Open(reference_raster_location)
output = gdal.GetDriverByName('GTiff').Create(output_land_location, reference_raster.RasterXSize,
                                              reference_raster.RasterYSize, 1, gdal.GDT_Byte)
output.GetRasterBand(1).WriteArray(land_matrix)
output.SetProjection(reference_raster.GetProjection())  # copy direction projection to output raster
output.SetGeoTransform(reference_raster.GetGeoTransform())  # copy direction resolution/location to output raster
output = None

# Move the WWTP load to a 'certain' land tile, in case it is currently in an ocean tile. This is done to ensure
# it falls inside a basin when the script 'zonal statistics' is run.

# find closest land
closest_land = shapefile_raster_functions.find_discharge(output_WWTP_location, output_land_location, maximum_radius=8, precision=0.5)
closest_land = closest_land.transpose()
rows, columns = numpy.shape(closest_land)

WWTP_raster = gdal.Open(output_WWTP_location)
discharge_matrix = WWTP_raster.GetRasterBand(1).ReadAsArray()
treated_load_matrix = WWTP_raster.GetRasterBand(2).ReadAsArray()

# write to closest land
treated_load_on_land_matrix = numpy.zeros([rows, columns])
for i in range(rows):
    for j in range(columns):
        if discharge_matrix[i, j] > 0:  # if there is a treatment plant in the cell
            location = closest_land[i,j]  # take the discharge point of the treatment plant in that cell
            if location != ['none found']:  # this is the entry if no discharge point is found
                location_i = location[1]  # take the second element since closest_land was transposed
                location_j = location[0]
                treated_load_on_land_matrix[location_i, location_j] += treated_load_matrix[i, j]

# create the new WWTP raster file, where the load is on certain land.
output = gdal.GetDriverByName('GTiff').Create(output_location, WWTP_raster.RasterXSize,
                                              WWTP_raster.RasterYSize, 1, gdal.GDT_Int32)
output.GetRasterBand(1).WriteArray(treated_load_on_land_matrix)

output.SetProjection(WWTP_raster.GetProjection())  # copy direction projection to output raster
output.SetGeoTransform(WWTP_raster.GetGeoTransform())  # copy direction resolution/location to output raster
output = None

WWTP_raster = None
os.remove(output_WWTP_location)  # remove the intermediary WWTP output.
