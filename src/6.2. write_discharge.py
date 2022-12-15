import geopandas
from osgeo import gdal, ogr
import numpy
import os

from src.library import shapefile_raster_functions



directory = os.path.join(os.path.dirname(os.getcwd()), 'data')

# input files
AGG_WWTP_location = rivers_shapefile = os.path.join(directory, "AGG_WWTP_adapted.shp")
reference_raster_location = os.path.join(directory, "reference_raster.tif")
river_raster_location = os.path.join(directory, "rivers_from_graph.tif")
land_shapefile_location = os.path.join(directory, "countries/Countries.shp")
land_raster_location = os.path.join(directory, "land.tif")
river_land_temp_raster_location = os.path.join(directory, "river_land.tif")


# temporary locations
temporary_AGG_WTTP_raster_location = os.path.join(directory, "AGG_WWTP.tif")
temporary_AGG_WWTP_location = os.path.join(directory, "AGG_WWTP_temp.shp")

# output location
output_location = os.path.join(directory, "AGG_WWTP_df_no_treatment.csv")

# create a raster that is land
land_raster = shapefile_raster_functions.shapefile_to_raster(land_shapefile_location, reference_raster_location, land_raster_location)

# open the shapefile of the contaminants
AGG_WWTP = ogr.Open(AGG_WWTP_location, 1)  # the 1 signifies editing mode
AGG_WWTP_layer = AGG_WWTP.GetLayer(0)  # the object only has a single layer
AGG_WWTP_layer.CreateField(ogr.FieldDefn('lat_pixel', ogr.OFTInteger64))  # stores the pixel number latitude
AGG_WWTP_layer.CreateField(ogr.FieldDefn('long_pixel', ogr.OFTInteger64))

# prepare indicator raster for the contaminant shapefile
reference_raster = gdal.Open(reference_raster_location)
point_matrix = numpy.zeros([reference_raster.RasterYSize, reference_raster.RasterXSize])

# In this step the shapefile is rasterized
rows, columns = numpy.shape(point_matrix)
for feature in AGG_WWTP_layer:
    pt = feature.geometry()  # returns a geometry object
    if pt is not None:  # do not consider empty geometry objects
        try:
            pt_x = float(feature.GetField('dcpLatitud'))
            pt_y = float(feature.GetField('dcpLongitu'))
        except TypeError:
            pt_x = None

        if pt_x is None or pt_y is None:
            pt_x = pt.GetY()
            pt_y = pt.GetX()
        # give_pixel converts the coordinates to the correct raster pixel numbers
        x_location, y_location = shapefile_raster_functions.give_pixel([pt_x, pt_y], reference_raster)
        # save the location of the shapefile geometry in the raster
        feature.SetField('lat_pixel', x_location)
        feature.SetField('long_pixel', y_location)
        AGG_WWTP_layer.SetFeature(feature)
        if 0 <= y_location <= columns:
            if 0 <= x_location <= rows:
                point_matrix[x_location, y_location] = 1  # rasterize the shapefile geometry
    else:
        AGG_WWTP_layer.DeleteFeature(feature.GetFID())  # if geometry object is empty

# save the indicator matrix as a raster such that it can be used in some standard functions
output = gdal.GetDriverByName('GTiff').Create(temporary_AGG_WTTP_raster_location, reference_raster.RasterXSize,
                                              reference_raster.RasterYSize, 1, gdal.GDT_Byte)
output.SetProjection(reference_raster.GetProjectionRef())
output.SetGeoTransform(reference_raster.GetGeoTransform())
output.GetRasterBand(1).WriteArray(point_matrix)
output = None

# create a raster of valid discharge points. A valid discharge point is river or ocean.
river_raster = gdal.Open(river_raster_location)
river_raster_matrix = river_raster.GetRasterBand(1).ReadAsArray()
ocean_raster_matrix = 1 - land_raster.GetRasterBand(1).ReadAsArray()  # e.g. value is 1 if ocean and 0 if land
land_raster = None
# gives an indicator raster for either ocean or river
ocean_and_river_matrix = river_raster_matrix + ocean_raster_matrix - river_raster_matrix * ocean_raster_matrix
river_raster = None
# Save the ocean or river matrix as a raster.
output = gdal.GetDriverByName('GTiff').Create(river_land_temp_raster_location, reference_raster.RasterXSize,
                                              reference_raster.RasterYSize, 1, gdal.GDT_Byte)
output.SetProjection(reference_raster.GetProjectionRef())
output.SetGeoTransform(reference_raster.GetGeoTransform())
output.GetRasterBand(1).WriteArray(ocean_and_river_matrix)
output = None

# the function below finds the closest discharge point for the AGG WWTP raster in the river_ocean raster.
discharge_matrix = shapefile_raster_functions.find_discharge(temporary_AGG_WTTP_raster_location, river_land_temp_raster_location,
                                     maximum_radius=18, precision=1)
discharge_matrix = discharge_matrix.transpose()
ocean_raster_matrix = numpy.transpose(ocean_raster_matrix)

# for each feature record the discharge point.
for feature in AGG_WWTP_layer:
    x_location = feature.GetField('lat_pixel')
    y_location = feature.GetField('long_pixel')
    point_in_raster = False
    if 0 <= x_location <= rows:
        if 0 <= y_location <= columns:
            point_in_raster = True
            discharge_point = discharge_matrix[x_location, y_location]  # contains the discharge point location

            # if a discharge point is found and it is not in the ocean:
            if discharge_point != ['none found'] and\
                    int(ocean_raster_matrix[discharge_point[0], discharge_point[1]]) == 0:
                # then adapt the feature location to the discharge location
                feature.SetField('long_pixel', discharge_point[0])
                feature.SetField('lat_pixel', discharge_point[1])
                AGG_WWTP_layer.SetFeature(feature)
            else:  # ignore if no discharge point is found, or if it is in the ocean (then it does not affect the model)
                AGG_WWTP_layer.DeleteFeature(feature.GetFID())

    if not point_in_raster:  # delete if point falls outside of the raster
        AGG_WWTP_layer.DeleteFeature(feature.GetFID())
# save file and empty memory
AGG_WWTP_layer = None
AGG_WWTP = None

# open the shapefile as a geopandas data frame
AGG_WWTP_df = geopandas.read_file(AGG_WWTP_location)

# ensure that the treatment fields are dummies
treatment_fields = ['Primary', 'Secondary', 'Other', 'NRemoval', 'PRemoval', 'UV', 'Chlorinati', 'Ozonation', 'Sand',
                    'MicroFiltr', 'uwwOther']
for field in treatment_fields:
    AGG_WWTP_df[field] = AGG_WWTP_df[field] == 1  # replaces non valid values with 0
    AGG_WWTP_df[field] = AGG_WWTP_df[field].astype(int)

# pixel number determines in which point of the river graph there is discharge.
AGG_WWTP_df["pixel_number"] = AGG_WWTP_df["long_pixel"] + AGG_WWTP_df["lat_pixel"] * columns

# fill nans and other invalid values
for i in range(len(AGG_WWTP_df)):
    if not AGG_WWTP_df['Filt_a'][i] > 0:
        AGG_WWTP_df['Filt_a'][i] = 0
AGG_WWTP_df = AGG_WWTP_df.fillna(0)

# save to disk
other_fields = ['countryID', 'dcpLatitud', 'dcpLongitu', 'Treat_a', 'Filt_a', 'Unfilt_a', 'pollution', 'pixel_number']
AGG_WWTP_df = AGG_WWTP_df[treatment_fields + other_fields]  # extracts only specified fields
AGG_WWTP_df.to_csv(output_location)

# remove temporary files
os.remove(temporary_AGG_WTTP_raster_location)
os.remove(river_land_temp_raster_location)
os.remove(land_raster_location)
