import contextlib
import shutil
from osgeo import ogr
import pandas as pd
import numpy as np
import pandas
import os


directory = os.path.join(os.path.dirname(os.getcwd()), 'data')

# load the inputs
agglomerations_location = os.path.join(directory, "Raw data/AGG.csv")
WWTP_location = os.path.join(directory, "Raw data/WWTPS.csv")

# open the data in panda dataframes
agglomerations_df = pandas.read_csv(agglomerations_location)
WWTP_df = pandas.read_csv(WWTP_location)  # used for getting consistent country codes

# outputs
WWTP_output_location = os.path.join(directory, "WWTP.shp")
AGG_output_location = os.path.join(directory, "AGG.shp")
country_id_table_location = os.path.join(directory, "Country id equivalence table.csv")

# add country ids

# obtain the first 2 characters of the WWTP id to have the country
for i in range(len(WWTP_df)):
    WWTP_df['dcpCode'][i] = str(WWTP_df['dcpCode'][i])[0:2]

# obtain the unique countries and build a dataframe with its equivalent cardinal
unique_id_2 = list(np.unique(WWTP_df["dcpCode"]))
unique_num = [i+1 for i in range(len(unique_id_2))]  # creates an id for each country code

# Save correspondence between id and country code
countryID_df = pd.DataFrame({'Country': unique_id_2,  'id_Country': unique_num})
countryID_df.to_csv(country_id_table_location, index=False)  # Save the country number equivalence

for i in range(len(agglomerations_df)):
    agglomerations_df['aggCode'][i] = str(agglomerations_df['aggCode'][i])[0:2]

# add the respective country id to the main dataset
agglomerations_df = pd.merge(countryID_df, agglomerations_df, left_on='Country', right_on='aggCode', how='outer')
WWTP_df = pd.merge(countryID_df, WWTP_df, left_on='Country', right_on='dcpCode', how='outer')

# remove unuseful columns
agglomerations_df.pop('aggCode')
agglomerations_df.pop('Country')
WWTP_df.pop('dcpCode')
WWTP_df.pop('Country')

# create the point shapefile
driver = ogr.GetDriverByName('Esri Shapefile')
ds = driver.CreateDataSource(WWTP_output_location)
layer = ds.CreateLayer('', None, ogr.wkbPoint)

# create fields. Fields contain a maximum of 10 characters, hence why some words are truncated.
layer.CreateField(ogr.FieldDefn('countryID', ogr.OFTInteger))
layer.CreateField(ogr.FieldDefn('treated PE', ogr.OFTInteger))
layer.CreateField(ogr.FieldDefn('Primary', ogr.OFTInteger))
layer.CreateField(ogr.FieldDefn('Secondary', ogr.OFTInteger))
layer.CreateField(ogr.FieldDefn('Other', ogr.OFTInteger))
layer.CreateField(ogr.FieldDefn('NRemoval', ogr.OFTInteger))
layer.CreateField(ogr.FieldDefn('PRemoval', ogr.OFTInteger))
layer.CreateField(ogr.FieldDefn('UV', ogr.OFTInteger))
layer.CreateField(ogr.FieldDefn('Chlorination', ogr.OFTInteger))
layer.CreateField(ogr.FieldDefn('Ozonation', ogr.OFTInteger))
layer.CreateField(ogr.FieldDefn('Sand', ogr.OFTInteger))
layer.CreateField(ogr.FieldDefn('MicroFiltr', ogr.OFTInteger))
layer.CreateField(ogr.FieldDefn('uwwOther', ogr.OFTInteger))
layer.CreateField(ogr.FieldDefn('Specification', ogr.OFTString))
layer.CreateField(ogr.FieldDefn('dcpLatitude', ogr.OFTReal))
layer.CreateField(ogr.FieldDefn('dcpLongitude', ogr.OFTReal))


defn = layer.GetLayerDefn()  # stores the type of features this layer has

errorList = []  # list to keep track of erroneous points
error = 0  # count of erroneous points

# defining a shapefile element for each treatment plant
for i in range(len(WWTP_df)):  # for all treatment plants.
    try:
        # store location
        latitude = float(WWTP_df.iat[i, 1])
        longitude = float(WWTP_df.iat[i, 2])

        if float(WWTP_df.iat[i, 3]) < 0:  # treated persons must be positive
            raise ValueError

    except ValueError:  # occurs if invalid latitude or longitude or if negative treatment
        errorList.append([WWTP_df.iat[i, 1], WWTP_df.iat[i, 2]])
        error += 1
        continue  # continue from the top of the for loop

    feature = ogr.Feature(defn)  # create the feature
    # next, link the point with the feature
    point = ogr.Geometry(ogr.wkbPoint)
    point.AddPoint(longitude, latitude)
    feature.SetGeometry(point)

    # Store the information of the treatment plant in the new feature
    feature.SetField('countryID', int(WWTP_df.iat[i, 0]))
    feature.SetField('treated PE', float(WWTP_df.iat[i, 3]))
    feature.SetField('Primary', float(WWTP_df.iat[i, 4]))
    feature.SetField('Secondary', float(WWTP_df.iat[i, 5]))
    feature.SetField('Other', float(WWTP_df.iat[i, 6]))
    feature.SetField('NRemoval', float(WWTP_df.iat[i, 7]))
    feature.SetField('PRemoval', float(WWTP_df.iat[i, 8]))
    feature.SetField('UV', float(WWTP_df.iat[i, 9]))
    feature.SetField('Chlorinati', float(WWTP_df.iat[i, 10]))
    feature.SetField('Ozonation', float(WWTP_df.iat[i, 11]))
    feature.SetField('Sand', float(WWTP_df.iat[i, 12]))
    feature.SetField('MicroFiltr', float(WWTP_df.iat[i, 13]))
    feature.SetField('uwwOther', float(WWTP_df.iat[i, 14]))
    feature.SetField('Specificat', str(WWTP_df.iat[i, 15]))
    feature.SetField('dcpLatitud', float(WWTP_df.iat[i, 16]))
    feature.SetField('dcpLongitu', float(WWTP_df.iat[i, 17]))

    layer.CreateFeature(feature)  # add the feature to the layer
    feature = point = None  # empty the feature and point for next iteration

ds = layer = None  # saves the layer and the shapefile

# continue by adding agglomerations
# create the point shapefile
driver = ogr.GetDriverByName('Esri Shapefile')
ds = driver.CreateDataSource(AGG_output_location)
layer = ds.CreateLayer('', None, ogr.wkbPoint)

# create fields. Fields contain a maximum of 10 characters, hence why some words are truncated.
layer.CreateField(ogr.FieldDefn('countryID', ogr.OFTInteger))
layer.CreateField(ogr.FieldDefn('treated PE', ogr.OFTInteger))
layer.CreateField(ogr.FieldDefn('Primary', ogr.OFTInteger))
layer.CreateField(ogr.FieldDefn('Secondary', ogr.OFTInteger))
layer.CreateField(ogr.FieldDefn('Other', ogr.OFTInteger))
layer.CreateField(ogr.FieldDefn('NRemoval', ogr.OFTInteger))
layer.CreateField(ogr.FieldDefn('PRemoval', ogr.OFTInteger))
layer.CreateField(ogr.FieldDefn('UV', ogr.OFTInteger))
layer.CreateField(ogr.FieldDefn('Chlorination', ogr.OFTInteger))
layer.CreateField(ogr.FieldDefn('Ozonation', ogr.OFTInteger))
layer.CreateField(ogr.FieldDefn('Sand', ogr.OFTInteger))
layer.CreateField(ogr.FieldDefn('MicroFiltr', ogr.OFTInteger))
layer.CreateField(ogr.FieldDefn('uwwOther', ogr.OFTInteger))
layer.CreateField(ogr.FieldDefn('Specification', ogr.OFTString))
layer.CreateField(ogr.FieldDefn('dcpLatitude', ogr.OFTReal))
layer.CreateField(ogr.FieldDefn('dcpLongitude', ogr.OFTReal))


defn = layer.GetLayerDefn()  # stores the type of features this layer has

# add all agglomerations
for i in range(len(agglomerations_df)):
    try:
        # store location
        latitude = float(agglomerations_df.iat[i, 1])
        longitude = float(agglomerations_df.iat[i, 2])

        feature = ogr.Feature(defn)  # create the feature if a location was found
    except ValueError:  # else, add to the error list for a null point
        errorList.append([agglomerations_df.iat[i, 1], agglomerations_df.iat[i, 2]])
        error += 1
        continue  # start from the top of the for loop

    # next, set the location of the shapefile element
    point = ogr.Geometry(ogr.wkbPoint)
    point.AddPoint(longitude, latitude)
    feature.SetGeometry(point)

    # convert string into a unique identifier for the country

    # Store the information of the treatment plant in the new feature
    feature.SetField('countryID', int(agglomerations_df.iat[i, 0]))
    feature.SetField('treated PE', 0)
    feature.SetField('Primary', 0)
    feature.SetField('Secondary', 0)
    feature.SetField('Other', 0)
    feature.SetField('NRemoval', 0)
    feature.SetField('PRemoval', 0)
    feature.SetField('UV', 0)
    feature.SetField('Chlorinati', 0)
    feature.SetField('Ozonation', 0)
    feature.SetField('Sand', 0)
    feature.SetField('MicroFiltr', 0)
    feature.SetField('uwwOther', 0)
    feature.SetField('Specificat', 0)
    feature.SetField('dcpLatitud', latitude)
    feature.SetField('dcpLongitu', longitude)

    layer.CreateFeature(feature)  # add the feature to the layer
    feature = point = None  # empty the feature and point for next iteration

ds = layer = None