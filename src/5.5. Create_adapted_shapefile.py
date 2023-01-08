import geopandas
from geopandas.tools import sjoin
import numpy
import pandas
import os

directory = os.path.join(os.getcwd(), 'data')

# inputs
AGG_location = os.path.join(directory, "AGG.shp")
basin_location = os.path.join(directory, "basin_with_filt_unfilt.shp")
WWTP_location = os.path.join(directory, "WWTP.shp")
countries_location = os.path.join(directory, "countries with multipliers/countries with multipliers.shp")
country_id_equivalent_location = os.path.join(directory, "Country id equivalence table.csv")

# output
contaminant_shapefile_location = os.path.join(directory, "AGG_WWTP.shp")

# Step 1: adapting the agglomerations
point_agg = geopandas.GeoDataFrame.from_file(AGG_location)
poly = geopandas.GeoDataFrame.from_file(basin_location)
poly['unique'] = numpy.array([i for i in range(len(poly))])  # create a new id for each basin
point_in_poly = sjoin(point_agg, poly, how='left')  # creates a point dataset with data for the polygon the point is in
country_id_equivalent = pandas.read_csv(country_id_equivalent_location)
country_id_equivalent.loc[len(country_id_equivalent)] = ["RS", len(country_id_equivalent) + 1]  # Add Serbia

# new data fields
point_agg["Treat_a"] = numpy.zeros(len(point_agg))
point_agg["Filt_a"] = numpy.zeros(len(point_agg))
point_agg["Unfilt_a"] = numpy.zeros(len(point_agg))

for row in range(len(poly)):
    indices = point_in_poly["unique"] == poly["unique"][row]
    count = sum(indices)  # gives the amount of agglomerations in the basin
    indices = numpy.array(indices)
    zeros = [0] * (len(point_agg) - len(indices))  # adapts indices to the changing length of point_agg
    indices = numpy.concatenate((indices, zeros), axis=0)
    if count > 0:  # write to agglomerations
        point_agg.loc[indices.astype(bool), "Filt_a"] = poly["filtered p"][row] / count
        point_agg.loc[indices.astype(bool), "Unfilt_a"] = poly["unfiltered"][row] / count
    else:  # if no point, define a new one
        new_point = poly["geometry"][row].representative_point()
        country_id = country_id_equivalent.loc[country_id_equivalent['Country'] == poly["CNTR_ID"][row],
                                                'id_Country'].values[0]
        Filt_a = poly["filtered p"][row]
        Unfilt_a = poly["unfiltered"][row]
        new_point = [country_id, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, new_point.y,
                     new_point.x, new_point, 0, Filt_a, Unfilt_a]
        point_agg.loc[len(point_agg)] = new_point  # add new point to the shapefile
# Step 2: adapting the WWTPs

point_WWTP = geopandas.GeoDataFrame.from_file(WWTP_location)
point_WWTP = pandas.merge(point_WWTP, country_id_equivalent, left_on='countryID', right_on='id_Country')
country_multipliers = geopandas.GeoDataFrame.from_file(countries_location)
country_multipliers = country_multipliers.drop(['geometry'], 1)
WWTP_with_multipliers = pandas.merge(point_WWTP, country_multipliers, left_on='Country', right_on='Country')
WWTP_with_multipliers['Treat_a'] = WWTP_with_multipliers['treated PE'] * WWTP_with_multipliers['treat_mult']
WWTP_with_multipliers['Filt_a'] = 0
WWTP_with_multipliers['Unfilt_a'] = 0


# Joining the shapefiles
point_WWTP = geopandas.read_file(WWTP_location)
# Merge of the shapefiles of the agglomerations and WWTP
contaminant_shapefile = geopandas.GeoDataFrame(pandas.concat([WWTP_with_multipliers, point_agg]))
contaminant_shapefile['pollution'] = 1  # sets all the pollution to 1
contaminant_shapefile['country'] = 0
for row in range(len(contaminant_shapefile)):
    country_id = contaminant_shapefile['countryID'].iloc[row]
    contaminant_shapefile['country'].iloc[row] = country_id_equivalent.iloc[country_id-1, 0]
# Write the merged shapefile
contaminant_shapefile.to_file(contaminant_shapefile_location)
