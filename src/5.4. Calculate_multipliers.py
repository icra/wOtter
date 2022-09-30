from fiona.errors import DriverIOError
import geopandas
import numpy
import pandas
import statsmodels.api as sm
import sys
import os

def merge_rows(array, row_index, basin_list, depth):
    """"Function takes the current row and combines it with a basin in the same parent-basin with the lowest
     population"""
    n = numpy.shape(array)[0]  # amount of rows
    array[:, 6] = numpy.array([i for i in range(n)])
    current_set = basin_list[row_index]  # will store what sub-basins are included in a new row
    array_id = array[:, 0]  # these are Pfafstetter codes
    found_basin = False
    count = 0
    while not found_basin:
        count += 1
        if count > depth:
            # found_basin remains false here, such that the outer function knows no basin was found.
            return array, basin_list, found_basin  # if no basin is found at the depth, return the inputs unchanged

        array_id = numpy.divide(array_id, 10).astype(int)  # increases the level of the basin resolution
        # 'indices' finds sub-basins that are part of the same basin at a larger resolution. It also requires that the
        # basins are in the same country.
        indices = numpy.where((array_id == array_id[row_index]) & (array[:, 5] ==
                                                                   array[row_index, 5]), 1, 0).astype(bool)
        indices[row_index] = False  # Basin cannot merge with itself
        found_basin = (sum(indices) > 0)

    # this section can only be entered if a mergeable basin was found
    eligible_basins = array[indices, :]

    # find eligible basin with least difference in the Pfafstetter codes
    eligible_basins[:, 0] = eligible_basins[:, 0] - array[row_index, 0]
    eligible_basins[:, 0] = numpy.abs(eligible_basins[:, 0])
    minimum = numpy.amin(eligible_basins[:, 0])
    minimum_row = numpy.where(eligible_basins[:, 0] == minimum)
    minimum_row_number = int(eligible_basins[minimum_row[0][0], 6])  # row of the picked basin in the original array

    # update the list of sub-basins that are included in the lower resolution basin.
    [current_set.append(item) for item in basin_list[minimum_row_number]]

    # create the newly merged basin
    new_row = array[row_index, :] + array[minimum_row_number, :]  # add the rows
    new_row[0] = array[minimum_row_number, 0]  # copy the Pfafstetter code of the merged basin
    new_row[5] = array[minimum_row_number, 5]  # set the country id as that of the merged basin
    new_row = new_row.reshape([1, len(new_row)])

    # delete the two rows and merge the combined row
    array = numpy.delete(array, [row_index, minimum_row_number], 0)
    array = numpy.concatenate([new_row, array], 0)

    # delete the two sub-basins from the basin list and add the newer larger basin
    if row_index > minimum_row_number:
        basin_list.pop(row_index), basin_list.pop(minimum_row_number), basin_list.insert(0, current_set)
    else:
        basin_list.pop(minimum_row_number), basin_list.pop(row_index), basin_list.insert(0, current_set)

    return array, basin_list, found_basin


directory = os.path.join(os.path.dirname(os.getcwd()), 'data')

# inputs
basins_location = os.path.join(directory, "basin_with_filt_unfilt.shp")
country_shapefile = os.path.join(directory, "countries/Countries.shp")
country_ids_location = os.path.join(directory, "Country id equivalence table.csv")

# outputs
output_shapefile_location = os.path.join(directory, "countries with multipliers/countries with multipliers.shp")

# Open and convert the data to numpy
data_basins = geopandas.read_file(basins_location)
data_basins = data_basins[['PFAF_ID', 'population', 'treated_p', 'filtered p', 'unfiltered', 'CNTR_ID']]
country_ids = pandas.read_csv(country_ids_location)
data_basins = pandas.merge(data_basins, country_ids, left_on="CNTR_ID", right_on="Country")
data_basins = data_basins.drop(columns=['CNTR_ID', 'Country'])
data_basins = data_basins.to_numpy()

# delete invalid data
data_basins = numpy.delete(data_basins, numpy.isnan(data_basins[:, 1]), axis=0)
data_basins = numpy.delete(data_basins, numpy.isnan(data_basins[:, 2]), axis=0)
data_basins = numpy.delete(data_basins, data_basins[:, 2] < 0, axis=0)

# add row number column in the data matrix. This is necessary for the 'merge-rows' function.
rows, columns = numpy.shape(data_basins)
numbers = numpy.array([i for i in range(rows)])
numbers = numbers.reshape([rows, 1])
data_basins = numpy.concatenate([data_basins, numbers], 1)

# list that contains the basins included in a row. Necessary input for the merge-rows function.
basin_set = [[element] for element in data_basins[:, 0]]

# This next piece of code uses the merge_rows function to ensure that all the basins are sufficiently large
too_small = True  # Boolean that determines whether there are still basins that are too small
too_small_value = 10000  # minimal population of an observation basin
i = -1

# depth determines what basins may be combined. If 1, level 7 basins may only be combined if they are in the same level
# 6 basin. If depth is 2, level 7 basins may be combined if they are in the same level 5 basin.
depth = 1
while too_small:  # while some basin is too small
    if depth > 100:
        sys.exit('infinite loop')
    too_small = False
    any_basins_found = False
    i = -1
    while i < rows - 1:  # go through all the basins
        i += 1
        # adapt a row if the population is too small, or treated is 0.
        if data_basins[i, 1] < too_small_value or data_basins[i, 2] == 0:
            # merge rows gives a new dataset with a merged basin, and a basin_set that shows what basins are included in
            # a row. Basins_found is a boolean that shows whether a suitable merge was found.
            data_basins, basin_set, basins_found = merge_rows(data_basins, i, basin_set, depth)
            rows, columns = numpy.shape(data_basins)
            too_small = True  # repeat the loop if a basin was found that was too small
            if basins_found:  # prevents depth from increasing
                any_basins_found = True
    if not any_basins_found:
        # if more basins need to be merged, but no suitable merges are found. Extent the possible depth for searching.
        depth += 1

# create a list that contains the country name in a list for each basin row in data_basins
country_id_vector = data_basins[:, 5]
country_names = []
for i in range(len(country_id_vector)):
    country_names.append(country_ids["Country"].iloc[int(country_id_vector[i]-1)])

# create dummies for each country
dummy_matrix = pandas.get_dummies(country_names)

# create the regressors and dependent variables
pop_treated = numpy.divide(data_basins[:, 1], data_basins[:, 2])
filtered_treated = numpy.divide(data_basins[:, 3], data_basins[:, 2])
unfiltered_treated = numpy.divide(data_basins[:, 4], data_basins[:, 2])
Y = pop_treated - filtered_treated - unfiltered_treated

# Do and print the regression
model = sm.WLS(Y, dummy_matrix, data_basins[:, 2])
result = model.fit()
print(result.summary())

# store the parameters in a new shapefile
parameters = result.params
parameters = pandas.DataFrame({'Country': parameters.index, 'treat_multiplier': parameters.values})
country_data = geopandas.read_file(country_shapefile)
country_data = pandas.merge(country_data, parameters, left_on="CNTR_ID", right_on="Country")
country_data = country_data[["Country", "CNTR_ID", "geometry", "treat_multiplier"]]
try:
    country_data.to_file(output_shapefile_location)
except (FileExistsError, DriverIOError):
    split_string = output_shapefile_location.split('/')
    os.makedirs(os.path.join(directory, split_string[0]))
    country_data.to_file(output_shapefile_location)
