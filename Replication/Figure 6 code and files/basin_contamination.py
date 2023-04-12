import os
import pickle
import networkx
import numpy
import pandas
import geopandas
from config.config import DATA_DIR
from src.library import graph_functions
import sys

directory = DATA_DIR

#set current directory as working directory
os.chdir(os.path.dirname(os.path.abspath(__file__)))

# inputs
graph_location = os.path.join(directory, 'river_graph.pkl')
sorted_river_list_location = os.path.join(directory, "sorted_river_list.pkl")
reference_raster_location = os.path.join(directory, 'reference_raster.tif')
contamination_df_location = os.path.join(directory, 'AGG_WWTP_df.csv')
shapefile = os.path.join(directory, 'hybas_eu_lev06_v1c.shp')

# output
output_name = 'contamination_pfafstetter_basin.shp'
# pick the scenario
scenario = ''  # standard scenario
if scenario == '':
    RT = "RT_HR"
    dis = "flow_HR"
else:
    RT = "RT " + scenario
    dis = "discharge " + scenario
dist = "cell_distance"

river_graph = graph_functions.simple_load(graph_location)
sorted_river_list = graph_functions.simple_load(sorted_river_list_location)
region_field = 'PFAF_ID'
graph_functions.absorb_shapefile(river_graph, shapefile, [region_field], reference_raster_location)

contamination_df = pandas.read_csv(contamination_df_location)

# parameters
attenuation = 0.0045
excretion = 1
primary_efficacy = 0.33
secondary_efficacy = 0.7
tertiary_efficacy = 0.92
filtered_efficacy = 1
parameters = [excretion, attenuation, filtered_efficacy, primary_efficacy, secondary_efficacy, tertiary_efficacy]

output_field_name = "basin contamination"
river_graph = graph_functions.run_model(river_graph, sorted_river_list, contamination_df, parameters,
                                        output_field_name=output_field_name)
# delete parts of the rivers that do not have a country

for node in networkx.topological_sort(river_graph):
    if river_graph.nodes[node]['id_Country'] == 0:
        river_graph.remove_node(node)
cont = "Contamination " + output_field_name
rel_cont = "Relative contaminant " + output_field_name
pixels_by_PFAF_id = []
field_ids = []
field_ids_mean_discharge = []
for n in river_graph:
    field = int(river_graph.nodes[n][region_field])
    if field not in field_ids:
        field_ids.append(field)
        field_ids_mean_discharge.append(river_graph.nodes[n][dis])
        pixels_by_PFAF_id.append([])
    index = field_ids.index(field)
    pixels_by_PFAF_id[index].append(n)
    field_ids_mean_discharge[index] += river_graph.nodes[n][dis]

aggregator_exponent = [1.2]
discharge_weight = [0.5]

combinations = [[exponent, weight] for exponent in aggregator_exponent for weight in discharge_weight]
field_names = [str(combination[0]) +" / " + str(combination[1])
               for combination in combinations]

contamination_by_field = []

for index in range(len(field_ids)):
    field_ids_mean_discharge[index] *= 1 / len(pixels_by_PFAF_id[index])
    weight_sum = [0] * len(combinations)
    cont_sum = [0] * len(combinations)
    for pixel in pixels_by_PFAF_id[index]:
        for combination_index in range(len(combinations)):
            aggregator_exponent = combinations[combination_index][0]
            discharge_weight = combinations[combination_index][1]
            distance_weight = river_graph.nodes[pixel][dist]
            weight = distance_weight * numpy.power(river_graph.nodes[pixel][dis] / field_ids_mean_discharge[index],
                                                   discharge_weight)
            weight_sum[combination_index] += weight
            cont_sum[combination_index] += weight * numpy.power(river_graph.nodes[pixel][rel_cont], aggregator_exponent)
    normalized_contamination = list(numpy.array(cont_sum)/numpy.array(weight_sum))
    contamination_by_field.append(normalized_contamination)

for i in range(len(contamination_by_field)):
    contamination_by_field[i] = [field_ids[i]] + contamination_by_field[i]

shapefile = geopandas.read_file(shapefile)
contamination_basin_df = pandas.DataFrame(contamination_by_field)
columns = ['field_ids'] + field_names
contamination_basin_df.columns = columns

shapefile = pandas.merge(shapefile, contamination_basin_df, left_on=region_field, right_on='field_ids')
shapefile = shapefile[field_names + ['geometry']]

# normalizing each field such that 1 is the average
for i in range(len(field_names)):
    shapefile[field_names[i]] = shapefile[field_names[i]] / sum(shapefile[field_names[i]]) * len(shapefile)

shapefile.to_file(output_name)

