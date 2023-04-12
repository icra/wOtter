import os
import pickle
import networkx
import pandas
import numpy
from time import time
from src.library import graph_functions
from src.library import shapefile_raster_functions
from src.library import matrix_functions
from src.library.graph_functions import simple_load
import sys
from config.config import DATA_DIR


directory = DATA_DIR

#set current directory as working directory
os.chdir(os.path.dirname(os.path.abspath(__file__)))


# input_files
graph_location = os.path.join(directory, "river_graph.pkl")
topological_sort_location = os.path.join(directory, "sorted_river_list.pkl")
reference_raster_location = os.path.join(directory, "reference_raster.tif")
contamination_df_location = os.path.join(directory, "AGG_WWTP_df.csv")
ordered_basins_location = os.path.join(directory, "ordered_basins.pkl")
wwtp_per_basin_location = os.path.join(directory, 'WWTP_per_basin.pkl')
basin_matrices_location = os.path.join(directory, "basin_matrices.pkl")

parameters = [1.927644446155748, 0.00995, 1, 0.3, 0.7, 0.9]
iterations = 5

runtimes = []
# overhead
current_time = time()
river_graph = graph_functions.load_selected_attributes_graph(graph_location, ['flow_HR', 'RT_HR'])
ts = list(networkx.topological_sort(river_graph))
cont_df = pandas.read_csv(contamination_df_location)
runtimes.append(time() - current_time)
# full implementation
current_time = time()
for i in range(iterations):
    river_graph = graph_functions.run_model(river_graph, ts, cont_df, parameters)
runtimes.append((time() - current_time)/iterations)

# partial implementation

# overhead
current_time = time()
ts = simple_load(topological_sort_location)
cont_df = pandas.read_csv(contamination_df_location)
basin_nodes, basin_ids = simple_load(ordered_basins_location)
WWTP_per_basin = simple_load(wwtp_per_basin_location)
river_graph = graph_functions.load_selected_attributes_graph(graph_location, ['flow_HR', 'RT_HR'])

basin_index = basin_ids.index(9722)  # 9722 is the danube
sub_graph = networkx.subgraph(river_graph, basin_nodes[basin_index])
river_graph = None
runtimes.append(time() - current_time)

cont_basin_df = pandas.merge(pandas.DataFrame(WWTP_per_basin[basin_index]), cont_df, left_on=0, right_on='pixel_number')
current_time = time()
for i in range(iterations):
    sub_graph = graph_functions.run_model(sub_graph, basin_nodes[basin_index], cont_basin_df, parameters)
runtimes.append((time() - current_time)/iterations)

# matrix implementation
sub_graph = None
# overhead
current_time = time()
basin_nodes, basin_ids = simple_load(ordered_basins_location)
WWTP_per_basin = simple_load(wwtp_per_basin_location)
cont_df = pandas.read_csv(contamination_df_location)

basin_index = basin_ids.index(9722)
basin_matrices, pixels_df = matrix_functions.matrix_subset_sparse(
    basin_matrices_location, [9722], basin_nodes[basin_index], basin_ids, cut_size=5000)
basin_matrices = matrix_functions.create_attenuation_matrices_sparse(basin_matrices, parameters[1])
ts = simple_load(topological_sort_location)
cont_basin_df = pandas.merge(pandas.DataFrame(WWTP_per_basin[basin_index]), cont_df, left_on=0, right_on='pixel_number')
pixels_df['init_cont'] = 0
for row in range(len(cont_basin_df)):
    pixels_df['init_cont'].loc[cont_basin_df['pixel_number'].iloc[row]] += cont_basin_df['pollution'].iloc[row]
initial_contaminant = numpy.array(pixels_df['init_cont'])
runtimes.append(time() - current_time)

# run
current_time = time()
for i in range(iterations):
    matrix_functions.run_basin_matrices_sparse(basin_matrices, basin_nodes[basin_index], initial_contaminant)
runtimes.append((time() - current_time)/iterations)

runtimes_df = pandas.DataFrame(runtimes)
runtimes_df.index = ['overhead_full', 'runtime_full', 'overhead_partial', 'runtime_partial', 'overhead_matrix',
                     'runtime_matrix']
runtimes_df.to_csv('runtimes.csv')