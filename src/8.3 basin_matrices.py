import os
import pickle
import pandas

from src.library import matrix_functions


directory = os.path.join(os.getcwd(), 'data')

# input files
graph_location = os.path.join(directory, 'river_graph.pkl')
ordered_basins_location = os.path.join(directory, "ordered_basins.pkl")

# output files
basin_matrices_location = os.path.join(directory, "basin_matrices.pkl")

open_graph = open(graph_location, "rb")
river_graph = pickle.load(open_graph)
open_graph.close()

open_file = open(ordered_basins_location, "rb")
basin_list, basin_ids = pickle.load(open_file)
open_file.close()

matrix_functions.graph_to_RT_matrix(river_graph, basin_list, cut_size=800, output_location=basin_matrices_location)
