import os
import graph_functions
import pickle
import pandas

directory = os.path.join(os.path.dirname(os.getcwd()), 'data')

# input files
graph_location = os.path.join(directory, 'river_graph.pkl')
cont_df_location = os.path.join(directory, 'AGG_WWTP_df.csv')

# output files
ordered_basins_location = os.path.join(directory, "ordered_basins.pkl")
wwtp_per_basin_location = os.path.join(directory, 'WWTP_per_basin.pkl')

river_graph = graph_functions.load_selected_attributes_graph(graph_location, ['basin'])

basins_list, basin_ids = graph_functions.create_basin_lists(river_graph, ordered_basins_location)
river_graph = None
cont_df = pandas.read_csv(cont_df_location)
graph_functions.discharge_per_basin(basins_list, cont_df, wwtp_per_basin_location)
