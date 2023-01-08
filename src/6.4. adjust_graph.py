import pickle
import pandas
import os

from src.library import graph_functions

directory = os.path.join(os.getcwd(), 'data')

graph_location = os.path.join(directory, "river_graph.pkl")
topological_sort_location = os.path.join(directory, "sorted_river_list.pkl")
adjustments_location = os.path.join(directory, "Raw data/water_fix.csv")
contamination_df_location = os.path.join(directory, "AGG_WWTP_df.csv")

open_graph = open(graph_location, "rb")
river_graph = pickle.load(open_graph)
open_graph.close()
open_ts = open(topological_sort_location, "rb")
sorted_river_list = pickle.load(open_ts)

contamination_df = pandas.read_csv(contamination_df_location)
adjustments_df = pandas.read_csv(adjustments_location)
contamination_df = pandas.merge(contamination_df, adjustments_df, left_on='dcpLongitu', right_on='longitude')
graph_functions.simulate_waste_water(river_graph, contamination_df, sorted_river_list)

# saving the output
open_file = open(graph_location, "wb")
pickle.dump(river_graph, open_file)
open_file.close()
