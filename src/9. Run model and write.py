import math
import os
import networkx
import numpy
import pandas
import pickle
import matplotlib.pyplot as plt
import time

# Load data
directory = os.path.join(os.path.dirname(os.getcwd()), 'data')
directory2 = os.path.join(os.path.dirname(os.getcwd()), 'results')

# output
written_model_location = os.path.join(directory2, "written_model.pkl")

# open the input files
graph_location = os.path.join(directory, "river_graph.pkl")
topological_sort_location = os.path.join(directory, "sorted_river_list.pkl")
reference_raster_location = os.path.join(directory, "reference_raster.tif")
contamination_df_location = os.path.join(directory, "AGG_WWTP_df_all_adapted.csv")
scenarios = ['']

contamination_df = pandas.read_csv(contamination_df_location)

open_graph = open(graph_location, "rb")
river_graph = pickle.load(open_graph)
river_graph = river_graph[0]
open_graph.close()

open_ts = open(topological_sort_location, "rb")
sorted_river_list = pickle.load(open_ts)
sorted_river_list = sorted_river_list[0]
open_ts.close()

k_dec, filt_eff, s_eff, t_eff = [0.00995203, 0.99993389, 0.63221308, 0.81011981]
beta_0 = 5.01466473

gr, rl, cont = [river_graph, sorted_river_list, contamination_df]

for scenario_number in scenarios:
    # contamination
    treatment_efficacy = 0.3 * (contamination_df["Treatment_level"] == 1) + s_eff * \
                         (contamination_df["Treatment_level"] == 2) + t_eff * \
                         (contamination_df["Treatment_level"] == 3)
    # formula
    contamination = (1 - treatment_efficacy) * contamination_df["Treat_a"] + (1 - filt_eff) \
                    * contamination_df["Filt_a"] + contamination_df["Unfilt_a"]
    contamination *= contamination_df["pollution"] * beta_0
    location = contamination_df["pixel_number"]

    cont = "Contaminant " + scenario_number
    RT = "RT " + scenario_number
    dis = "discharge " + scenario_number
    rel_cont = "Relative contaminant " + scenario_number
    initial_cont = "Initial contaminant"

    if scenario_number == "":
        cont = "Contaminant"
        RT = "RT_HR"
        dis = "flow_HR"
        rel_cont = "Relative Contaminant"
    networkx.set_node_attributes(river_graph, 0, name=cont)
    networkx.set_node_attributes(river_graph, 0, name=initial_cont)

    for i in range(len(contamination_df)):
        pixel_number = location[i]
        river_graph.nodes[pixel_number][cont] += contamination[i]
        river_graph.nodes[pixel_number][initial_cont] = river_graph.nodes[pixel_number][cont]

    for n in sorted_river_list:  # for all river pixels
        # Add the contamination of the parent cells
        parents = list(river_graph.predecessors(n))
        for k in parents:
            river_graph.nodes[n][cont] += river_graph.nodes[k][cont]
            try:
                river_graph.nodes[n][cont] *= math.exp(-k_dec * river_graph.nodes[n][RT])
            except:
                river_graph.nodes[n][cont] = 0
        river_graph.nodes[n][rel_cont] = river_graph.nodes[n][cont] / river_graph.nodes[n][dis]

open_file = open(written_model_location, "wb")
pickle.dump(river_graph, open_file)
open_file.close()
