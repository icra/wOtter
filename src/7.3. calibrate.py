import networkx
import numpy
import pandas
import pickle
import os
import time
from return_contaminant import simulated_contaminants, error_formulae
from scipy.optimize import minimize

from src.library import graph_functions


# Load data
directory = os.path.join(os.getcwd(), 'data')
directory2 = os.path.join(os.getcwd(), 'results')

# outputs
last_location = os.path.join(directory2, 'last.npy')

def objective_function(params, fixed):
    k, s_eff = params
    beta_0 = 1
    f_eff = 1
    t_eff = 0.8

    gr, rl, cont, sc_nr, locs, obs, last_location = fixed
    sim_results, discharges = simulated_contaminants(f_eff, s_eff, t_eff, k, beta_0, gr, rl, cont, sc_nr,
                                                                   locs)
    excretion = numpy.sum(sim_results*obs*discharges)/numpy.sum(sim_results**2 * discharges)
    sim_results *= excretion
    error, mean_error = error_formulae(obs, sim_results, discharges, option=0, weighted=1)
    numpy.save(last_location, params)
    print(error/mean_error)
    print(params)
    if t_eff < s_eff:
        error *= 1000 * (s_eff - t_eff)
    return error


# open the input files
graph_location = os.path.join(directory, "river_graph.pkl")
reference_raster_location = os.path.join(directory, "reference_raster.tif")
contamination_df_location = os.path.join(directory, "AGG_WWTP_df.csv")
observed_df_location = os.path.join(directory, "pollution_observed.csv")
scenario_number = ""

contamination_df = pandas.read_csv(contamination_df_location)
observed_df = pandas.read_csv(observed_df_location)

datapoint_locations = observed_df["locations"]
observed_values = observed_df["contaminant"].to_numpy()
datapoint_count = len(observed_values)

RT = "RT_" + scenario_number
dis = scenario_number

if scenario_number == "":
    RT = "RT_HR"
    dis = "flow_HR"

river_graph = graph_functions.load_selected_attributes_graph(graph_location, [RT, dis, 'basin'])
bnds = ((0, 0.05), (0.4, 1))
starting_param = [0.0056, 0.5]

# cut river graph
basins = set([])
for location in datapoint_locations:
    basins.add(river_graph.nodes[location]['basin'])

node_list = []
for node in river_graph:
    if river_graph.nodes[node]['basin'] in basins:
        node_list.append(node)

river_graph = networkx.subgraph(river_graph, node_list)
node_df = pandas.DataFrame(node_list)
contamination_df = pandas.merge(contamination_df, node_df, left_on='pixel_number', right_on=0)
node_df = node_list = None
sorted_river_list = list(networkx.topological_sort(river_graph))

res = minimize(objective_function, numpy.array(starting_param),
           args=([river_graph, sorted_river_list, contamination_df, scenario_number, datapoint_locations,
                  observed_values, last_location]), method="Nelder-Mead", bounds=bnds)
