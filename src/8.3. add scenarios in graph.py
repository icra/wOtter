import pickle
from time import time
import os
import geopandas as gpd

from src.library import graph_functions


# Load data
directory = os.path.join(os.getcwd(), 'data')

# input
river_graph_location = os.path.join(directory, "river_graph.pkl")

# parameters
# model = 'csc_mpi'  # csc_mpi, rac_hadgem, rca4_hadgem, rca4_mpi
rivers_shapefile = os.path.join(directory, "Raw data", "Copernicus scenarios shp", "HR and csc_mpi scenarios.shp")
all = 0  # 0 is yes, 1 if no
scen_names = ['flow1', 'flow2', 'flow3', 'flow5', 'flow6', 'flow8', 'flow9']

# output
graph_scenarios_location = os.path.join(directory, "graph with scenarios csc_mpi 1 2 3 5 6 8 9.pkl")


starting_time = time()

# Load river shapefile with scenarios.

# model = input("Write the scenario model name you want to add (csc_mpi, rac_hadgem, rca4_hadgem, rca4_mpi).")
# if model == "csc_mpi": rivers_shapefile = os.path.join(directory, "Copernicus scenarios shp", "HR and csc_mpi scenarios.shp")
# if model == "rac_hadgem": rivers_shapefile = os.path.join(directory, "Copernicus scenarios shp", "HR and rac_hadgem scenarios.shp")
# if model == "rca4_hadgem": rivers_shapefile = os.path.join(directory, "Copernicus scenarios shp", "HR and rca4_hadgem scenarios.shp")
# if model == "rca4_mpi": rivers_shapefile = os.path.join(directory, "Copernicus scenarios shp", "HR and rca4_mpi scenarios.shp")


rd_time0 = time()
shp_rivers = gpd.read_file(rivers_shapefile)
rd_timeF = time()

print("Time river shapefile reading (min): {}".format((rd_timeF - rd_time0)/60))

# Load graph and topological sort data
try:
    rd_time0 = time()

    open_graph = open(os.path.join(directory, river_graph_location), "rb")
    loaded_graph = pickle.load(open_graph)
    open_graph.close()

    river_graph = loaded_graph[0]

    rd_timeF = time()
    print("Time graph reading (min): {}".format((rd_timeF - rd_time0)/60))

except FileNotFoundError:
    print("Graph is not created.")
    loaded_graph = False

##############################################################################################################
# Add scenarios
##############################################################################################################

if loaded_graph:
    # get inputs number of scenarios and their names.
    # all = int(input("Do you want to add all the scenarios? Write 0 if the answer is not, write 1 if yes."))
    # if all == 0:
    #     scen_names = []
    #     num_scenarios = int(input("Enter the number of scenarios that you want to add (e.g. 3)."))
    #     for i in range(num_scenarios):
    #         new_name = input("Write the new scenario name (e.g. flow1).")
    #         scen_names.append(new_name)
    #     print("The new scenarios added will be:")
    #     print(scen_names)
    # else:
    #     scen_names = []

    print("Loading scenarios...")

    graph_functions.add_scenario(river_graph, shp_rivers, scen_names=scen_names, all_scenarios=all)

    print("Copernicus scenarios with RT added in graph.")

    open_file = open(os.path.join(directory, graph_scenarios_location), "wb")
    pickle.dump(river_graph, open_file)
    open_file.close()

    print("Graph has been saved in graph with scenarios.pkl file.")
