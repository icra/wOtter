import numpy
import pandas
import pickle
from return_contaminant import simulated_contaminants, error_formulae
import matplotlib.pyplot as plt
import os
import time
from src.library import shapefile_raster_functions


# Load data
directory = os.path.join(os.path.dirname(os.getcwd()), 'data')
directory2 = os.path.join(os.path.dirname(os.getcwd()), 'results')

# output
calibrated_nonadapted_location = os.path.join(directory2, "calibrated_nonadapted.shp")

# open the input files
graph_location = os.path.join(directory, "river_graph.pkl")
topological_sort_location = os.path.join(directory, "sorted_river_list.pkl")
reference_raster_location = os.path.join(directory, "reference_raster.tif")
contamination_df_location = os.path.join(directory, "AGG_WWTP_df_all_adapted.csv")
scenario_number = ""
observed_df_location = os.path.join(directory, "pollution_observed.csv")

contamination_df = pandas.read_csv(contamination_df_location)
observed_df = pandas.read_csv(observed_df_location)

datapoint_locations = observed_df["locations"]
observed_values = observed_df["contaminant"].to_numpy()
datapoint_count = len(observed_values)
open_graph = open(graph_location, "rb")
river_graph = pickle.load(open_graph)
river_graph = river_graph[0]
open_graph.close()
open_ts = open(topological_sort_location, "rb")
sorted_river_list = pickle.load(open_ts)
sorted_river_list = sorted_river_list[0]
open_ts.close()

k, filt_eff, s_eff, t_eff = [0.00995203, 0.99993389, 0.63221308, 0.81011981]

beta_0 = 1
gr, rl, cont, sc_nr, locs, obs = [river_graph, sorted_river_list, contamination_df, '', datapoint_locations,
                                  observed_values]
sim_results, discharges = simulated_contaminants(filt_eff, s_eff, t_eff, k, beta_0, gr, rl, cont, sc_nr,
                                                               locs, obs)
discharges = discharges/numpy.mean(discharges)
discharge_weights = numpy.sqrt(discharges)
implied_excretion = numpy.mean(obs)/numpy.mean(sim_results)
implied_excretion = numpy.mean(obs*discharge_weights) / numpy.mean(sim_results*discharge_weights)
#implied_excretion = numpy.sum(sim_results*obs*discharges) / numpy.sum(sim_results**2 * discharges)
sim_results *= implied_excretion
print(implied_excretion)
error, mean_error = error_formulae(obs, sim_results, discharges, option=0, weighted=1)
print(1 - error / mean_error)

plt.plot(obs, sim_results, 'o')
plt.plot(sim_results, sim_results)
plt.xlabel("Prediction")
plt.ylabel("Outcome")
plt.show()
plt.clf()

plt.plot(numpy.log(1+obs), numpy.log(1+sim_results), 'o')
plt.plot(numpy.log(1+sim_results), numpy.log(1+sim_results))
plt.xlabel("Prediction")
plt.ylabel("Outcome")
plt.show()
plt.clf()

plt.plot(obs*discharge_weights, sim_results*discharge_weights, 'o')
plt.plot(sim_results*discharge_weights, sim_results*discharge_weights)
plt.xlabel("Prediction")
plt.ylabel("Outcome")
plt.show()
plt.clf()

# log error
sim_results_log = numpy.log(1 + sim_results*discharges)
obs_log = numpy.log(1 + obs*discharges)
obs_log_mean = numpy.mean(obs_log)
log_error = numpy.square(sim_results_log - obs_log)
log_var = numpy.square(obs_log - obs_log_mean)
log_R_sqr = 1 - numpy.sum(log_error)/numpy.sum(log_var)
print("log squared error: " + str(log_R_sqr))

normal_error = numpy.square(sim_results-obs)
mean_error = numpy.square(obs-numpy.mean(obs))
normal_R_sqr = 1 - numpy.sum(normal_error)/numpy.sum(mean_error)
print("normal squared error: " + str(normal_R_sqr))

distance_error = numpy.abs(sim_results-obs)
absolute_dev = numpy.abs(obs-numpy.mean(obs))
distance_explained = 1 - numpy.sum(distance_error)/numpy.sum(absolute_dev)
print("distance to mean error: " + str(distance_explained))

absolute_dev = numpy.abs(obs-numpy.median(obs))
distance_explained = 1 - numpy.sum(distance_error)/numpy.sum(absolute_dev)
print("distance to median error: " + str(distance_explained))

percentage_error = numpy.abs(sim_results-obs)/(obs+1)
percentage_error = numpy.mean(percentage_error)
print("percentual error: " + str(percentage_error))

error_of_sum = 1 - numpy.sum(distance_error)/numpy.sum(obs)
print("error of sum: " + str(error_of_sum))

# create shapefile
discharges = discharges / numpy.mean(discharges)
dataframe = pandas.DataFrame()
dataframe['locations'] = locs
dataframe['Prediction'] = sim_results
dataframe['Observations'] = obs
dataframe['Error'] = (dataframe['Prediction'] - dataframe['Observations'])**2
dataframe['Error'] = dataframe['Error'] / numpy.mean(dataframe['Error'])
dataframe['weighted error'] = dataframe['Error'] * numpy.sqrt(discharges)
dataframe['weighted error'] = dataframe['weighted error'] / numpy.mean(dataframe['weighted error'])
shapefile_raster_functions.contaminant_to_shapefile(dataframe, reference_raster_location, output_name=calibrated_nonadapted_location)
