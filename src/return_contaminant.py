import numpy
import math
import networkx


def simulated_contaminants(filtered_efficacy, secondary_efficacy, tertiary_efficacy, k_dec, b_0, river_graph,
                           sorted_river_list, contamination_df, scenario_number, datapoint_locations):

    # contamination
    treatment_efficacy = 0.3 * (contamination_df["Treatment_level"] == 1) + secondary_efficacy * \
                         (contamination_df["Treatment_level"] == 2) + tertiary_efficacy * \
                         (contamination_df["Treatment_level"] == 3)
    # formula
    contamination = (1 - treatment_efficacy) * contamination_df["Treat_a"] + (1 - filtered_efficacy)\
                    * contamination_df["Filt_a"] + contamination_df["Unfilt_a"]
    contamination *= contamination_df["pollution"]*b_0
    location = contamination_df["pixel_number"]

    cont = "Contaminant " + scenario_number
    RT = "RT_" + scenario_number
    dis = scenario_number
    rel_cont = "Relative contaminant " + scenario_number

    if scenario_number == "":
        cont = "Contaminant"
        RT = "RT_HR"
        dis = "flow_HR"
        rel_cont = "Relative Contaminant"
    networkx.set_node_attributes(river_graph, 0, name=cont)

    for i in range(len(contamination_df)):
        pixel_number = location[i]
        river_graph.nodes[pixel_number][cont] += contamination[i]

    for n in sorted_river_list:  # for all river pixels
        # Add the contamination of the parent cells
        parents = list(river_graph.predecessors(n))
        for k in parents:
            river_graph.nodes[n][cont] += river_graph.nodes[k][cont]

        river_graph.nodes[n][cont] *= math.exp((-k_dec * river_graph.nodes[n][RT]))
        river_graph.nodes[n][rel_cont] = river_graph.nodes[n][cont] / river_graph.nodes[n][dis]
    results = numpy.zeros([len(datapoint_locations)])
    discharges = numpy.zeros([len(datapoint_locations)])
    for i in range(len(datapoint_locations)):
        results[i] = river_graph.nodes[datapoint_locations[i]][rel_cont]
        discharges[i] = river_graph.nodes[datapoint_locations[i]][dis]

    return results, discharges


def error_formulae(observations, model_outcomes, discharges, option, weighted):
    if weighted == 1:
        observations = observations * numpy.sqrt(discharges)
        model_outcomes = model_outcomes * numpy.sqrt(discharges)
    if option == 0:  # squared error
        error_list = numpy.square(observations - model_outcomes)
        error = numpy.mean(error_list)
        mean_obs = numpy.mean(observations)
        mean_error_list = numpy.square(observations - mean_obs)
        mean_error = numpy.mean(mean_error_list)
        return [error, mean_error]

    if option == 1:  # absolute error
        error_list = numpy.sort(numpy.abs(observations - model_outcomes))
        error = numpy.mean(error_list)
        median_error_list = numpy.sort(numpy.abs(observations - numpy.median(observations)))
        median_error = numpy.mean(median_error_list)
        return [error, median_error]

    if option == 2:  # mean log error
        observations = numpy.log(observations + 1)
        model_outcomes = numpy.log(model_outcomes + 1)
        error_list = numpy.sort(numpy.square(observations - model_outcomes))
        error = numpy.mean(error_list)
        mean_error_list = numpy.sort(numpy.square(observations - numpy.mean(observations)))
        mean_error = numpy.mean(mean_error_list)
        return [error, mean_error]

