import networkx
import numpy
import pandas
import os
from src.library import shapefile_raster_functions
from return_contaminant import simulated_contaminants, error_formulae
from scipy.optimize import minimize
from src.library import graph_functions
from config.config import DATA_DIR
from osgeo import gdal



#set current directory as working directory
os.chdir(os.path.dirname(os.path.abspath(__file__)))


# Load data
directory = DATA_DIR
primary_eff = 0.33
secondary_eff = 0.7
tertiary_eff = 0.92
filter_eff = 1  # excludes in-situ load.



def objective_function(params, fixed):
    k, nothing = params

    # from European scale assessment of the potential of ozonation and activated carbon treatment to reduce
    # micropollutant emissions with wastewater. Pistocchi 2022
    gr, rl, cont, sc_nr, locs, obs, save_location, p_eff, s_eff, t_eff, f_eff = fixed
    sim_results, discharges = simulated_contaminants(f_eff, p_eff, s_eff, t_eff, k, gr, rl, cont, sc_nr,
                                                                   locs)
    excretion = numpy.sum(sim_results*obs)/numpy.sum(sim_results**2)
    sim_results *= excretion
    error, mean_error = error_formulae(obs, sim_results, discharges, option=0, weighted=0)
    print(error/mean_error)
    print(params)
    graph_functions.simple_save([error/mean_error, [k, excretion], sim_results, discharges], save_location)
    return error




# open the input files
graph_location = os.path.join(directory, "river_graph.pkl")
reference_raster_location = os.path.join(directory, "reference_raster.tif")
contamination_df_location = os.path.join(directory, "AGG_WWTP_df.csv")
observed_df_location = os.path.join(directory,"Occurences_in_river.csv")
scenario_number = ""

contamination_df = pandas.read_csv(contamination_df_location)
observed_df = pandas.read_csv(observed_df_location)

reference_raster = gdal.Open(reference_raster_location)
shapefile_raster_functions.give_pixel(observed_df['locations'][0], reference_raster, reverse=True)
# output name
output_name = 'results_calibration.csv'
## The code also generates shapefile and raster outputs for each contaminant. These files get the names from the columns
## in "Occurences_in_river.csv"
suppress_shapefile_raster_creation = False

datapoint_locations = observed_df["locations"]
datapoint_count = len(datapoint_locations)

RT = "RT_" + scenario_number
dis = scenario_number
if scenario_number == "":
    RT = "RT_HR"
    dis = "flow_HR"
try:
    river_graph, sorted_river_list, contamination_df = graph_functions.simple_load('graph_and_sort.pkl')
    print('the code managed to recover an older river graph. If this gives an error,'
          ' delete the file \'graph and sort\'.pkl')
except:
    river_graph = graph_functions.load_selected_attributes_graph(graph_location, [RT, dis, 'basin', 'x', 'y'])

    # cut river graph
    basins = set([])
    for location in datapoint_locations:
        basins.add(river_graph.nodes[location]['basin'])

    node_list = []
    for node in river_graph:
        if river_graph.nodes[node]['basin'] in basins:
            node_list.append(node)

    river_graph = networkx.subgraph(river_graph, node_list)
    sorted_river_list = list(networkx.topological_sort(river_graph))

    # cut the contamination_df to only leave those that (might) affect the observations.
    node_df = pandas.DataFrame(node_list)
    contamination_df = pandas.merge(contamination_df, node_df, left_on='pixel_number', right_on=0)
    node_df = node_list = None
    graph_functions.simple_save([river_graph.copy(), sorted_river_list, contamination_df], 'graph_and_sort.pkl')

bnds = ((0, 0.05), (0, 0))
starting_param = [0.01, 0]

contaminant_list = ['Lumped 14',	'Lumped 21',	'Lumped 61',
                    'Atenolol', 'Carbamazepine', 'Caffeine', 'Cetirizine',
                    'Citalopram', 'Codeine', 'Cotinine', 'Desvenlafaxine', 'Diltiazem', 'Fexofenadine', 'Gabapentin',
                    'Lidocaine', 'Metformin', 'Nicotine', 'Paracetamol', 'Propranolol', 'Ranitidine', 'Sitagliptin',
                    'Sulfamethoxazole', 'Trimethoprim', 'Venlafaxine']

result_dataframe = pandas.DataFrame()
result_dataframe.index = ['R^2', 'excretion', 'attenuation']
for contaminant in contaminant_list:
    save_location = contaminant + '.pkl'
    observed_values = observed_df[contaminant]
    fixed_args = [river_graph, sorted_river_list, contamination_df, scenario_number, datapoint_locations,
                  observed_values, save_location, primary_eff, secondary_eff, tertiary_eff, filter_eff]
    res = minimize(objective_function, numpy.array(starting_param),
           args=(fixed_args), method="Nelder-Mead", bounds=bnds)
    calibration_information = graph_functions.simple_load(save_location)
    # create shapefile
    if not suppress_shapefile_raster_creation:
        discharges = calibration_information[3]
        discharges_norm = discharges / numpy.mean(discharges)
        dataframe = pandas.DataFrame()
        dataframe['locations'] = observed_df['locations']
        dataframe['Prediction'] = calibration_information[2]
        dataframe['Observations'] = observed_values
        dataframe['discharge'] = discharges
        dataframe['Longitude'] = observed_df['longitude']
        dataframe['Latitude'] = observed_df['latitude']
        dataframe['Error'] = (dataframe['Prediction'] - dataframe['Observations']) ** 2
        dataframe['Error'] = dataframe['Error'] / numpy.mean(dataframe['Error'])
        dataframe['weighted error'] = dataframe['Error'] * numpy.sqrt(discharges_norm)
        dataframe['weighted error'] = dataframe['weighted error'] / numpy.mean(dataframe['weighted error'])
        shapefile_raster_functions.csv_to_shapefile(dataframe, reference_raster_location,
                                                            output_name=contaminant+".shp", options=False)
        
        if contaminant == 'Lumped 14':
            shapefile_raster_functions.csv_to_shapefile(dataframe, reference_raster_location, output_name=os.path.join(DATA_DIR, contaminant, contaminant + ".shp"), options=False)
        
        run_parameters = [calibration_information[1][1], calibration_information[1][0], filter_eff, primary_eff,
                          secondary_eff, tertiary_eff]
        river_graph = graph_functions.run_model(river_graph, sorted_river_list, contamination_df, run_parameters)
        graph_functions.print_graph(river_graph, ["concentration", "flow_HR"], reference_raster_location,
                    contaminant + ".tif")

    result_dataframe[contaminant] = [1 - calibration_information[0], calibration_information[1][1],
                                         calibration_information[1][0]]
    os.remove(save_location)


result_dataframe.to_csv(output_name)