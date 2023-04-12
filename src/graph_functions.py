import os
import progressbar
import numpy
import math
import networkx
import pandas
#import src.library.shapefile_raster_functions as shapefile_raster_functions
from osgeo import gdal
from time import time
import pickle


def test():
    """
    test function for printing
    """
    print('aaaaaaa')
    pass

def add_RT_lakes(river_graph: networkx.DiGraph, sorted_graph: object, RT_name: str, discharge_name: str = 'flow_HR'
                 ) -> None:
    """
    add_RT_lakes replaces the residence time of the graph according to the lakes and their volumes.
    :rtype: None
    :river_graph: networkx.DiGraph: the river graph with fields 'lakes', 'volume', 'flow_HR' and some field RT_name.
    :sorted_graph: Iterator: a topological sort of the nodes in river_graph.
    :RT_name: string: The field to which the residence times should be written.
    :return: Saves the calculated residence times in the graph. Nothing is returned.
    """
    lakes_done = []
    sorted_graph_reverse = list(sorted_graph)
    sorted_graph_reverse.reverse()
    for i in sorted_graph_reverse:
        if river_graph.nodes[i]["lakes"] > 1:
            if river_graph.nodes[i]["lakes"] not in lakes_done:
                lakes_done.append(river_graph.nodes[i]["lakes"])
                river_graph.nodes[i][RT_name] = river_graph.nodes[i]["volume"] * 1000 ** 2 / river_graph.nodes[i][
                    discharge_name]
            else:
                river_graph.nodes[i][RT_name] = 0
    pass


def calculate_residence_time(discharge: float, slope: float, distance: float) -> float:
    """
    Calculate_residence_time gives the residence time of a cell according to discharge, slope and distance. The
    calculation uses coefficients from the following paper:
    https://adgeo.copernicus.org/articles/5/133/2005/adgeo-5-133-2005.pdf
    Schulze, Kerstin, Martin Hunger, and Petra Döll. "Simulating river flow velocity on global scale."
     Advances in Geosciences 5 (2005): 133-136.

    :rtype: float
    :discharge: float: the discharge of the cell in m^3/s
    :slope: float: The slope of the cell in m/m
    :distance: float: The length of the cell (depending on the direction of the river) in meters.
    :return: returns the residence time in hours.
    """
    # Estimation of velocities according to the cited paper
    inverse_manning_coefficient = 22.7
    a = 2.71
    b = 0.557
    c = 0.349
    d = 0.341

    w = a * math.pow(discharge, b)
    h = c * math.pow(discharge, d)

    r = (w * h) / (2 * h + w)
    v = inverse_manning_coefficient * math.pow(r, 2 / 3) * math.pow(slope, 1 / 2)

    # calculation of residence time in hours
    RT = distance / v
    RT_hours = 0.00027777777777 * RT

    return RT_hours


def simulate_waste_water(graph_location: object, contamination_df: pandas.DataFrame, sorted_river: object,
                         liter_per_equivalent: int = 150) -> networkx.DiGraph:
    """
    simulate_waste_water takes a contamination dataframe and propagates water coming from the treatment plants. It
    adjusts the discharge according to these water outflows.
    :rtype: networkx.DiGraph.
    :graph_location: string or networkx.DiGraph: the river_graph whose discharge will be adjusted.
    :contamination_df: string or pandas.DataFrame: the dataframe with the treatment plants.
    :sorted_river: iterator_object: A topological sort of the graph.
    :liter_per_equivalent: int: the amount of liters excreted per (adjusted) person equivalent.
    :return: returns the river graph with the adjusted discharges.
    """
    # preliminaries to check data types and load the data if file locations are supplied
    if isinstance(graph_location, str):
        open_graph = open(graph_location, "rb")
        river_graph = pickle.load(open_graph)

        open_graph.close()
    elif isinstance(graph_location, networkx.DiGraph):
        river_graph = graph_location
    else:
        raise TypeError('Pass either a Digraph or a string. your input was ' + str(type(graph_location)))

    if isinstance(contamination_df, str):
        contamination_df = pandas.read_csv(contamination_df)
    elif isinstance(contamination_df, pandas.DataFrame):
        pass
    else:
        raise TypeError("contamination_df needs to be either a string with the file location or a dataframe")

    liter_to_cube_meter = 1000
    hours_in_day = 24
    cubic_meter_hour = liter_per_equivalent / (liter_to_cube_meter * hours_in_day)

    # add water to the initial discharge points
    person_equivalents = contamination_df['Treat_a'] + contamination_df['Unfilt_a']
    networkx.set_node_attributes(river_graph, 0, name='additional water')

    for i in range(len(contamination_df)):
        river_graph.nodes[contamination_df['pixel_number'][i]]['additional water'] += \
            cubic_meter_hour * person_equivalents[i]

    # propogate the water
    for n in sorted_river:
        parents = list(river_graph.predecessors(n))
        for k in parents:
            river_graph.nodes[n]['additional water'] += river_graph.nodes[k]['additional water']

    # calculate the total discharge
    additional_water = list(networkx.get_node_attributes(river_graph, 'additional water').items())
    additional_water = numpy.array(additional_water)
    additional_water = additional_water[:, 1]
    discharge = list(networkx.get_node_attributes(river_graph, 'flow_HR').items())
    discharge = numpy.array(discharge)
    discharge = discharge[:, 1]
    total_discharge = discharge + additional_water

    # update the total discharge for the changed cells
    keys = numpy.array(list(networkx.get_node_attributes(river_graph, 'additional water').keys()))
    discharge = {}
    for key, value in zip(keys, total_discharge):
        discharge[key] = value
    networkx.set_node_attributes(river_graph, discharge, 'flow_HR')

    return river_graph


def extract_river_network(river_graph: networkx.DiGraph, pixel_number: int, include_progenitors: bool = 0) -> list:
    """
    extract_river_network takes a pixel_number, which must be a river, and extracts the entire upstream and downstream.
    :rtype: list
    :river_graph: networkx.DiGraph: The river_graph from which to extract the river network.
    :pixel_number: int: The pixel number of the river cel that is to be extracted.
    :include_progenitors: bool: if set to 1, this extracts the other upstream cells of downstream cells
     (like including the in-laws in a family tree)
    :return: a list of the pixel numbers in the river network.
    """
    if pixel_number not in river_graph:
        raise IndexError("The pixel_number is not in the river graph")

    # if include_progenitors is selected as an option, pick the most downstream cell as the pixel_number. This ensures
    # that the code that follows will extract the entire network, as the in-law cells are upstream cells for the most
    # downstream cell.
    if include_progenitors:
        child = list(river_graph.successors(pixel_number))
        while len(child) > 0:
            pixel_number = child[0]
            child = list(river_graph.successors(pixel_number))

    # extracting the river network
    node_network_ordered = [pixel_number]  # the list of the cells included in the network
    cells_to_consider = [pixel_number]  # a list that keeps track of the cells that are to be included

    # First, the parents are extracted
    while len(cells_to_consider) > 0:  # while there are parent cells
        element = cells_to_consider[0]  # take a cell that is to be included in the river network
        node_network_ordered.append(element)  # add the element to the network
        parents = list(river_graph.predecessors(element))
        for parent in parents:  # add parents to the cells to consider
            cells_to_consider.append(parent)

        del cells_to_consider[0]  # delete the cell that was considered
    node_network_ordered.reverse()  # reverse the network such that parent cells come first.

    # second, append all the children of the pixel number.
    child = list(river_graph.successors(pixel_number))
    while len(child) > 0:  # while the current cell has a child
        node_network_ordered.append(child[0])  # add current child to the network
        child = list(river_graph.successors(child[0]))  # update the current child

    return node_network_ordered

def describe_contamination_network(river_graph: networkx.DiGraph, ordered_node_network: list, reference_raster_location:
str = 'reference_raster.tif', scenario: str = '', output_name:
str = 'contaminant_network.tif') -> None:
    """
    describe_contamination_network takes a river network that must be ordered from parent to child and prints the
    network with detailed information on contamination.
    :rtype: None
    :river_graph: networkx.DiGraph: The river_graph from which to extract the river network.
    :ordered_node_network: list: List that contains the pixel numbers that form up a river network. The pixel numbers
    must be in topological order.
    :reference_raster_location: str: The reference raster that concords with the pixel numbers of river_graph.
    :scenario: str:
    :return: a list of the pixel numbers in the river network.
    """
    sub_river_graph = river_graph.subgraph(ordered_node_network).copy()
    networkx.set_node_attributes(sub_river_graph, 0, name='first_contaminator')
    networkx.set_node_attributes(sub_river_graph, 0, name='first_contaminator_value')
    networkx.set_node_attributes(sub_river_graph, 0, name='second_contaminator')
    networkx.set_node_attributes(sub_river_graph, 0, name='second_contaminator_value')
    networkx.set_node_attributes(sub_river_graph, 0, name='third_contaminator')
    networkx.set_node_attributes(sub_river_graph, 0, name='third_contaminator_value')

    if scenario == '':
        cont = "Contaminant"
        RT = "RT_HR"
        rel_cont = "Relative Contaminant"
    else:
        cont = "Contaminant " + scenario
        RT = "RT " + scenario
        rel_cont = "Relative Contaminant " + scenario

    # find some value for decay
    decay_rate = 0
    found = False
    for node in ordered_node_network:
        if sub_river_graph.nodes[node][cont] > 0:
            if sub_river_graph.nodes[node]['Initial contaminant'] == 0:
                predecessors = list(sub_river_graph.predecessors(node))
                if len(predecessors) == 1:
                    decay_factor = sub_river_graph.nodes[node][cont] / sub_river_graph.nodes[predecessors[0]][cont]
                    decay_rate = -numpy.log(decay_factor) / sub_river_graph.nodes[node][RT]
                    found = True
                    break
    if not found:
        print("no decay rate was found, the model must be empty")
        return None
    dtype = [('origin', int), ('contamination', float)]
    for node in ordered_node_network:  # for all river pixels
        # Add the contamination of the parent cells
        parents = list(sub_river_graph.predecessors(node))
        candidates = []
        initial_contaminant = sub_river_graph.nodes[node]['Initial contaminant']
        if initial_contaminant > 0:
            candidates.append((node, sub_river_graph.nodes[node]['Initial contaminant']))

        for k in parents:
            candidates.append((sub_river_graph.nodes[k]['first_contaminator'],
                               sub_river_graph.nodes[k]['first_contaminator_value']))
            candidates.append((sub_river_graph.nodes[k]['second_contaminator'],
                               sub_river_graph.nodes[k]['second_contaminator_value']))
            candidates.append((sub_river_graph.nodes[k]['third_contaminator'],
                               sub_river_graph.nodes[k]['third_contaminator_value']))

        candidates = numpy.array(candidates, dtype=dtype)
        candidates = numpy.sort(candidates, order='contamination')
        candidate_count = len(candidates)

        decay_factor = math.exp(-decay_rate * sub_river_graph.nodes[node][RT])
        try:
            sub_river_graph.nodes[node]['first_contaminator'] = candidates[candidate_count - 1]['origin']
            sub_river_graph.nodes[node]['first_contaminator_value'] = \
                candidates[candidate_count - 1]['contamination'] * decay_factor
        except IndexError:
            pass
        try:
            sub_river_graph.nodes[node]['second_contaminator'] = candidates[candidate_count - 2]['origin']
            sub_river_graph.nodes[node]['second_contaminator_value'] = \
                candidates[candidate_count - 2]['contamination'] * decay_factor
        except IndexError:
            pass
        try:
            sub_river_graph.nodes[node]['third_contaminator'] = candidates[candidate_count - 3]['origin']
            sub_river_graph.nodes[node]['third_contaminator_value'] = \
                candidates[candidate_count - 3]['contamination'] * decay_factor
        except IndexError:
            pass

    reference_raster = gdal.Open(reference_raster_location)
    contaminator_strings = ['first_contaminator', 'second_contaminator', 'third_contaminator']
    contaminator_value_strings = ['first_contaminator_value', 'second_contaminator_value',
                                  'third_contaminator_value']
    for node in ordered_node_network:  # for all river pixels
        for i in range(len(contaminator_strings)):
            if sub_river_graph.nodes[node][contaminator_strings[i]] != 0:
                [lat, long] = shapefile_raster_functions.give_pixel(
                    sub_river_graph.nodes[node][contaminator_strings[i]],
                    reference_raster, reverse=True)
                long = 100000 * int(long * 1000 + 0.5)
                long_lat_code = long + int(lat * 1000 + 0.5)
                sub_river_graph.nodes[node][contaminator_strings[i]] = long_lat_code

            try:
                sub_river_graph.nodes[node][contaminator_value_strings[i]] = 0.5 + 1000 * \
                                                                             sub_river_graph.nodes[node][
                                                                                 contaminator_value_strings[i]] / \
                                                                             (0.00001 + sub_river_graph.nodes[node][
                                                                                 cont])
            except (ZeroDivisionError, RuntimeWarning):
                sub_river_graph.nodes[node][contaminator_value_strings[i]] = 0
            sub_river_graph.nodes[node]['Relative contaminant x1000'] = sub_river_graph.nodes[node][rel_cont] * 1000

    attributes = ['Relative contaminant x1000', 'Initial contaminant', cont, 'first_contaminator',
                  'first_contaminator_value', 'second_contaminator', 'second_contaminator_value', 'third_contaminator',
                  'third_contaminator_value']
    print_sub_graph(sub_river_graph, attributes, reference_raster_location, output_name,
                    datatype=gdal.GDT_Int32)
    pass


def print_graph(graph_location: str, attribute_list: list, reference_raster_location: str, output_name: str,
                datatype=gdal.GDT_Float32):
    """
    This function takes in a graph, and converts it to a raster
    :rtype: gdal raster
    :graph_location: str: location of the raster
    :attribute_list: list: list of strings that specify the attributes to be copied
    :reference_raster_location: str: location of the reference raster
    :output_name: str: the location name that the output raster should have
    :datatype: gdal.Datatype: the datatype of the raster
    """
    if isinstance(graph_location, str):
        open_graph = open(graph_location, "rb")
        graph = pickle.load(open_graph)
        open_graph.close()
    elif isinstance(graph_location, networkx.DiGraph):
        graph = graph_location
    else:
        raise TypeError('Pass either a Digraph or a string. your input was ' + str(type(graph_location)))
    reference_raster = gdal.Open(reference_raster_location)
    count = len(attribute_list)
    indicator = False

    if count == 0:
        count = 1
        indicator = True

    no_data_val = -523521
    raster_matrix = numpy.zeros([reference_raster.RasterYSize, reference_raster.RasterXSize, count]) + no_data_val
    iterations_number = 0
    tot = graph.number_of_nodes()
    if indicator:
        for node_id in graph:
            [i, j] = [graph.nodes[node_id]["x"], graph.nodes[node_id]["y"]]
            raster_matrix[i, j] = 1
    else:
        for node_id in graph:
            iterations_number +=1
            for index in range(count):
                [i, j] = [graph.nodes[node_id]["x"], graph.nodes[node_id]["y"]]  # recover pixel coordinates from node
                raster_matrix[i, j, index] = graph.nodes[node_id][attribute_list[index]]

    gtiff_driver = gdal.GetDriverByName('GTiff')
    if not output_name.endswith('.tif'):
        output_name += '.tif'
    out_ds = gtiff_driver.Create(output_name, reference_raster.RasterXSize, reference_raster.RasterYSize,
                                 count, datatype)  # this creates a raster document with dimensions, bands, datatype

    out_ds.SetProjection(reference_raster.GetProjection())  # copy direction projection to output raster
    out_ds.SetGeoTransform(reference_raster.GetGeoTransform())  # copy direction resolution/location to output raster

    for index in range(1, count + 1):
        out_ds.GetRasterBand(index).WriteArray(raster_matrix[:, :, index - 1])
        try:
            out_ds.GetRasterBand(index).SetDescription(attribute_list[index - 1])
        except IndexError:
            out_ds.GetRasterBand(index).SetDescription('indicator')

        out_ds.GetRasterBand(index).SetNoDataValue(no_data_val)
    out_ds = None
    pass


def print_sub_graph(graph_location: str, attribute_list: list, reference_raster_location: str, output_name: str,
                    datatype=gdal.GDT_Float32):
    """
    This function takes in a graph, and converts it to a raster of minimum size.
    :rtype: gdal raster
    :graph_location: str: location of the raster
    :attribute_list: list: list of strings that specify the attributes to be copied
    :reference_raster_location: str: location of the reference raster
    :output_name: str: the location name that the output raster should have
    :datatype: gdal.Datatype: the datatype of the raster
    """
    if isinstance(graph_location, str):
        open_graph = open(graph_location, "rb")
        graph = pickle.load(open_graph)
        open_graph.close()
    elif isinstance(graph_location, networkx.DiGraph):
        graph = graph_location
    else:
        raise TypeError('Pass either a Digraph or a string. your input was ' + str(type(graph_location)))

    reference_raster = gdal.Open(reference_raster_location)
    count = len(attribute_list)
    indicator = False

    if count == 0:
        count = 1
        indicator = True

    no_data_val = -523521
    raster_matrix = numpy.zeros([reference_raster.RasterYSize, reference_raster.RasterXSize, count]) + no_data_val
    if indicator:
        for node_id in graph:
            [i, j] = [graph.nodes[node_id]["x"], graph.nodes[node_id]["y"]]
            raster_matrix[i, j] = 1
    else:
        for node_id in graph:
            for index in range(count):
                [i, j] = [graph.nodes[node_id]["x"], graph.nodes[node_id]["y"]]  # recover pixel coordinates from node
                raster_matrix[i, j, index] = graph.nodes[node_id][attribute_list[index]]

    # find minimum and maximum 0 index:

    non_zeros = raster_matrix[:, :, 0] - no_data_val
    non_zero_rows = numpy.flatnonzero(non_zeros)
    min_index = numpy.min(non_zero_rows)
    max_index = numpy.max(non_zero_rows)
    min_row = int(min_index / reference_raster.RasterXSize)
    max_row = int(max_index / reference_raster.RasterXSize) + 1

    non_zero_columns = numpy.flatnonzero(numpy.transpose(non_zeros))
    min_index = numpy.min(non_zero_columns)
    max_index = numpy.max(non_zero_columns)
    min_column = int(min_index / reference_raster.RasterYSize)
    max_column = int(max_index / reference_raster.RasterYSize) + 1

    raster_matrix = raster_matrix[min_row:max_row:1, min_column:max_column:1, :]  # cut the matrix
    row_count = max_row - min_row
    column_count = max_column - min_column
    column_count = int(column_count)
    gtiff_driver = gdal.GetDriverByName('GTiff')
    out_ds = gtiff_driver.Create(output_name, column_count, row_count,
                                 count, datatype)  # this creates a raster document with dimensions, bands, datatype

    out_ds.SetProjection(reference_raster.GetProjection())  # copy direction projection to output raster
    geo_transform = reference_raster.GetGeoTransform()
    upper_left_pixel = min_column + min_row * reference_raster.RasterXSize
    upper_left_pixel_coord = shapefile_raster_functions.give_pixel(upper_left_pixel, reference_raster, reverse=True)
    geo_transform = list(geo_transform)
    geo_transform[3] = upper_left_pixel_coord[0]
    geo_transform[0] = upper_left_pixel_coord[1]

    out_ds.SetGeoTransform(geo_transform)  # copy direction resolution/location to output raster

    for index in range(1, count + 1):
        out_ds.GetRasterBand(index).WriteArray(raster_matrix[:, :, index - 1])
        try:
            out_ds.GetRasterBand(index).SetDescription(attribute_list[index - 1])
        except IndexError:
            out_ds.GetRasterBand(index).SetDescription('indicator')

        out_ds.GetRasterBand(index).SetNoDataValue(no_data_val)
    out_ds = None


def extract_river_network_by_shapefile(graph_location: object, shapefile_location: str, reference_raster_location: str,
                                       include_children: bool = 0) -> list:
    """
    This function takes a shapefile and extracts all the river pixels within that shapefile. It also includes all the
    parents of the river cells. An additional option is provided to inlude the entire river network of all rivers that
     enter into the shapefile.
    This means that all parent cells of chilren that are not in the shapefile are included.
    :rtype: list
    :graph_location: string or networkx.DiGraph: location of the river graph or the river graph DiGraph object.
    :shapefile_location: string: the location of the shapefile used to extract the river network. The polygons in the
    shapefile must cover which rivers to extract.
    :reference_raster_location: string: A reference_raster that concords with the graph.
    :option: boolean: If this is set to 1, it also includes parents of children cells outside of the shapefile.
    :return: an ordered list of all nodes that are included in a river networks within the shapefile polygons.
    """

    # Load inputs
    if isinstance(graph_location, str):
        open_graph = open(graph_location, "rb")
        river_graph = pickle.load(open_graph)
        open_graph.close()
    elif isinstance(graph_location, networkx.DiGraph):
        river_graph = graph_location
    else:
        raise TypeError('Pass either a Digraph or a string. your input was ' + str(type(graph_location)))

    # extract an indicator raster from the polygons
    area_raster_matrix = shapefile_raster_functions.shapefile_to_raster(shapefile_location, reference_raster_location,
                                                                        "temporary.tif", option=0)
    # first, save all nodes that lie in the network
    rows, columns = numpy.shape(area_raster_matrix)
    area_raster = area_raster_matrix.flatten()  # with a correct reference raster, the indices coincide with the node id
    inclusion_list = []  # nodes which to include in the output
    check_list = []  # nodes which can have parents not within the bounds of the shapefile.
    for node in river_graph:
        if area_raster[node] == 1:  # if node is in shapefile
            inclusion_list.append(node)  # keep the node

            # the next code creates a list of nodes whose parents may have to be included in the graph.
            i = int(node / columns)
            j = node - i * columns
            if i == 0 or j == 0 or i == rows or j == columns:  # ensures the next if statement does not raise an error
                check_list.append(node)
            elif numpy.sum(area_raster_matrix[i - 1:i + 2:1,
                           j - 1:j + 2:1]) != 9:  # if not all adajacent nodes are in the graph
                check_list.append(node)  # add the node to check since its parents may not be in the shapefile.

    # next include parents of nodes
    for node in check_list:
        parents_to_consider = list(river_graph.predecessors(node))
        while len(parents_to_consider) > 0:  # while there are still parents to consider
            parent = parents_to_consider[0]  # take the first parent in the list
            del parents_to_consider[0]  # and delete that first parent
            if parent not in inclusion_list:
                inclusion_list.append(parent)  # include the parent
                parents_to_consider += list(river_graph.predecessors(parent))  # add the parents of the parents
    if include_children:  # includes children of all the river cells if set to true.
        for node in check_list:
            child = list(river_graph.successors(node))
            not_in_list = True
            while not_in_list:
                not_in_list = False
                if child[0] not in inclusion_list:
                    not_in_list = True
                    inclusion_list.append(child[0])
                    child = list(river_graph.successors(child[0]))

    # order the list
    sub_river_graph = networkx.subgraph(river_graph, inclusion_list)  # first define graph with the list
    ordered_inclusion_list = list(networkx.topological_sort(sub_river_graph))  # then order the graph
    return ordered_inclusion_list


def extract_plant_network(graph_location: object, contamination_location: object, longitude: float,
                          model_parameters: list) -> list:
    """
     This function takes a treatment plant, given by the contamination dataframe and the longitude, and extracts
     the river network that the treatment plant is in. It supplies the river network with various indicators of
     contamination.
     :rtype: list
     :graph_location: string or networkx.DiGraph: location of the river graph or the river graph DiGraph object.
     :contamination_location: string: the location of the contamination dataframe on the computer.
     :longitude: float: the longitude of the treatment plant; serves as the unique identifier for the treatment plant.
     :model_parameters: list: Contains all the parameters that determine the model, such as excretion, attenuation etc.
     :return: A list that contains a sub-graph of the supplied graph as the first element. The second element is a list
     of attributes that are included in the sub-graph.
     """
    # Load inputs
    if isinstance(graph_location, str):
        open_graph = open(graph_location, "rb")
        river_graph = pickle.load(open_graph)
        river_graph = graph[0]
        open_graph.close()
    elif isinstance(graph_location, networkx.DiGraph):
        river_graph = graph_location
    else:
        raise TypeError('Pass either a Digraph or a string. your input was ' + str(type(graph_location)))

    if isinstance(contamination_location, str):
        contamination_df = pandas.read_csv(contamination_location)
    elif isinstance(contamination_location, pandas.DataFrame):
        contamination_df = contamination_location
    else:
        print("the contamination location object was neither a string nor a dataframe")
        raise ValueError

    # find treatment plants by longitude
    longitude_str = str(longitude)
    whole, decimals = longitude_str.split('.')
    multiplier = 10 ** len(decimals)
    candidates = (multiplier * contamination_df['dcpLongitu'] + 0.5)
    indices = (candidates.to_numpy(int) - 1 <= int(multiplier * longitude)) & \
              (int(multiplier * longitude) <= candidates.to_numpy(int))  # allow for imprecise longitudes
    treatment_plant = contamination_df[indices]

    if len(treatment_plant) > 1:
        print("more than one treatment plant has this longitude up to the specified decimal. Provide more decimals")
        raise Exception
    elif len(treatment_plant) == 0:
        print("No treatment plant had the longitude specified")

    # extract the downstream of the plant
    starting_pixel = treatment_plant['pixel_number'].iloc[0]
    river_list = []
    child = list(river_graph.successors(starting_pixel))
    while len(child) > 0:
        river_list.append(child[0])
        child = list(river_graph.successors(child[0]))
    full_river = [starting_pixel] + river_list
    sub_graph = networkx.subgraph(river_graph, full_river)  # make a new graph containing the downstream of the plant.

    excretion, k, primary, secondary, tertiary, filt = model_parameters

    # derive plant contamination and the removal percentage of an upgrade (e.g. from 0.6 to 0.8 means an additional 50
    # percent of treatment.
    if treatment_plant['Treatment_level'].iloc[0] == 1:
        sub_graph.nodes[starting_pixel]["Plant contamination"] = treatment_plant['Treat_a'].iloc[0] * (1 - primary)
        upgrade_multiplier = (1 - secondary) / (1 - primary)
    elif treatment_plant['Treatment_level'].iloc[0] == 2:
        sub_graph.nodes[starting_pixel]["Plant contamination"] = treatment_plant['Treat_a'].iloc[0] * (1 - secondary)
        upgrade_multiplier = (1 - tertiary) / (1 - secondary)
    elif treatment_plant['Treatment_level'].iloc[0] == 3:
        sub_graph.nodes[starting_pixel]["Plant contamination"] = treatment_plant['Treat_a'].iloc[0] * (1 - tertiary)
        upgrade_multiplier = None
    else:  # do not accept treatment level 0, since this implies that there is no treatment plant.
        print("Treatment level was not 1, 2 or 3.")
        raise ValueError

    # calculate additional fields for the graph that contain information on the treatment plant
    # plant concentration calculated the contamination concentration that is due to the plant.
    # % contamination remaining calculates the percentage of the contamination of the plant still remaining.
    # % from plant calculates the percentage of contamination that is due to the plant.
    # if an upgrade is possible, upgrade plant concentration is the amount of contaminant concentration due to the plant
    # if it were upgraded
    #
    # if an upgrade is possible, upgrade concentration gives the total contaminant concentration if the plant were
    # upgraded
    sub_graph.nodes[starting_pixel]["Plant concentration"] = \
        sub_graph.nodes[starting_pixel]["Plant contamination"] / sub_graph.nodes[starting_pixel]["flow_HR"]
    sub_graph.nodes[starting_pixel]["% Contamination remaining"] = 1
    sub_graph.nodes[starting_pixel]["% from plant"] = \
        sub_graph.nodes[starting_pixel]["Plant concentration"] / sub_graph.nodes[starting_pixel]["Relative Contaminant"]
    if upgrade_multiplier is not None:
        sub_graph.nodes[starting_pixel]["Upgrade plant concentration"] = \
            upgrade_multiplier * sub_graph.nodes[starting_pixel]["Plant concentration"]

        sub_graph.nodes[starting_pixel]["Upgrade concentration"] = \
            sub_graph.nodes[starting_pixel]["Relative Contaminant"] + \
            sub_graph.nodes[starting_pixel]["Upgrade plant concentration"] - \
            sub_graph.nodes[starting_pixel]["Plant concentration"]

    for node in river_list:
        predecessor = list(sub_graph.predecessors(node))
        sub_graph.nodes[node]["Plant contamination"] = sub_graph.nodes[predecessor[0]]["Plant contamination"]
        sub_graph.nodes[node]["Plant contamination"] *= math.exp(-k * sub_graph.nodes[node]["RT_HR"])
        sub_graph.nodes[node]["Plant concentration"] =\
            sub_graph.nodes[node]["Plant contamination"] / sub_graph.nodes[node]["flow_HR"]

        sub_graph.nodes[node]["% Contamination remaining"] =\
            sub_graph.nodes[node]["Plant contamination"] / sub_graph.nodes[starting_pixel]["Plant contamination"]

        sub_graph.nodes[node]["% from plant"] =\
            sub_graph.nodes[node]["Plant concentration"] / sub_graph.nodes[node]["Relative Contaminant"]

        if upgrade_multiplier is not None:
            sub_graph.nodes[node]["Upgrade plant concentration"] = upgrade_multiplier * \
                                                                   sub_graph.nodes[node]["Plant concentration"]
            sub_graph.nodes[node]["Upgrade concentration"] =\
                sub_graph.nodes[node]["Relative Contaminant"] + sub_graph.nodes[node]["Upgrade plant concentration"] - \
                sub_graph.nodes[node]["Plant concentration"]

    # define additional output and return
    fields = ["Plant contamination", "Plant concentration", "% Contamination remaining", "% from plant", "Contaminant",
              "Relative Contaminant", "Upgrade plant concentration", "Upgrade concentration"]

    return [sub_graph, fields]


def find_largest_contamination_sites(graph_location: object, count: int, minimal_length: int,
                                     dilution_floor: float) -> list:
    """
     This function takes a river graph and derives the river networks of those rivers with the largest contamination.
     the size of those river networks depends on the minimum length. Those river networks that that are diluted by 1/
     dilution_floor before the minimum length in pixels, are not counted. (e.g. if dilution_floor is 0.2, the initial
     contamination may not be diluted by more than 500% before minimal_length pixels. Otherwise, another contamination
     site is chosen).
     :rtype: list
     :graph_location: string or networkx.DiGraph: location of the river graph or the river graph DiGraph object.
     :count: int: the amount of river networks to return
     :minimal_length: int: The minimal length of a contamination site (as defined by the dilution_floor) to be included.
     :dilution_floor: float: Determines the length by requiring that the initial contaminant may not dilute by more than
     1/dilution_floor.
     :return: A list of highly contaminated river networks.
     """
    # Load inputs
    if isinstance(graph_location, str):
        open_graph = open(graph_location, "rb")
        river_graph = pickle.load(open_graph)
        open_graph.close()
    elif isinstance(graph_location, networkx.DiGraph):
        river_graph = graph_location
    else:
        raise TypeError('Pass either a Digraph or a string. your input was ' + str(type(graph_location)))

    # obtain a pandas dataframe from the graph and sort it.
    sub_graph_df = pandas.DataFrame.from_dict(dict(river_graph.nodes(data=True)), orient='index')
    sub_graph_df = sub_graph_df.sort_values(by=['Relative Contaminant'], ascending=False)
    row_names = list(sub_graph_df.index)
    sub_graph_df['pixel_number'] = row_names

    # extracting the rivers
    river_networks_extracted = 0
    candidate_index = 0
    nodes_accepted = []
    considered_nodes = []  # prevents nodes to be considered that were already checked
    while river_networks_extracted < count:  # while there remain rivers to be extracted
        pixel_number = sub_graph_df['pixel_number'].iloc[candidate_index]
        if pixel_number in considered_nodes:  # ignore those numbers that were considered
            candidate_index += 1
            continue
        contamination_concentration = sub_graph_df['Relative Contaminant'].iloc[candidate_index]
        river_network_list = []
        child = [pixel_number]
        remainder = 1
        while remainder > dilution_floor and len(child) > 0:  # continue network if there is still contamination:
            river_network_list.append(child[0])
            remainder = river_graph.nodes[child[0]]['Relative Contaminant'] / contamination_concentration
            child = list(river_graph.successors(child[0]))

        if len(river_network_list) > minimal_length:  # this means the network is accepted. Include the parent cells.
            inclusion_list = []
            parents_to_consider = list(river_graph.predecessors(pixel_number))
            while len(parents_to_consider) > 0:
                parent = parents_to_consider[0]
                del parents_to_consider[0]
                if parent not in nodes_accepted:
                    inclusion_list.append(parent)
                    parents_to_consider += list(river_graph.predecessors(parent))
            inclusion_list.reverse()  # ensures topological order
            river_network_list = inclusion_list + river_network_list  # add to the rivers we wish to extract
            nodes_accepted += river_network_list
            river_networks_extracted += 1  # update the current count of rivers extracted

        candidate_index += 1  # update the index that is tried, even if no network is extracted
        considered_nodes += river_network_list  # ensures nodes are not reconsidered if they already were once.

    return nodes_accepted


def add_scenario(river_graph, shp_rivers, scen_names: list = [], all_scenarios: int = 0) -> None:
    """
     This function takes a river graph, given by the HydroRIVER id, and adds the equivalent flow values of the
     scenarios required which are extracted from the shapefile. It also computes the RT equivalent value.
     :rtype: None
     :river_graph: networkx.DiGraph: river graph DiGraph object.
     :shp_river: geopandas GeoDataFrame: Shapefile of the river network with the scenarios values.
     :scen_names: list: list of the scenario names.
     :all_scenarios: int: 1 if all the scenarios needs to be added, 0 otherwise.
     :return: None: The graph is updated with scenarios.
     """
    # progress bar while adding scenarios
    bar = progressbar.ProgressBar(maxval=len(river_graph.nodes), widgets=[progressbar.Bar('=', '[', ']'), ' ',
                                                                          progressbar.Percentage()])
    bar.start()

    # Transform graph to csv to merge with shapefile.
    river_csv = graph_to_csv(river_graph, "temp.csv")
    # HYRIV_ID from graph are shapefile HYRIV_ID - 1.
    river_csv["HYRIV_ID"] += 1
    # Merge shapefile and graph data to optimize the code.
    merged_df = pandas.merge(river_csv, shp_rivers, left_on="HYRIV_ID", right_on="HYRIV_ID")

    # Add scenarios
    if all_scenarios == 0:
        num_scenarios = len(scen_names)
    else:
        # in our case the variables starting with f contains the scenario flow values.
        names = shp_rivers.columns.tolist()
        scen_names = list(filter(lambda x: x.startswith("f"), names))
        num_scenarios = len(scen_names)

    for i in range(len(merged_df)):
        bar.update(i + 1)
        pixel_number = merged_df["pixel_number"][i]
        for k in range(num_scenarios):
            flow = scen_names[k]  # scenario name.
            flow_value = merged_df[flow][i]  # river scenario flow.
            # Avoid None flow problems.
            try:
                flow_value = float(flow_value)
            except (ValueError, TypeError):
                flow_value = river_graph.nodes[pixel_number]['flow_HR']

            # Minimum flow value 0.01 to allow calculation of RT.
            if flow_value < 0.01:
                flow_value = 0.01
            # Assign flow value.
            river_graph.nodes[pixel_number][flow] = flow_value
            RT = "RT_" + flow
            # Compute and assign RT value.
            slope = river_graph.nodes[pixel_number]['slope']
            distance = river_graph.nodes[pixel_number]['cell_distance']
            river_graph.nodes[pixel_number][RT] = calculate_residence_time(flow_value, slope, distance)
    bar.finish()
    pass


def absorb_raster(river_graph: networkx.DiGraph, raster_location: str, attribute_names: list = [],
                  datatype=numpy.double) -> networkx.DiGraph:
    """
     This code takes a raster and a river graph, and puts the attributes of the raster specified by attribute_names into
     the river_graph.
     :rtype: networkx.DiGraph
     :river_graph: networkx.DiGraph: river graph DiGraph object.
     :raster_location: str: Location of the raster on the computer
     :attribute_names: list: list of the names of the attribute that are to be included. These should be the
     descriptions of the band as seen for instance in QGIS.
     :return: river graph with new attributes according to the attribute_names and the raster.
     """

    list_of_raster_matrices = []  # matrices with values to absorb
    ordered_names = []  # orders the attribute_names in order of appearance.
    raster = gdal.Open(raster_location)

    # extract raster matrices according to their name and store the order.
    if len(attribute_names) > 0:
        for band in range(raster.RasterCount):
            if raster.GetRasterBand(band + 1).GetDescription() in attribute_names:
                ordered_names.append(raster.GetRasterBand(band + 1).GetDescription())
                # include matrix of data. The flattened entry must concord with the pixel number.
                list_of_raster_matrices.append(raster.GetRasterBand(band + 1).ReadAsArray().flatten())

        if len(list_of_raster_matrices) != len(attribute_names):
            raise ValueError('Some attribute names were not found')

    # write raster values to the graph
    for node in river_graph:
        for i in range(len(ordered_names)):
            array = numpy.array([1], dtype=datatype)
            array[0] = numpy.asarray(list_of_raster_matrices[i][node], dtype=datatype)
            river_graph.nodes[node][ordered_names[i]] = array[0]

    return river_graph


def absorb_shapefile(river_graph: networkx.DiGraph, shapefile_location: str, attribute_names: list = [],
                     reference_raster_location: str = "reference_raster.tif", datatype=numpy.double)\
        -> networkx.DiGraph:
    """
     This code takes a shapefile and a river graph, and puts the attributes of the shapefile specified by
     attribute_names into the river_graph.
     :rtype: networkx.DiGraph
     :river_graph: networkx.DiGraph: river graph DiGraph object.
     :shapefile_location: str: Location of the shapefile on the computer
     :attribute_names: list: list of the names of the attribute that are to be included. These should be the
     descriptions of the shapefile layers as seen for instance in QGIS.
     :return: river graph with new attributes according to the attribute_names and the shapefile.
     """

    temp_output_name = str(time()) + '.tif'  # ensures an untaken output name

    # create a raster from the shapefile
    shapefile_raster_functions.shapefile_to_raster(shapefile_location, reference_raster_location, temp_output_name,
                                                   attribute_names, include_ind=0)
    # Write the raster to the shapefile
    river_graph = absorb_raster(river_graph, temp_output_name, attribute_names, datatype)
    os.remove(temp_output_name)  # remove the temporary raster
    return river_graph


def graph_to_csv(river_graph: object, output_name: str = "", save: bool = True) -> pandas.DataFrame:
    """
     This code takes a river graph and converts it to a dataframe (without topological information) with output_name.
     :rtype: pandas.DataFrame
     :river_graph: networkx.DiGraph or str: the river graph or the location on the computer
     :output_name: str: name to which to write the csv of the river_graph
     :return: dataframe with information on the nodes of the river graph. Does not contain topological information.
     """
    # interpreting the input
    if isinstance(river_graph, str):
        river_graph_string = river_graph.split('.')
        output_name = river_graph_string[0]
        open_graph = open(river_graph, "rb")
        graph = pickle.load(open_graph)
        try:
            graph = graph[0]
        except KeyError:
            pass
        open_graph.close()
    elif isinstance(river_graph, networkx.DiGraph):
        graph = river_graph
    else:
        raise TypeError('Pass either a Digraph or a string. your input was ' + str(type(river_graph)))

    # converting to a dataframe
    dataframe = pandas.DataFrame.from_dict(dict(graph.nodes(data=True)), orient='index')
    row_names = list(dataframe.index)  # node pixel numbers appear as row names
    dataframe['pixel_number'] = row_names  # give pixel_numbers their own field

    # rename the rows to 1, 2, 3, etc.
    rename_dictionary = {row_names[i]: i for i in range(len(row_names))}
    dataframe.rename(mapper=rename_dictionary, inplace=True)

    if output_name == "" and save:  # provides functionally random output name.
        output_name = time()

    if not output_name.endswith('.csv'):  # adds csv to the name if forgotten.
        output_name += '.csv'
    dataframe.to_csv(output_name)  # save the dataframe
    return dataframe


def contaminators_benefit_cost(basin_graph: networkx.DiGraph, nodes: list, contamination_df: pandas.DataFrame,
                               parameters: list, cost_df: pandas.DataFrame, scenario_number: str =''):
    """
     this code takes a basin with a collection of nodes and calculates the relationship between those nodes and the
     contaminators. It calculates the benefits of upgrading the discharge points with respect to the nodes and it
     calculates the costs.
     :rtype: pandas.DataFrame
     :basin_graph: networkx.DiGraph: the river graph of only the basin.
     :nodes: list: a list of nodes (integer ids) that need to be solved. The overarching idea is that decreasing
     contamination in these nodes solves contamination for the whole basin network.
     :contamination_df: pandas.DataFrame: a dataframe that contains information on all the discharge points for the
     current basin. The dataframe should not contain discharge points that do not belong to the basin.
     :parameters: list: the list of the parameters. they are [excretion, attenuation, filt_efficacy, primary_efficacy,
     secondary_efficacy, tertiary_efficacy]
     :cost_df: A dataframe that contains information on the cost of the upgrades.
     :scenario_number: A string that contains the scenario that is under consideration.
     :return: for each node a list of benefits of upgrading each treatment plant. a list of costs for upgrading each
     treatment plant. The (water) discharges for each point. A changed contamination_df that contains upgrades instead
     of treatment plants.
     """
    excretion, attenuation, filt_eff, primary_eff, secondary_eff, tertiary_eff = parameters
    rem = 'remainder'  # records the distance measured in attenuation between pixels
    lvl = 'attenuation level'  # integer id that resets if two parts are too distant in terms of residence time
    path_info = 'path'
    networkx.set_node_attributes(basin_graph, 0, name=rem)
    networkx.set_node_attributes(basin_graph, 0, name=lvl)
    networkx.set_node_attributes(basin_graph, [], name=path_info)

    sorted_river_list = list(networkx.topological_sort(basin_graph))

    # Write path values to each river node, such that it can be identified whether nodes are connected (in a directed
    # sense)
    unique_number = 0
    for n in sorted_river_list:  # consider nodes in topological order
        parents = list(basin_graph.predecessors(n))
        basin_graph.nodes[n][path_info] = []  # needed to ensure list is empty at start, it randomly fills and I do not
        # know why.

        if parents: # if it has parents, put parent ids into the list
            for parent in parents:
                basin_graph.nodes[n][path_info] += basin_graph.nodes[parent][path_info]
        else:  # otherwise give the list a new unique value
            basin_graph.nodes[n][path_info] = [unique_number]
            unique_number += 1

    # Next, run the river graph in reverse such that the attenuation between any two points may be derived.

    sorted_river_list.reverse()  # reverse list to run in reverse topological order (downstream to upstream)
    basin_graph.nodes[sorted_river_list[0]][rem] = 1  # most downstream cell has value 1
    RT = "RT_HR"
    dis = "RT_Flow"
    if scenario_number != '':
        RT = "RT_" + scenario_number
        dis = scenario_number

    for n in sorted_river_list:  # run the river graph in reverse
        parents = list(basin_graph.predecessors(n))
        for parent in parents:
            # calculate what remains of the contamination in the cell
            basin_graph.nodes[parent][rem] = basin_graph.nodes[n][rem] * \
                                        math.exp(-attenuation * basin_graph.nodes[n][RT])
            basin_graph.nodes[parent][lvl] = basin_graph.nodes[n][lvl]  # give parent and son the same zone
            if basin_graph.nodes[parent][rem] < 0.00001:  # if the remainder is negligible
                basin_graph.nodes[parent][rem] = 1  # reset the count to 1
                basin_graph.nodes[parent][lvl] += 1  # and increment the zone with respect to the zone of the parent

    # Here the remainder from a discharge point to a pixel point is calculated
    remainder = [] # dimension is number of nodes
    for node in nodes:  # for all nodes
        remainder_for_node = []  # dimension is amount of discharge points
        for pixel_number in contamination_df['pixel_number']:  # for all discharges
            # equal level determines whether the discharge and the node are in the same zone
            equal_level = basin_graph.nodes[pixel_number][lvl] == basin_graph.nodes[node][lvl]

            path_value = basin_graph.nodes[pixel_number][path_info][0]

            # if the nodes are upstream of one another and in the same zone, calculate the remainder
            if path_value in basin_graph.nodes[node][path_info] and equal_level:
                remainder_value = basin_graph.nodes[pixel_number][rem] / basin_graph.nodes[node][rem]
                remainder_value *= math.exp(-attenuation * basin_graph.nodes[pixel_number][RT])
                if remainder_value <= 1:
                    remainder_for_node.append(remainder_value)
                else:  # if remainder is larger than 1, then discharge point is downstream, so no connection.
                    remainder_for_node.append(0)
            else:  # if not upstream of one another or if not in the same zone.
                remainder_for_node.append(0)
        remainder.append(remainder_for_node)

    # triple the dataframe for all possible upgrades (to primary, to secondary, to tertiary)
    dataframe_list = [contamination_df.copy(), contamination_df.copy(), contamination_df]
    remainder = [remainder_list + remainder_list + remainder_list for remainder_list in remainder]
    for i in range(len(dataframe_list)):
        dataframe_list[i]['upgrade'] = i + 1  # gives each discharge point a value of 1, 2 and 3.
    contamination_df = pandas.concat(dataframe_list, axis=0)
    contamination_df.reset_index(inplace=True)

    # defines efficacies for upgrading once, twice or three times. if upgrading once, twice or thrice is not possible,
    # then the upgrade efficiency is given as 0.
    upgrade_one = primary_eff * (contamination_df["Treatment_level"] == 0) + \
                  (secondary_eff - primary_eff) * (contamination_df["Treatment_level"] == 1) + \
                  (tertiary_eff - secondary_eff) * (contamination_df["Treatment_level"] == 2)
    upgrade_two = secondary_eff * (contamination_df["Treatment_level"] == 0) + \
                  (tertiary_eff - primary_eff) * (contamination_df["Treatment_level"] == 1)
    upgrade_three = tertiary_eff * (contamination_df["Treatment_level"] == 0)
    upgrade_eff = upgrade_one * (contamination_df['upgrade'] == 1) + upgrade_two * (contamination_df['upgrade'] == 2)\
                  + upgrade_three * (contamination_df['upgrade'] == 3)

    # give benefit of an upgrade as increase in treatment efficacy and persons treated.
    benefit = [val * upgrade_eff * (contamination_df['Treat_a'] + contamination_df['Unfilt_a']) for val in remainder]
    discharges = [basin_graph.nodes[node][dis] for node in nodes]  # put discharges in a list to return
    cost = 0 * contamination_df['Treat_a']  # gives empty series with correct length and numbering.

    # fill in the costs
    for i in range(len(contamination_df)):
        current_treatment = contamination_df['Treatment_level'][i]
        end_treatment = current_treatment + contamination_df['upgrade'][i]
        if end_treatment > 3:  # set cost to 0 if amount of upgrades not possible (later deleted)
            cost[i] = 0
            continue
        else:
            person_equivalents = contamination_df['Treat_a'][i] + contamination_df['Unfilt_a'][i]
            country_ID = contamination_df['countryID'][i]
            # the function below calculates the cost of upgrading
            cost[i] = calc_cost(current_treatment, end_treatment, person_equivalents, country_ID, cost_df)

    # thin out the options by deleting some unfeasible ones or some without benefit
    drop_list = []
    for i in range(len(cost)):
        if cost[i] == 0:  # means the upgrade is impossible
            drop_list.append(i)
        else:  # if upgrade is possible, let us see if there is some benefit
            some_benefit = False
            for benefit_list in benefit:
                if benefit_list[i] > 0:  # if it benefits at least one point
                    some_benefit = True
            if not some_benefit:  # if no benefit, drop.
                drop_list.append(i)
    indexes_to_keep = list(set(range(len(cost))) - set(drop_list))  # determines which upgrades are unviable or useless

    index = 0
    for benefit_list in benefit: # benefit list is the benefit to a single point.
        benefit[index] = benefit_list.take(indexes_to_keep)  # keep only viable or useful upgrades
        index += 1
    cost = cost.take(indexes_to_keep)
    contamination_df = contamination_df.take(indexes_to_keep)
    return benefit, cost, discharges, contamination_df

def calc_runoff_max_after_parents(basin_graph: networkx.DiGraph, contamination_df: pandas.DataFrame, parameters: list,
                                  limit_cont: float, global_min: float, scenario: str ='', children_number: int = 0)\
                                  -> list:
    """
     This code takes a river graph that contains a unique basin and a contamination dataframe that only has discharge
     points related to that basin. It then calculates contamination within the basin and identifies a few problematic
     nodes. These nodes are a minimal sufficient set of nodes in the sense that if these nodes are brought under the
     contamination limit, all nodes will satisfy the contamination limit.
     :rtype: list
     :basin_graph: networkx.DiGraph: the river graph of only the basin.
     :contamination_df: pandas.DataFrame: a list of nodes (integer ids) that need to be solved. The overarching idea is that decreasing
     contamination in these nodes solves contamination for the whole basin network.
     :contamination_df: pandas.DataFrame: a dataframe that contains information on all the discharge points for the
     current basin. The dataframe should not contain discharge points that do not belong to the basin.
     :parameters: list: the list of the parameters. they are [excretion, attenuation, filt_efficacy, primary_efficacy,
     secondary_efficacy, tertiary_efficacy]
     :cost_df: A dataframe that contains information on the cost of the upgrades.
     :scenario_number: A string that contains the scenario that is under consideration.
     :return: for each node a list of benefits of upgrading each treatment plant. a list of costs for upgrading each
     treatment plant. The (water) discharges for each point. A changed contamination_df that contains upgrades instead
     of treatment plants.
     """
    excretion, attenuation, filt_eff, primary_eff, secondary_eff, tertiary_eff = parameters
    pixel_list = list(networkx.topological_sort(basin_graph))
    if scenario == '':
        cont = "Contaminant"
        RT = "RT_HR"
        dis = "flow_HR"
        rel_cont = "Relative Contaminant"
        initial_cont = "Initial contaminant"
    else:
        cont = "Contaminant " + scenario
        RT = "RT_" + scenario
        #dis = "discharge " + scenario
        dis = scenario
        rel_cont = "Relative contaminant " + scenario
        initial_cont = "Initial contaminant " + scenario

    # contamination
    treatment_efficacy = primary_eff * (contamination_df["Treatment_level"] == 1) + secondary_eff * \
                     (contamination_df["Treatment_level"] == 2) + tertiary_eff * \
                     (contamination_df["Treatment_level"] == 3)
    contamination = (1 - treatment_efficacy) * contamination_df["Treat_a"] + (1 - filt_eff) \
                * contamination_df["Filt_a"] + contamination_df["Unfilt_a"]
    contamination *= contamination_df["pollution"] * excretion
    location = contamination_df["pixel_number"]

    # define names and fields
    sewage_volume = 'sewage_volume'  # will contain water generated from discharge point
    add_volume = 'add_volume'  # will contain water added to the rivers from the treatment plant
    min_water_volume = 'min_water_volume'  # the minimum amount of water required to satisfy the contamination limit
    potential_problem = "potential_problem"
    networkx.set_node_attributes(river_graph, 0, name=cont)
    networkx.set_node_attributes(river_graph, 0, name=sewage_volume)
    networkx.set_node_attributes(river_graph, 0, name=add_volume)
    networkx.set_node_attributes(river_graph, 0, name=min_water_volume)
    networkx.set_node_attributes(river_graph, False, name=potential_problem)

    cbm_per_equivalent = 0.0065  # cubic meters of water per person equivalent, a little over 150l

    for i in range(len(contamination_df)):
        pixel_number = location[i]
        river_graph.nodes[pixel_number][cont] += contamination[i]
        river_graph.nodes[pixel_number][sewage_volume] += (contamination_df['Treat_a'][i] +\
                                                          contamination_df['Unfilt_a'][i]) * cbm_per_equivalent
        river_graph.nodes[pixel_number][min_water_volume] += (contamination_df['Treat_a'][i] + \
                                                           contamination_df['Unfilt_a'][i]) * excretion\
                                                           * (1 - tertiary_eff) / (0.99 * global_min)
        river_graph.nodes[pixel_number][potential_problem] = True

    cont_above_limit = []
    pixel_number_above_limit = []
    unfeasibility = False  # will keep track if the basin is feasible
    for n in pixel_list:  # for all river pixels
            # Add the contamination of the parent cells
        parents = list(river_graph.predecessors(n))
        discharge_parents = 0
        for k in parents:
            river_graph.nodes[n][cont] += river_graph.nodes[k][cont]
            river_graph.nodes[n][min_water_volume] += river_graph.nodes[k][min_water_volume]
            river_graph.nodes[n][add_volume] += river_graph.nodes[k][add_volume]
            discharge_parents += river_graph.nodes[k][dis]
        river_graph.nodes[n][dis] += river_graph.nodes[n][add_volume]


        if river_graph.nodes[n][min_water_volume] > river_graph.nodes[n][dis]:  # point cannot satisfy cont limit
            if river_graph.nodes[n][sewage_volume] < river_graph.nodes[
                n][min_water_volume] - river_graph.nodes[n][dis]: # if even with sewage water can't satisfy, unfeasible
                unfeasibility = True
            # make the problem feasible (even if it is unfeasible)
            river_graph.nodes[n][add_volume] += river_graph.nodes[n][min_water_volume] - river_graph.nodes[n][dis]
            river_graph.nodes[n][dis] = river_graph.nodes[n][min_water_volume]

        if river_graph.nodes[n][dis] < discharge_parents:  # additional check to ensure conservation of water
            river_graph.nodes[n][dis] = discharge_parents

        # calculate the contamination and relative contamination
        river_graph.nodes[n][cont] *= math.exp(-attenuation * river_graph.nodes[n][RT])
        river_graph.nodes[n][rel_cont] = river_graph.nodes[n][cont] / river_graph.nodes[n][dis]

        # if contamination is over the limit, see if it is a focus-node
        if river_graph.nodes[n][rel_cont] > limit_cont:
            node_problem = True

            # if the node is not the source of pollution, it is not a problematic node
            if not river_graph.nodes[n][potential_problem]:
                node_problem = False

            if node_problem:  # if node is problematic, add it to the list of nodes to be solved
                cont_above_limit.append(river_graph.nodes[n][rel_cont])
                pixel_number_above_limit.append(n)

    # Here some leeway is provided. Instead of solving the problematic points, we give a certain pixel amount of leeway.
    # For isntance, the source node may be contaminated but 2 pixels down the river it has to satisfy the contamination
    # limit.
    number_of_pixels = len(pixel_number_above_limit)
    if children_number > 0:
        for j in range(number_of_pixels):
            # in reverse order such that elements in the list can be deleted without affecting the next index
            i = number_of_pixels - 1 - j
            child = pixel_number_above_limit[i]
            try:
                for m in range(children_number):
                    children = list(river_graph.successors(child))
                    child = children[0]
            except IndexError:
                continue

            contamination_value = river_graph.nodes[child][rel_cont]
            if contamination_value > limit_cont:
                pixel_number_above_limit[i] = child  # set the progeny as the pixel above limit
                cont_above_limit[i] = contamination_value
            else:
                del pixel_number_above_limit[i]
                del cont_above_limit[i]

    return cont_above_limit, pixel_number_above_limit, unfeasibility

def calc_cost(current_treatment: int, end_treatment: int, person_equivalents: float, countryID: int, cost_df:
              pandas.DataFrame) -> float:
    """
       This code takes an upgrade and calculates the cost of that upgrade.
       :rtype: float
       :current_treatment: int: the current level of treatment at the discharge point, excluding the upgrade.
       :contamination_df: pandas.DataFrame: a list of nodes (integer ids) that need to be solved. The overarching idea is that decreasing
       contamination in these nodes solves contamination for the whole basin network.
       :contamination_df: pandas.DataFrame: a dataframe that contains information on all the discharge points for the
       current basin. The dataframe should not contain discharge points that do not belong to the basin.
       :parameters: list: the list of the parameters. they are [excretion, attenuation, filt_efficacy, primary_efficacy,
       secondary_efficacy, tertiary_efficacy]
       :cost_df: A dataframe that contains information on the cost of the upgrades.
       :scenario_number: A string that contains the scenario that is under consideration.
       :return: for each node a list of benefits of upgrading each treatment plant. a list of costs for upgrading each
       treatment plant. The (water) discharges for each point. A changed contamination_df that contains upgrades instead
       of treatment plants.
       """
    upgrade_rows = [(countryID -1) * 3 + i - 1 for i in range(current_treatment, end_treatment)]
    cost = 0
    for upgrade_row in upgrade_rows:
        cost += cost_df.iloc[upgrade_row, 0] * math.pow(person_equivalents, cost_df.iloc[upgrade_row, 1])
        cost += cost_df.iloc[upgrade_row, 2] * math.pow(person_equivalents, cost_df.iloc[upgrade_row, 3])
        cost += cost_df.iloc[upgrade_row, 4] * math.pow(person_equivalents, cost_df.iloc[upgrade_row, 5])
        cost += cost_df.iloc[upgrade_row, 6] * math.pow(person_equivalents, cost_df.iloc[upgrade_row, 7])
    return cost


def upgrade_length_cont_benefit(graph_location: object, upgrade_location: object,
                          model_parameters: list, scenario: str = '') -> None:
    """
     This function takes a dataframe with upgraded treatment plants and calculates for all the treatment plants the
     benefit and adds this to the dataframe. The benefit is defined as the fall in the contamination times the length.
     :rtype: None; output is saved.
     :graph_location: string or networkx.DiGraph: location of the river graph or the river graph DiGraph object.
     :upgrade_location: string: the location of the upgrade dataframe on the computer.
     :model_parameters: list: Contains all the parameters that determine the model, such as excretion, attenuation etc.
     :scenario: string: gives the name of the scenario. The standard scenario is that of hydrorivers.
     :return: A list that contains a sub-graph of the supplied graph as the first element. The second element is a list
     of attributes that are included in the sub-graph.
     """
    # Load inputs
    if isinstance(graph_location, str):
        open_graph = open(graph_location, "rb")
        river_graph = pickle.load(open_graph)
        open_graph.close()
    elif isinstance(graph_location, networkx.DiGraph):
        river_graph = graph_location
    else:
        raise TypeError('Pass either a Digraph or a string. your input was ' + str(type(graph_location)))

    if isinstance(upgrade_location, str):
        contamination_df = pandas.read_csv(upgrade_location)
    elif isinstance(contamination_location, pandas.DataFrame):
        contamination_df = upgrade_location
    else:
        print("the contamination location object was neither a string nor a dataframe")
        raise ValueError

    # extract the downstream of the plant

    excretion, k, primary, secondary, tertiary, filt = model_parameters
    treated_adapted = treatment_plant['Treat_a'] + treatment_plant['Unfilt_a']
    treatment_efficacies = [0, primary, secondary, tertiary]
    current_level = treatment_plant['Treatment_level']
    next_level = current_level + treatment_plant['upgrade']
    treatment_difference = 0 * treatment_plant['upgrade']  # gives correct dimensions immediately
    for i in range(len(treatment_plant)):
        treatment_difference[i] = treatment_efficacies[next_level[i]][i] -\
                                                     treatment_efficacies[current_level[i]][i]
    saved_contamination_list = treated_adapted * treatment_difference

    if scenario == '':
        RT = "RT_HR"
    else:
        RT = "RT " + scenario_number

    treatment_plant['length_cont_benefit'] = 0
    for i in range(len(treatment_plant)):
        saved_contamination = saved_contamination_list[i]
        next_cell = [treatment_plant['pixel_number'][i]] # list that contains next cells
        while next_cell:
            next_cell = next_cell[0]
            saved_contamiantion *= math.exp(-k * river_graph[next_cell][RT])
            length_cont_benefit += saved_contamiantion * river_graph[next_cell]['cell_distance'] / 1000
        treatment_plant['log_length_cont_benefit'][i] = math.log(length_cont_benefit)
    treatment_plant.to_csv(upgrade_location)
    pass


def load_selected_attributes_graph(graph_location: str, attribute_names: list) -> networkx.DiGraph:
    """
     This functions loads a graph from a file location. The graph loaded contains only the 'attribute_names' of
      interest.
     :rtype: networkx.DiGraph
     :graph_location: string: location of the river graph on the computer.
     :attribute_names: list: list of strings that are names of the attributes that need to be in the returned graph.
     :return: networkx.DiGraph: A copy of the river graph with only the attributes of interest.
     """

    # load the supplied graph
    open_graph = open(graph_location, "rb")
    graph = pickle.load(open_graph)
    open_graph.close()

    # make a copy of the structure of the graph.
    copy = networkx.DiGraph()
    copy.add_nodes_from(graph)
    copy.add_edges_from(graph.edges)
    # put the attributes of interest into the copied graph.
    for attribute_name in attribute_names:
        attribute = networkx.get_node_attributes(graph, attribute_name)
        if not attribute:
            graph = None
            raise AttributeError('the attribute is not in the graph')
        networkx.set_node_attributes(copy, attribute, attribute_name)

    graph = None
    return copy


def save_to_main_graph(copy_graph: networkx.DiGraph, graph_location: str, attribute_names: list) -> None:
    """
     This functions takes a copy of the main graph and saves attribute fields from the copy to the main graph.
     :rtype: None
     :copy_graph: The copy of the main graph that contains some new attribute fields.
     :graph_location: string: location of the main river graph on the computer.
     :attribute_names: list: list of strings that are names of the attributes that need to be in the main graph.
     :return: saves the new attributes into the main graph.
     """

    # load the main graph
    open_graph = open(graph_location, "rb")
    main_graph = pickle.load(open_graph)
    open_graph.close()


    # put the attributes of interest into the main graph.
    for attribute_name in attribute_names:
        attribute = networkx.get_node_attributes(copy_graph, attribute_name)
        networkx.set_node_attributes(main_graph, attribute, attribute_name)

    # saving the output
    open_file = open(graph_location, "wb")
    pickle.dump(main_graph, open_file)
    open_file.close()

    return None

def run_model(river_graph: networkx.DiGraph, sorted_river_list: list, contamination_df: pandas.DataFrame,
              parameters: list, scenario: str ='', output_field_name: str = '')\
        -> networkx.DiGraph:
    """
     This functions takes a graph, a contamination dataframe and parameters to calculate contamination in the river
     network.
     :rtype: networkx.DiGraph
     :river_graph: networkx.DiGraph: The river network for which the contamination is calculated
     :sorted_river_list: list: An ordering of the river network according to river hierarchy.
     :contamination_df: pandas.DataFrame: A dataframe of contamination containing the discharge points.
     :parameters: list: A list of the parameters in the model, including excretion, attenuation, filtered efficacy,
                        primary efficacy, secondary efficacy and tertiary efficacy.
     :scenario: string: the scenario to be run. The default gives the hydroRIVERS scenario.
     :output_field_name: string: the output field of the graph. If this field is left empty, the output field will be
     set by default.
     :return: networkx.DiGraph: the river graph with two additional fields, containing the contamination and
                                contaminant concentration.
     """
    excretion, attenuation, filt_eff, primary_eff, secondary_eff, tertiary_eff = parameters

    if scenario == '':
        cont = "Contaminant"
        RT = "RT_HR"
        dis = "flow_HR"
        rel_cont = "Relative Contaminant"
    else:
        cont = "Contaminant " + scenario
        RT = "RT_" + scenario
        dis = scenario
        rel_cont = "Relative contaminant " + scenario

    if output_field_name != '':
        cont = "Contaminant " + output_field_name
        rel_cont = "Relative contaminant " + output_field_name

    # contamination
    treatment_efficacy = primary_eff * (contamination_df["Treatment_level"] == 1) + secondary_eff * \
                     (contamination_df["Treatment_level"] == 2) + tertiary_eff * \
                     (contamination_df["Treatment_level"] == 3)
    # formula
    contamination = (1 - treatment_efficacy) * contamination_df["Treat_a"] + (1 - filt_eff) \
                * contamination_df["Filt_a"] + contamination_df["Unfilt_a"]
    contamination *= contamination_df["pollution"] * excretion
    location = contamination_df["pixel_number"]

    networkx.set_node_attributes(river_graph, 0, name=cont)
    networkx.set_node_attributes(river_graph, 0, name=rel_cont)

    for i in range(len(contamination_df)):
        pixel_number = location[i]
        river_graph.nodes[pixel_number][cont] += contamination[i]

    for n in sorted_river_list:  # for all river pixels
            # Add the contamination of the parent cells
        parents = list(river_graph.predecessors(n))
        for k in parents:
            river_graph.nodes[n][cont] += river_graph.nodes[k][cont]

        river_graph.nodes[n][cont] *= math.exp(-attenuation * river_graph.nodes[n][RT])
        river_graph.nodes[n][rel_cont] = river_graph.nodes[n][cont] / river_graph.nodes[n][dis]

    return river_graph

def graph_to_df_with_tp(river_graph: networkx.DiGraph, topological_sort: list, attributes: list):
    """
     This code takes a river graph and converts it to a dataframe (WITH topological information).
     :rtype: pandas.DataFrame
     :river_graph: networkx.DiGraph or str: the river graph or the location on the computer
     :topological_sort: list of ints: list of pixel numbers that represent the hierarchy of the graph
     :attributes: list of str: the attributes that you wish to be in the dataframe.
     :return: dataframe with information on the nodes of the river graph. Contains topological information.
     """
    first_row = (len(attributes) + 2) * [0]
    graph_df = pandas.DataFrame(first_row)
    graph_df = graph_df.transpose()
    graph_df.columns = ['pixel_number'] + attributes + ['parents']
    graph_df['parents'] = graph_df['parents'].astype(object)
    for i in range(len(topological_sort)):
        new_row = []
        new_row.append(topological_sort[i])
        for attribute in attributes:
            new_row.append(river_graph.nodes[topological_sort[i]][attribute])
        parents = river_graph.predecessors(topological_sort[i])
        new_row.append(list(parents))
        graph_df.loc[i] = new_row
    graph_df = graph_df.set_index('pixel_number')

    return graph_df

def simple_load(file_location):
    open_file = open(file_location, "rb")
    file = pickle.load(open_file)
    open_file.close()
    return file

def create_basin_lists(river_graph, ordered_basins_location='', save=True):
    if ordered_basins_location =='' and save:
        raise ValueError('If you want to save, specify the saving location')
    # creates a list of basins and their ids.
    basins = []
    basins_id = []
    tp = list(networkx.topological_sort(river_graph))
    try:
        river_graph.nodes[tp[0]]['basin']
    except:
        raise KeyError('The river_graph does not contain basin as a field. This might be because you did not load it')
    for n in tp:
        basin_value = river_graph.nodes[n]['basin']
        try:
            index = basins_id.index(basin_value)
            basins[index].append(n)
        except ValueError:
            basins_id.append(basin_value)
            basins.append([n])

    if save:
        open_ts = open(ordered_basins_location, "wb")
        pickle.dump([basins, basins_id], open_ts)
        open_ts.close()
    return basins, basins_id

def discharge_per_basin(basins, contamination_df, discharge_per_basin_location='', save=True):
    if discharge_per_basin_location =='' and save:
        raise ValueError('If you want to save, specify the saving location')

    WWTP_per_basin = [list(pandas.merge(pandas.DataFrame(basins[i]), contamination_df, left_on=0,
                                        right_on='pixel_number')["pixel_number"]) for i in range(len(basins))]
    if save:
        open_ts = open(discharge_per_basin_location, "wb")
        pickle.dump(WWTP_per_basin, open_ts)
        open_ts.close()
    return WWTP_per_basin
