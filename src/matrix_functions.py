import numpy
from scipy import sparse
import networkx
import pandas
import pickle


def sub_basins(river_graph: networkx.DiGraph, ordered_basin_list: list, cut_size: int) -> list:
    """
     this code divides a set of basins into smaller sub-basins. If a basin is split, a connecting cell is shared between
     the upper basin and the lower basin.
     :rtype: list of ints
     :river_graph: networkx.DiGraph: the river graph or the location on the computer
     :ordered_basin_list: list of ints: list of pixel numbers that represent a basin. The basin must be ordered.
     :cut_size: float: maximum size of a sub-basin. (Sometimes this is exceeded if the child cell of the sub-basin would
     have two parents).
     :return: list of sub-basins.
     """
    # preparation
    networkx.set_node_attributes(river_graph, 0, name='parent_count')
    networkx.set_node_attributes(river_graph, 0, name='mini_basin')
    current_basin = 0
    mini_basins = []
    sorted_river_list = []

    # create a topological sort for all the basins
    for basin in ordered_basin_list:
        sorted_river_list += basin

    # create a field within the graph that gives the amount of ancestors that a cell has in a sub-basin.
    for node in sorted_river_list:
        parents = list(river_graph.predecessors(node))
        for parent in parents:
            river_graph.nodes[node]['parent_count'] += river_graph.nodes[parent]['parent_count'] + 1

        # reset the count if cut size is exceeded. (here we will split the basins)
        if river_graph.nodes[node]['parent_count'] > cut_size:
            parents = list(parents)
            if len(parents) == 1:
                river_graph.nodes[node]['parent_count'] = -1

    # Now determine the sub-basin that the pixels belong to.
    sorted_river_list.reverse()
    current_basin = 0
    for node in sorted_river_list:  # note that we are running in reverse. Children come before parents.
        children = list(river_graph.successors(node))

        # if no cut, set basin parent cell to that off child cell.
        if children and river_graph.nodes[node]['parent_count'] > -1:
            river_graph.nodes[node]['mini_basin'] = river_graph.nodes[children[0]]['mini_basin']
        else:  # either if no children, or there is a cut, then set the basin to a new basin id and increment the count.
            river_graph.nodes[node]['mini_basin'] = current_basin
            current_basin += 1

    # collect the nodes into the mini-basins lists.
    sorted_river_list.reverse()
    mini_basins_nodes = [None] * current_basin
    for node in sorted_river_list:

        if mini_basins_nodes[river_graph.nodes[node]['mini_basin']] is None:  # if mini_basin collection undefined yet
            # define the collection {nodes, child, index of child). child and index of child may remain empty.
            mini_basins_nodes[river_graph.nodes[node]['mini_basin']] = [[], None, None]

        mini_basins_nodes[river_graph.nodes[node]['mini_basin']][0].append(node)  # add node to the collection

        # and ensure the child if it is a cut basin, is in the collection.
        child = list(river_graph.successors(node))
        if child:
            child = child[0]
            if river_graph.nodes[node]['mini_basin'] != river_graph.nodes[child]['mini_basin']:
                mini_basins_nodes[river_graph.nodes[node]['mini_basin']][2] = sorted_river_list.index(child)
                mini_basins_nodes[river_graph.nodes[node]['mini_basin']][1] = child
                mini_basins_nodes[river_graph.nodes[node]['mini_basin']][0].append(child)

    # now join small basins
    j = 0
    while j < len(mini_basins_nodes) - 1:
        current_mini_basin_len = len(mini_basins_nodes[j][0])
        next_mini_basin_len = len(mini_basins_nodes[j+1][0])
        if current_mini_basin_len + next_mini_basin_len < cut_size:  # if basins can be joined without exceeding size
            mini_basins_nodes[j + 1][0] = mini_basins_nodes[j][0] + mini_basins_nodes[j+1][0]  # join them.
            del mini_basins_nodes[j]
        else:
            j += 1

    mini_basins_nodes.reverse()

    return mini_basins_nodes


def graph_to_RT_matrix(river_graph: networkx.DiGraph, ordered_basins: list, cut_size: int, output_location: str = '')\
        -> list:
    """
     This code takes a river_graph and creates the matrix A. The matrix A solves the model according to the matrix
     equation y= Ax, where x is the initial state of the model and y is the result. The rows of A represent the order
     that is given by the sorted_river_list. The sorted_river list must be a topological sort. The inputs x and y must
     have rows (e.g. entries since they are vectors) that represent the same nodes as the sorted_river list.
     :rtype: numpy.array
     :river_graph: networkx.DiGraph: the river graph or the location on the computer
     :topological_sort: list of ints: list of pixel numbers that represent the hierarchy of the graph
     :attenuation: float: the speed by which contaminants degrade per hour
     :return: returns the matrix A, that solves the model with equation y = Ax, where x is the initial state and y is
     the result.
     """
    mini_basins_nodes = sub_basins(river_graph, ordered_basins, cut_size)
    basin_matrices = []
    for basin in mini_basins_nodes:
        basin_ids = set()
        basin_nodes = basin[0]
        B = numpy.zeros((len(basin_nodes), len(basin_nodes)), dtype=numpy.float32)

        for j in reversed(range(len(basin_nodes))):  # the very last one can only have a false child!
            node = basin_nodes[j]
            RT = river_graph.nodes[node]["RT_HR"]
            child = list(river_graph.successors(node))

            # this recursion determines the residence time between all the children and the parent j
            if child and j != len(basin_nodes) - 1:
                child_index = basin_nodes.index(child[0], j)
                column = B[:, child_index]
                B[:, j] = (column > 0) * RT + column
            # 0 is reserved for no connection. If there is an RT of 0 and there is a connection (with lakes), then
            # we take a value very close to 0 instead.
            if RT == 0:
                RT = 0.000001
            B[j, j] = RT
            basin_ids.add(river_graph.nodes[node]["basin"])
        df = pandas.DataFrame(basin_nodes)
        df['indicators'] = 0
        df = df.set_index(0)

        basin_matrices.append([B, basin[1], df, basin_ids])
    if output_location != '':
        open_ts = open(output_location, "wb")
        pickle.dump(basin_matrices, open_ts)
        open_ts.close()

    return basin_matrices


def matrix_subset(matrices_source, basin_ids_picked, basin_nodes, basin_ids_locations, cut_size):
    # open the matrices
    open_ts = open(matrices_source, "rb")
    basin_matrices = pickle.load(open_ts)
    open_ts.close()

    # step 1; clean space
    for j in reversed(range(len(basin_matrices))):
        lets_delete = True
        for basin_id in basin_matrices[j][3]:
            if basin_id in basin_ids_picked:
                lets_delete = False
                break
        if lets_delete:
            del basin_matrices[j]

    # step 2, reduce matrices that are only partly necessary
    for j in range(len(basin_matrices)):
        delete_some = False
        for basin_id in basin_matrices[j][3]:
            if basin_id not in basin_ids_picked:
                delete_some = True
                break
        if delete_some:
            nodes_to_keep = []
            for basin_id in basin_matrices[j][3]:
                if basin_id in basin_ids_picked:
                    basin_id_index = basin_ids_locations.index(basin_id)
                    nodes_to_keep += basin_nodes[basin_id_index]
            try:
                basin_matrices[j][2]['indicators'].loc[nodes_to_keep] = 1
                indices = numpy.array(basin_matrices[j][2]['indicators'], bool)
                basin_matrices[j][0] = basin_matrices[j][0][indices, :]
                basin_matrices[j][0] = basin_matrices[j][0][:, indices]
                basin_matrices[j][2] = basin_matrices[j][2][indices == 1]
                basin_matrices[j][2]['indicators'] = 0

            except KeyError:  # should occur if the entire matrix must be kept
                pass

    # step 3, connect small matrices within same basin

    # collect the pixels that the matrices refer to
    final_row = [0]  # this is a list of the final row number of the matrices (actually starting)
    df = basin_matrices[0][2]  # this contains all the pixels in order of the matrices

    for j in range(len(basin_matrices)):
        if basin_matrices[j][1] is not None:
            basin_matrices[j][2].drop(basin_matrices[j][1], axis=0, inplace=True)
        final_row.append(final_row[len(final_row) - 1] + len(basin_matrices[j][2]))
        if j > 0:
            df = pandas.concat([df, basin_matrices[j][2]], axis=0)

    nodes_df = df.reset_index()
    node_list = list(nodes_df[0])  # now we have the pixels in the order of the matrices

    #  We continue to try to merge matrices that are in the same basin
    merge_found = True
    while merge_found:
        merge_found = False
        j = 0
        while j < len(basin_matrices):  # while statement since basin_matrices may decrease in length.
            if basin_matrices[j][1] is not None:  # only consider a merge if the basin is a sub-basin
                child_index = node_list.index(basin_matrices[j][1])
                # find the matrix that contains this child_index
                for i in range(len(final_row)):
                    if child_index > final_row[i]:
                        index = i
                current_matrix = basin_matrices[j][0]
                next_matrix = basin_matrices[index][0]
                size_current = len(current_matrix) - 1
                size_next = len(next_matrix)
                if size_current + size_next < cut_size:  # merge the matrices if the cut_size is acceptable
                    merge_found = True
                    # last_row current forms the link between the two matrices.
                    last_row_current = current_matrix[size_current][0:size_current] - \
                                       current_matrix[size_current][size_current]
                    # delete last row and column from the current matrix
                    current_matrix = current_matrix[:, 0:size_current]
                    current_matrix = current_matrix[0:size_current][:]

                    # now create the new matrix that links the sub-basins
                    last_row_indicator = last_row_current > 0
                    C = numpy.zeros([size_next, size_current])
                    D = numpy.zeros([size_current, size_next])
                    child_index_next_matrix = child_index - final_row[index]
                    for row_number in range(size_next):
                        C[row_number][:] = last_row_current * (next_matrix[row_number][child_index_next_matrix] > 0) \
                                           + next_matrix[row_number][child_index_next_matrix]
                    new_matrix = numpy.asarray(numpy.bmat([[current_matrix, D], [C, next_matrix]]))
                    basin_matrices[index][0] = new_matrix
                    basin_matrices[index][2] = pandas.concat([basin_matrices[j][2], basin_matrices[index][2]], axis=0)
                    new_node_list = node_list[0:final_row[j]]
                    new_node_list += node_list[final_row[j + 1]:final_row[index]]
                    new_nodes = basin_matrices[index][
                        2].reset_index()  # the index is one to the left due to delete of j
                    new_node_list += list(new_nodes[0])
                    new_node_list += node_list[final_row[index + 1]:len(node_list)]
                    node_list = new_node_list
                    del basin_matrices[j]
                    del final_row[j]
                    final_row[j:index] = [final_row[n] - size_current for n in range(j, index)]

                else:
                    j += 1
            else:
                j += 1

    # step 3; rejoin some very small matrices.
    j = 0
    while j < len(basin_matrices) - 1:
        current_mini_basin_len = len(basin_matrices[j][0])
        next_mini_basin_len = len(basin_matrices[j + 1][0])
        if basin_matrices[j][1] is None:
            if current_mini_basin_len + next_mini_basin_len < cut_size:
                Y = numpy.zeros([current_mini_basin_len, next_mini_basin_len])
                Z = numpy.zeros([next_mini_basin_len, current_mini_basin_len])
                basin_matrices[j + 1][0] = numpy.asarray(
                    numpy.bmat([[basin_matrices[j][0], Y], [Z, basin_matrices[j + 1][0]]]))
                basin_matrices[j + 1][2] = pandas.concat([basin_matrices[j][2], basin_matrices[j + 1][2]], axis=0)
                del basin_matrices[j]
                j -= 1
        j +=1


    df = basin_matrices[0][2]  # this contains all the pixels in order of the matrices

    for j in range(1, len(basin_matrices)):
        df = pandas.concat([df, basin_matrices[j][2]], axis=0)

    return basin_matrices, df


def create_attenuation_matrices(basin_matrices, attenuation):
    for i in range(len(basin_matrices)):
        rows, columns = numpy.shape(basin_matrices[i][0])
        for j in range(rows):
           current_row = basin_matrices[i][0][j,:]
           basin_matrices[i][0][j,:] = (current_row > 0) * numpy.exp(-attenuation * current_row)
    return basin_matrices


def run_basin_matrices(basin_matrices, pixel_order, initial_contamination):
    final_contamination = numpy.zeros([len(initial_contamination)])
    indicators = numpy.zeros([len(initial_contamination)], bool)
    row_count = 0

    for j in range(len(basin_matrices)):
        rows, columns = numpy.shape(basin_matrices[j][0])
        last_row_count = row_count
        row_count += rows
        if basin_matrices[j][1] is None:
            final_contamination[last_row_count:row_count] = \
                numpy.matmul(basin_matrices[j][0], initial_contamination[last_row_count:row_count])
        else:
            row_count -= 1
            child_index = pixel_order.index(basin_matrices[j][1], row_count)
            indicators[last_row_count:row_count] = True
            indicators[child_index] = True
            final_contamination[indicators] = numpy.matmul(basin_matrices[j][0], initial_contamination[indicators])
            indicators = numpy.zeros([len(initial_contamination)], bool)
            initial_contamination[child_index] = final_contamination[child_index] / basin_matrices[j][0][
                rows - 1, rows - 1]

    return final_contamination


def matrix_subset_sparse(matrices_source, basin_ids_picked, basin_nodes, basin_ids_locations, cut_size, opt='rows',
                         cut_minimum=0):
    # open the matrices
    if isinstance(matrices_source, str):
        open_ts = open(matrices_source, "rb")
        basin_matrices = pickle.load(open_ts)
        open_ts.close()
    else:
        basin_matrices = matrices_source

    # step 1; clean space
    for j in reversed(range(len(basin_matrices))):
        lets_delete = True
        for basin_id in basin_matrices[j][3]:
            if basin_id in basin_ids_picked:
                lets_delete = False
                break
        if lets_delete:
            del basin_matrices[j]

    # step 2, reduce matrices that are only partly necessary
    for j in range(len(basin_matrices)):
        delete_some = False
        for basin_id in basin_matrices[j][3]:
            if basin_id not in basin_ids_picked:
                delete_some = True
                break
        if delete_some:
            nodes_to_keep = []
            for basin_id in basin_matrices[j][3]:
                if basin_id in basin_ids_picked:
                    basin_id_index = basin_ids_locations.index(basin_id)
                    nodes_to_keep += basin_nodes[basin_id_index]
            try:
                basin_matrices[j][2]['indicators'].loc[nodes_to_keep] = 1
                indices = numpy.array(basin_matrices[j][2]['indicators'], bool)
                basin_matrices[j][0] = basin_matrices[j][0][indices, :]
                basin_matrices[j][0] = basin_matrices[j][0][:, indices]
                basin_matrices[j][2] = basin_matrices[j][2][indices == 1]
                basin_matrices[j][2]['indicators'] = 0

            except KeyError:  # should occur if the entire matrix must be kept
                pass

    # step 3, connect small matrices within same basin

    # collect the pixels that the matrices refer to
    final_row = [0]  # this is a list of the final row number of the matrices (actually starting)
    df = basin_matrices[0][2]  # this contains all the pixels in order of the matrices

    for j in range(len(basin_matrices)):
        if basin_matrices[j][1] is not None:
            basin_matrices[j][2].drop(basin_matrices[j][1], axis=0, inplace=True)
        final_row.append(final_row[len(final_row) - 1] + len(basin_matrices[j][2]))
        if j > 0:
            df = pandas.concat([df, basin_matrices[j][2]], axis=0)

    nodes_df = df.reset_index()
    node_list = list(nodes_df[0])  # now we have the pixels in the order of the matrices

    #  We continue to try to merge matrices that are in the same basin
    merge_found = True
    while merge_found:
        merge_found = False
        j = 0
        while j < len(basin_matrices):  # while statement since basin_matrices may decrease in length.
            if isinstance(basin_matrices[j][0], numpy.ndarray):
                if len(basin_matrices[j][0]) > cut_minimum:
                    if opt=='rows':
                        basin_matrices[j][0] =  sparse.csr_matrix(basin_matrices[j][0])
                    elif opt=='columns':
                        basin_matrices[j][0] =  sparse.csc_matrix(basin_matrices[j][0])
                    j += 1
                    continue

            if basin_matrices[j][1] is not None:  # only consider a merge if the basin is a sub-basin
                child_index = node_list.index(basin_matrices[j][1])
                # find the matrix that contains this child_index
                for i in range(len(final_row)):
                    if child_index > final_row[i]:
                        index = i
                current_matrix = basin_matrices[j][0]
                next_matrix = basin_matrices[index][0]
                if not isinstance(next_matrix, numpy.ndarray) or not isinstance(current_matrix, numpy.ndarray):
                    j += 1
                    continue
                size_current = len(current_matrix) - 1
                size_next = len(next_matrix)
                if size_current + size_next < cut_size:  # merge the matrices if the cut_size is acceptable
                    merge_found = True
                    # last_row current forms the link between the two matrices.
                    last_row_current = current_matrix[size_current][0:size_current] - \
                                       current_matrix[size_current][size_current]
                    # delete last row and column from the current matrix
                    current_matrix = current_matrix[:, 0:size_current]
                    current_matrix = current_matrix[0:size_current][:]

                    # now create the new matrix that links the sub-basins
                    last_row_indicator = last_row_current > 0
                    C = numpy.zeros([size_next, size_current])
                    D = numpy.zeros([size_current, size_next])
                    child_index_next_matrix = child_index - final_row[index]
                    for row_number in range(size_next):
                        C[row_number][:] = last_row_current * (next_matrix[row_number][child_index_next_matrix] > 0) \
                                           + next_matrix[row_number][child_index_next_matrix]
                    new_matrix = numpy.asarray(numpy.bmat([[current_matrix, D], [C, next_matrix]]))
                    basin_matrices[index][0] = new_matrix
                    basin_matrices[index][2] = pandas.concat([basin_matrices[j][2], basin_matrices[index][2]], axis=0)
                    new_node_list = node_list[0:final_row[j]]
                    new_node_list += node_list[final_row[j + 1]:final_row[index]]
                    new_nodes = basin_matrices[index][
                        2].reset_index()  # the index is one to the left due to delete of j
                    new_node_list += list(new_nodes[0])
                    new_node_list += node_list[final_row[index + 1]:len(node_list)]
                    node_list = new_node_list
                    del basin_matrices[j]
                    del final_row[j]
                    final_row[j:index] = [final_row[n] - size_current for n in range(j, index)]

                else:
                    j += 1
            else:
                j += 1

    # step 4; rejoin some very small matrices.
    while j < len(basin_matrices) - 1:
        current_mini_basin_len = len(basin_matrices[j][0])
        next_mini_basin_len = len(basin_matrices[j+1][0])
        if basin_matrices[j][1] is None:
            if current_mini_basin_len + next_mini_basin_len < cut_size:
                Y = numpy.zeros([current_mini_basin_len, next_mini_basin_len])
                Z = numpy.zeros([next_mini_basin_len, current_mini_basin_len])
                basin_matrices[j+1][0] = numpy.asarray(
                    numpy.bmat([[basin_matrices[j][0], Y], [Z, basin_matrices[j+1][0]]]))
                basin_matrices[j+1][2] = pandas.concat([basin_matrices[j][2], basin_matrices[j+1][2]], axis=0)
                del basin_matrices[j]
                j -= 1
        j +=1


    df = basin_matrices[0][2]  # this contains all the pixels in order of the matrices

    for j in range(1, len(basin_matrices)):
        df = pandas.concat([df, basin_matrices[j][2]], axis=0)

    # make the matrices sparse:
    if opt == 'rows':
        for j in range(len(basin_matrices)):
            if isinstance(basin_matrices[j][0], numpy.ndarray):
                basin_matrices[j][0] = sparse.csr_matrix(basin_matrices[j][0])
    elif opt == 'columns':
        for j in range(len(basin_matrices)):
            if isinstance(basin_matrices[j][0], numpy.ndarray):
               basin_matrices[j][0] = sparse.csc_matrix(basin_matrices[j][0])

    return basin_matrices, df


def create_attenuation_matrices_sparse(basin_matrices, attenuation):
    for i in range(len(basin_matrices)):
        basin_matrices[i][0].data = numpy.exp(-attenuation * basin_matrices[i][0].data)

    return basin_matrices


def run_basin_matrices_sparse(basin_matrices, pixel_order, init_cont):
    initial_contamination = init_cont.copy()
    final_contamination = numpy.zeros([len(initial_contamination)])
    indicators = numpy.zeros([len(initial_contamination)], bool)
    row_count = 0

    for j in range(len(basin_matrices)):
        rows, columns = numpy.shape(basin_matrices[j][0])
        last_row_count = row_count
        row_count += rows
        if basin_matrices[j][1] is None:
            final_contamination[last_row_count:row_count] = \
                basin_matrices[j][0] @ initial_contamination[last_row_count:row_count]
        else:
            row_count -= 1
            child_index = pixel_order.index(basin_matrices[j][1], row_count)
            indicators[last_row_count:row_count] = True
            indicators[child_index] = True
            final_contamination[indicators] = basin_matrices[j][0] @ initial_contamination[indicators]
            indicators = numpy.zeros([len(initial_contamination)], bool)
            initial_contamination[child_index] = final_contamination[child_index] / basin_matrices[j][0][
               rows - 1, rows - 1]

    return final_contamination



def get_plant_column_sparse(basin_matrices, pixel_order, plant_pixel_number):

    plant_column = numpy.zeros([len(pixel_order)])
    indicators = numpy.zeros([len(pixel_order)], bool)
    row_count = 0

    # these values adapt in case the initial point of contamination works through more matrices.
    current_child_index = pixel_order.index(plant_pixel_number)
    multiply_value = 1

    for j in range(len(basin_matrices)):
        rows, columns = numpy.shape(basin_matrices[j][0])
        last_row_count = row_count
        row_count += rows
        if row_count <= current_child_index:
            if basin_matrices[j][1] is not None:
                row_count -= 1
            continue
        column = basin_matrices[j][0][:, current_child_index - last_row_count] * multiply_value
        column = column.toarray()
        column.shape = (rows)
        if basin_matrices[j][1] is None:  # means this is the last iteration.
            plant_column[last_row_count:row_count] = column
            break
        else:
            row_count -= 1
            plant_column[last_row_count:row_count] = column[0:rows-1]
            if column[current_child_index - last_row_count] < 0.00001:
                multiply_value = 0
            else:
                multiply_value *= column[rows-1] / column[current_child_index - last_row_count]
            current_child_index = pixel_order.index(basin_matrices[j][1], row_count)

    return plant_column


def simple_reduce(matrices_source, basin_ids_picked, basin_nodes, basin_ids_locations):
    # open the matrices
    if isinstance(matrices_source, str):
        open_ts = open(matrices_source, "rb")
        basin_matrices, df = pickle.load(open_ts)
        open_ts.close()
    else:
        basin_matrices = matrices_source

    # step 1; clean space
    for j in reversed(range(len(basin_matrices))):
        lets_delete = True
        for basin_id in basin_matrices[j][3]:
            if basin_id in basin_ids_picked:
                lets_delete = False
                break
        if lets_delete:
            del basin_matrices[j]

    # step 2, reduce matrices that are only partly necessary
    for j in range(len(basin_matrices)):
        delete_some = False
        for basin_id in basin_matrices[j][3]:
            if basin_id not in basin_ids_picked:
                delete_some = True
                break
        if delete_some:
            nodes_to_keep = []
            for basin_id in basin_matrices[j][3]:
                if basin_id in basin_ids_picked:
                    basin_id_index = basin_ids_locations.index(basin_id)
                    nodes_to_keep += basin_nodes[basin_id_index]
            try:
                basin_matrices[j][2]['indicators'].loc[nodes_to_keep] = 1
                indices = numpy.array(basin_matrices[j][2]['indicators'], bool)
                basin_matrices[j][0] = basin_matrices[j][0][indices, :]
                basin_matrices[j][0] = basin_matrices[j][0][:, indices]
                basin_matrices[j][2] = basin_matrices[j][2][indices == 1]
                basin_matrices[j][2]['indicators'] = 0

            except KeyError:  # should occur if the entire matrix must be kept
                pass
    df = basin_matrices[0][2]  # this contains all the pixels in order of the matrices

    for j in range(1, len(basin_matrices)):
        df = pandas.concat([df, basin_matrices[j][2]], axis=0)

    return basin_matrices, df


def get_plant_column_sparse(basin_matrices, pixel_order, plant_pixel_number):

    plant_column = numpy.zeros([len(pixel_order)])
    indicators = numpy.zeros([len(pixel_order)], bool)
    row_count = 0

    # these values adapt in case the initial point of contamination works through more matrices.
    current_child_index = pixel_order.index(plant_pixel_number)
    multiply_value = 1

    for j in range(len(basin_matrices)):
        rows, columns = numpy.shape(basin_matrices[j][0])
        last_row_count = row_count
        row_count += rows
        if row_count <= current_child_index:
            if basin_matrices[j][1] is not None:
                row_count -= 1
            continue
        column = basin_matrices[j][0][:, current_child_index - last_row_count] * multiply_value
        column = column.toarray()
        column.shape = (rows)
        if basin_matrices[j][1] is None:  # means this is the last iteration.
            plant_column[last_row_count:row_count] = column
            break
        else:
            row_count -= 1
            plant_column[last_row_count:row_count] = column[0:rows-1]
            if column[current_child_index - last_row_count] < 0.00001:
                multiply_value = 0
            else:
                multiply_value *= column[rows-1] / column[current_child_index - last_row_count]
            current_child_index = pixel_order.index(basin_matrices[j][1], row_count)

    return plant_column


def simple_reduce(matrices_source, basin_ids_picked, basin_nodes, basin_ids_locations):
    # open the matrices
    if isinstance(matrices_source, str):
        open_ts = open(matrices_source, "rb")
        basin_matrices, df = pickle.load(open_ts)
        open_ts.close()
    else:
        basin_matrices = matrices_source

    # step 1; clean space
    for j in reversed(range(len(basin_matrices))):
        lets_delete = True
        for basin_id in basin_matrices[j][3]:
            if basin_id in basin_ids_picked:
                lets_delete = False
                break
        if lets_delete:
            del basin_matrices[j]

    # step 2, reduce matrices that are only partly necessary
    for j in range(len(basin_matrices)):
        delete_some = False
        for basin_id in basin_matrices[j][3]:
            if basin_id not in basin_ids_picked:
                delete_some = True
                break
        if delete_some:
            nodes_to_keep = []
            for basin_id in basin_matrices[j][3]:
                if basin_id in basin_ids_picked:
                    basin_id_index = basin_ids_locations.index(basin_id)
                    nodes_to_keep += basin_nodes[basin_id_index]
            try:
                basin_matrices[j][2]['indicators'].loc[nodes_to_keep] = 1
                indices = numpy.array(basin_matrices[j][2]['indicators'], bool)
                basin_matrices[j][0] = basin_matrices[j][0][indices, :]
                basin_matrices[j][0] = basin_matrices[j][0][:, indices]
                basin_matrices[j][2] = basin_matrices[j][2][indices == 1]
                basin_matrices[j][2]['indicators'] = 0

            except KeyError:  # should occur if the entire matrix must be kept
                pass
    df = basin_matrices[0][2]  # this contains all the pixels in order of the matrices

    for j in range(1, len(basin_matrices)):
        df = pandas.concat([df, basin_matrices[j][2]], axis=0)

    return basin_matrices, df

def get_plant_column_sparse(basin_matrices, pixel_order, plant_pixel_number):

    plant_column = numpy.zeros([len(pixel_order)])
    indicators = numpy.zeros([len(pixel_order)], bool)
    row_count = 0

    # these values adapt in case the initial point of contamination works through more matrices.
    current_child_index = pixel_order.index(plant_pixel_number)
    multiply_value = 1

    for j in range(len(basin_matrices)):
        rows, columns = numpy.shape(basin_matrices[j][0])
        last_row_count = row_count
        row_count += rows
        if row_count <= current_child_index:
            if basin_matrices[j][1] is not None:
                row_count -= 1
            continue
        column = basin_matrices[j][0][:, current_child_index - last_row_count] * multiply_value
        column = column.toarray()
        column.shape = (rows)
        if basin_matrices[j][1] is None:  # means this is the last iteration.
            plant_column[last_row_count:row_count] = column
            break
        else:
            row_count -= 1
            plant_column[last_row_count:row_count] = column[0:rows-1]
            if column[current_child_index - last_row_count] < 0.00001:
                multiply_value = 0
            else:
                multiply_value *= column[rows-1] / column[current_child_index - last_row_count]
            current_child_index = pixel_order.index(basin_matrices[j][1], row_count)

    return plant_column


def simple_reduce(matrices_source, basin_ids_picked, basin_nodes, basin_ids_locations):
    # open the matrices
    if isinstance(matrices_source, str):
        open_ts = open(matrices_source, "rb")
        basin_matrices, df = pickle.load(open_ts)
        open_ts.close()
    else:
        basin_matrices = matrices_source

    # step 1; clean space
    for j in reversed(range(len(basin_matrices))):
        lets_delete = True
        for basin_id in basin_matrices[j][3]:
            if basin_id in basin_ids_picked:
                lets_delete = False
                break
        if lets_delete:
            del basin_matrices[j]

    # step 2, reduce matrices that are only partly necessary
    for j in range(len(basin_matrices)):
        delete_some = False
        for basin_id in basin_matrices[j][3]:
            if basin_id not in basin_ids_picked:
                delete_some = True
                break
        if delete_some:
            nodes_to_keep = []
            for basin_id in basin_matrices[j][3]:
                if basin_id in basin_ids_picked:
                    basin_id_index = basin_ids_locations.index(basin_id)
                    nodes_to_keep += basin_nodes[basin_id_index]
            try:
                basin_matrices[j][2]['indicators'].loc[nodes_to_keep] = 1
                indices = numpy.array(basin_matrices[j][2]['indicators'], bool)
                basin_matrices[j][0] = basin_matrices[j][0][indices, :]
                basin_matrices[j][0] = basin_matrices[j][0][:, indices]
                basin_matrices[j][2] = basin_matrices[j][2][indices == 1]
                basin_matrices[j][2]['indicators'] = 0

            except KeyError:  # should occur if the entire matrix must be kept
                pass
    df = basin_matrices[0][2]  # this contains all the pixels in order of the matrices

    for j in range(1, len(basin_matrices)):
        df = pandas.concat([df, basin_matrices[j][2]], axis=0)

    return basin_matrices, df


def get_plant_column(basin_matrices, pixel_order, plant_pixel_number):

    plant_column = numpy.zeros([len(pixel_order)])
    indicators = numpy.zeros([len(pixel_order)], bool)
    row_count = 0

    # these values adapt in case the initial point of contamination works through more matrices.
    current_child_index = pixel_order.index(plant_pixel_number)
    multiply_value = 1

    for j in range(len(basin_matrices)):
        rows, columns = numpy.shape(basin_matrices[j][0])
        last_row_count = row_count
        row_count += rows
        if row_count <= current_child_index:
            if basin_matrices[j][1] is not None:
                row_count -= 1
            continue
        column = basin_matrices[j][0][:, current_child_index - last_row_count] * multiply_value
        column.shape = (rows)
        if basin_matrices[j][1] is None:  # means this is the last iteration.
            plant_column[last_row_count:row_count] = column
            break
        else:
            row_count -= 1
            plant_column[last_row_count:row_count] = column[0:rows-1]
            if column[current_child_index - last_row_count] < 0.00001:
                multiply_value = 0
            else:
                multiply_value *= column[rows-1] / column[current_child_index - last_row_count]
            current_child_index = pixel_order.index(basin_matrices[j][1], row_count)

    return plant_column

def add_initial_contaminant_to_df(cont_df, pixels_df, all_df, parameters):
    excretion, attenuation, filt_eff, primary_eff, secondary_eff, tertiary_eff = parameters
    pixels_df = pandas.merge(pixels_df, all_df, left_index=True, right_index=True)

    # contamination
    treatment_efficacy = primary_eff * (cont_df["Treatment_level"] == 1) + secondary_eff * \
                         (cont_df["Treatment_level"] == 2) + tertiary_eff * \
                         (cont_df["Treatment_level"] == 3)
    # formula
    cont_df['init_cont'] = (1 - treatment_efficacy) * cont_df["Treat_a"] + (1 - filt_eff) \
                                    * cont_df["Filt_a"] + cont_df["Unfilt_a"]
    cont_df['init_cont'] *= cont_df["pollution"] * excretion
    cont_df = cont_df[['pixel_number', 'init_cont']]
    cont_df = cont_df.groupby('pixel_number').sum()
    cont_df = cont_df.reset_index()
    pixels_df['indicators'][cont_df['pixel_number']] = cont_df['init_cont']
    pixels_df.columns = ['initial_contaminant', 'weight', 'id_Country']
    return pixels_df
