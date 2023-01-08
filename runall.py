import os
from time import time
import numpy
import warnings
import pandas
import zipfile
import shutil
from src.library import shapefile_raster_functions

def run_script(script_to_run, supress_warnings):
    if supress_warnings:
        # silence command-line output temporarily
        warnings.filterwarnings("ignore")
        exec(open(script_to_run).read())
    else:
        exec(open(script_to_run).read())


import networkx as nx
g = nx.read_gpickle("./data/river_graph.pkl")


script_path = "./src/"
script_list = ["1.1 calculate slope.py", "2. adjust_hydrorivers.py", "3. shapefile to graph.py",
               "4.1. AGG_WWTP_to_shapefile.py", "4.2. Create_land_and_WWTP_rasters.py",
               "5.1. Cut_basins.py", "5.2. Zonal_statistics.py", "5.3. Assign_filtered_and_unfiltered.py",
               "5.4. Calculate_multipliers.py", "5.5. Create_adapted_shapefile.py",
               "6.1. adjust discharge.py", "6.2. write_discharge.py", "6.3. write_treatment_levels.py",
               "6.4. adjust_graph.py", "7.1. adjust observation.py", "7.2. create_observed_cont.py"]
last_step_location = 'last_step.npy'
recorded_runtimes_location = 'runtimes.npy'
supress_output = True  # swap for False if you want to see errors or warnings
general_name_extension = 'runtimes_intent'  # name of the maps for results and intermediate outputs
compress_all = False  # compresses even the final outputs if True
intermediate_output = False  # yields intermediate output
include_advanced_output = False
include_basin_matrices = False
advanced_scripts = ['8.1 basins_shapefile.py', '8.2 basins_and_WWTP lists.py']
basin_matrices_script = ['8.3 basin_matrices.py']
if include_advanced_output:
    script_list += advanced_scripts
if include_basin_matrices:
    script_list += basin_matrices_script
try:
    last_step = numpy.load(last_step_location)
    #last_step = numpy.load(last_step_location, allow_pickle=True)
    recorded_runtimes = numpy.load(recorded_runtimes_location)
except FileNotFoundError:
    last_step = 0
    recorded_runtimes = numpy.zeros([len(script_list)])

starting_step = last_step
for script_index in range(starting_step, len(script_list)):
    full_script_location = os.path.join(script_path, script_list[script_index])
    start_time = time()
    print(full_script_location)
    run_script(full_script_location, supress_output)
    recorded_runtimes[script_index] = time() - start_time
    print("step " + str(script_index) + " complete")
    last_step += 1
    numpy.save(last_step_location, last_step)
    numpy.save(recorded_runtimes_location, recorded_runtimes)

runtime_df = pandas.DataFrame(recorded_runtimes)
runtime_df.to_csv('runtimes.csv')

directory = os.path.join(os.getcwd(), 'data')

result_location = os.path.join(os.getcwd(), 'results', general_name_extension)
os.mkdir(result_location)
main_outputs = ['AGG_WWTP_df.csv', 'reference_raster.tif', 'rivers_from_graph.tif', 'river_graph.pkl',
                'sorted_river_list.pkl']
advanced_outputs = ['river_basins.shx', 'river_basins.prj', 'river_basins.dbf', 'river_basins.shp',
                    'ordered_basins.pkl', 'WWTP_per_basin.pkl']
matrix_output = ['basin_matrices.pkl']

if include_advanced_output:
    main_outputs += advanced_outputs
if include_basin_matrices:
    main_outputs += matrix_output
for output in main_outputs:
    root_location = directory + '/' + output
    result_name = result_location + '/' + output

if compress_all:
    zipfile_name = os.join(result_location, 'final_results.zip')
    compression = zipfile.ZIP_DEFLATED
    zf = zipfile.ZipFile(zipfile_name, mode="w")
    for output_name in os.listdir(result_location):
        complete_output_name = os.path.join()
        zf.write(complete_output_name, complete_output_name, compress_type=compression)
    zf.close()

intermediate_files = os.listdir(directory)
"""
intermediate_files.remove('Raw data')
if intermediate_output:
    zipfile_name = 'zipped_intermediate_results'
    compression = zipfile.ZIP_DEFLATED
    zf = zipfile.ZipFile(zipfile_name, mode="w")
    for file in intermediate_files:
        complete_file_name = os.path.join(directory, file)
        zf.write(complete_file_name, file, compress_type=compression)
    zf.close()
    os.rename(zipfile_name, os.path.join(result_location, 'intermediate_results.zip'))

for file in intermediate_files:
    complete_file_name = os.path.join(directory, file)
    try:
        os.remove(complete_file_name)
    except:
        try:
            shutil.rmtree(complete_file_name)  # works if a directory
        except:
            pass  # will occur if the file does not exist or is in use

"""
