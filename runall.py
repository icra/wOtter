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


