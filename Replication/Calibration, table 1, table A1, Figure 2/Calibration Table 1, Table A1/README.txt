Independent inputs are:
-

Dependent inputs are:
Occurences_in_river.csv
AGG_WWTP_df.csv
river_graph.pkl OR graph_and_sort.pkl*

Output files are:
shapefiles for each contaminant (e.g. Lumped 14.shp, Lumped 14.dbf, Lumped 14.shx)
Rasters for each contaminant (e.g. Lumped 14.tif)
results_calibration.csv
graph_and_sort.pkl (also an input) *

The dependent inputs are provided, but may be replicated.

Occurences_in_river.csv is created in: Calibration, table 1, table A1, Figure 2\Create observed occurences

AGG_WWTP_df.csv, river_graph.pkl, reference_raster.tif are created by installing the model from scratch. 
This requires running "runall.py".

The shapefiles and the rasters are for visualization and may be opened in the program Qgis(among others).
If there is little space on the drive, one may disable these outputs by adjusting the code. In particular,
set the variable suppress_shapefile_raster_creation to 'True'.
The table "results_calibration.csv" is the basis for table 1 and appendix table A1.

* Graph_and_sort.pkl is a reduced version of river_graph.pkl. Once this output is created, when the code
is rerun it will read graph_and_sort.pkl, which greatly increases speed.
