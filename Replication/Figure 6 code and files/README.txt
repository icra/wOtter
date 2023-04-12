Independent inputs are:
hybas_eu_lev06_v1c.shp
hybas_eu_lev06_v1c.dbf
hybas_eu_lev06_v1c.prj
hybas_eu_lev06_v1c.sbn
hybas_eu_lev06_v1c.sbx
hybas_eu_lev06_v1c.shp.xml
hybas_eu_lev06_v1c.shx

Dependent inputs are:
AGG_WWTP_df.csv
reference_raster.tif
river_graph.pkl
sorted_river_list.pkl

Outputs are:
contamination_pfafstetter_basin.shp

The dependent inputs are provided, but may be replicated.
All dependent inputs are created by installing the model. This may be done by running "runall.py"

The output is a shapefile that is used to create figure 5 in the paper. The figure was created by opening
contamination_pfafstetter_basin.shp in Qgis and exporting it as a map with leaflet (plugin qgis2web).
