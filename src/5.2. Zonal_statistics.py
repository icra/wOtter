from osgeo import ogr
from rasterstats import zonal_stats
import fiona
from src.library import area
import os

directory = os.path.join(os.getcwd(), 'data')

# inputs
basin_location = os.path.join(directory, "Basins_cut.shp")  # this will also become the output
pop_raster = os.path.join(directory, "Raw data/Population.tif")
estimator_pop_raster = os.path.join(directory, "Population_treated.tif")

# find the zonal statistics for the basin geometries
basin_shapefile = fiona.open(basin_location)

# calculate the sum of the population and treated population in each polygon of the basin shapefile
stats_pop = zonal_stats(basin_shapefile, pop_raster, stats=["sum"])
stats_pop = [element.get('sum') for element in stats_pop]  # gives a list of populations for the basins

stats_treated = zonal_stats(basin_shapefile, estimator_pop_raster, stats="sum")
stats_treated = [element.get('sum') for element in stats_treated]  # gives a list of treated populations for the basins

basin_shapefile = None

# write the zonal statistics for level 7
basin_shapefile = ogr.Open(basin_location, 1)
basin_layer = basin_shapefile.GetLayer(0)
basin_layer.CreateField(ogr.FieldDefn('population', ogr.OFTInteger64))
basin_layer.CreateField(ogr.FieldDefn('treated_p', ogr.OFTInteger64))
basin_layer.CreateField(ogr.FieldDefn('Area', ogr.OFTReal))
basin_layer.CreateField(ogr.FieldDefn('pop_dens', ogr.OFTReal))

basin_id = 0
for basin in basin_layer:
    basin.SetField('population', stats_pop[basin_id])
    basin.SetField('treated_p', stats_treated[basin_id])
    geom = basin.geometry()
    area_km = 0.000001 * area.area(geom.ExportToJson())
    basin.SetField("Area", area_km)
    try:
        pop_dens = stats_pop[basin_id] / area_km
    except:
        pop_dens = 0
    basin.SetField("Pop_dens", pop_dens)
    basin_layer.SetFeature(basin)
    basin_id += 1
    basin = None
basin_layer = None
basin_shapefile = None
