from osgeo import ogr
from shapely.geometry import shape, mapping
import fiona
import shutil
import os


directory = os.path.join(os.getcwd(), 'data')

# inputs
basins_location = os.path.join(directory, "Raw data/hydroBASINS/Basins.shp")
all_countries = os.path.join(directory, "Raw data/countries original shp/CNTR_RG_01M_2020_4326")
countries_included_location = os.path.join(directory, "Raw data/Countries_included.txt")

# outputs
countries = os.path.join(directory, "countries/Countries")
basins_cut_location = os.path.join(directory, "Basins_cut.shp")

# make a copy
try:
    os.makedirs(os.path.join(directory, "countries"))
except FileExistsError:
    pass

for element in ('.cpg', '.dbf', '.prj', '.shp', '.shx'):
    shutil.copyfile(all_countries + element, countries + element)
all_countries += '.shp'  # to open the file the extension needs to be included
countries += '.shp'

# Remove all the countries that are not considered
country_polygons = ogr.Open(all_countries, 1)
country_polygons_layer = country_polygons.GetLayer()
countries_text = open(countries_included_location)
countries_included = countries_text.read()
countries_included = set(countries_included.split(", "))

for feature in country_polygons_layer:
    if feature.GetField("CNTR_ID") not in countries_included:
        country_polygons_layer.DeleteFeature(feature.GetFID())  # removes the country from the shapefile

country_polygons_layer = None
country_polygons = None

#  open the geometries of the basins and make sure they are valid (otherwise the next step causes an error).
#  An invalid geometry is one that self-intersects, for instance. The input data has a few.
basins = ogr.Open(basins_location, 1)
basins_layer = basins.GetLayer(0)

for feature in basins_layer:
    geometry = feature.geometry()
    geometry = geometry.MakeValid() # this makes the geometry object valid
    feature.SetGeometry(geometry)
    basins_layer.SetFeature(feature)
    feature = None
basins_layer = None
basins = None

# create a new shapefile that contains the intersection between the basins and the countries.
basins = fiona.open(basins_location)
countries_included = fiona.open(all_countries)
# creation of the new shapefile with the intersection
schema = basins.schema
schema['properties']['CNTR_ID'] = 'str:2'
with fiona.open(basins_cut_location, 'w', driver='ESRI Shapefile', schema=schema) as output:
    for country in countries_included:
        for basin in basins:
            if shape(basin['geometry']).intersects(shape(country['geometry'])): # if the country and basin intersect
                # create a new polygon that is the intersection, with basin data
                properties = basin['properties']
                properties['CNTR_ID'] = country['properties'].get('CNTR_ID')
                output.write({'geometry': mapping(shape(basin['geometry']).intersection(shape(country['geometry']))),
                              'properties': properties})
