import geopandas
import pandas
import os


directory = os.path.join(os.getcwd(), 'data')

# changes the discharge of some rivers as specified by the hydrorivers flow excel file.
hydrorivers_location = os.path.join(directory, "Raw data/HydroSHEDS/HydroRIVERS_v10_eu.shp")
corrections_location = os.path.join(directory, "Raw data/hydrorivers flow.csv")
output_name = os.path.join(directory, "hydro_rivers_adapted.shp")

hydro_rivers = geopandas.GeoDataFrame.from_file(hydrorivers_location)  # open hydrorivers as a dataframe
corrections_df = pandas.read_csv(corrections_location)
hydro_ids_all = hydro_rivers["HYRIV_ID"]
hydro_ids = corrections_df.iloc[:, 0]
hydro_discharges = corrections_df.iloc[:, 1]
for i in range(len(corrections_df)):
    # apply the changes to the discharge
    hydro_rivers["DIS_AV_CMS"][hydro_ids_all == hydro_ids[i]] = hydro_discharges[i]

hydro_rivers.to_file(output_name)
