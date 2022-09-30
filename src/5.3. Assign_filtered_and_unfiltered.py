import geopandas
import pandas
import os
import numpy

directory = os.path.join(os.path.dirname(os.getcwd()), 'data')

# input files
country_basins_location = os.path.join(directory, "Basins_cut.shp")
filtered_unfiltered_location = os.path.join(directory, "Raw data/filtered and unfiltered.csv")

# output
country_basins_with_populations = os.path.join(directory, "basin_with_filt_unfilt.shp")

# Opening the basins
country_basins = geopandas.read_file(country_basins_location)
country_basins['filtered population'] = 0  # output field
country_basins['unfiltered population'] = 0  # output field
filtered_unfiltered = pandas.read_csv(filtered_unfiltered_location)

# Obtain the total population of a country with filtered and unfiltered.
total_populations = country_basins[["CNTR_ID", "population"]].groupby("CNTR_ID").sum()
total_populations = pandas.merge(total_populations, filtered_unfiltered, left_on="CNTR_ID", right_on="Country")
total_populations['Filtered_pop'] = total_populations['population'] * total_populations['Filtered'] / 100
total_populations['Unfiltered_pop'] = total_populations['population'] * total_populations['Unfiltered'] / 100
country_names = set(total_populations['Country'])

# Find a suitable cutoff point between urban and rural per country
standard_rural_urban_distinction = 100
countries_dataframe = []  # not used, but stores some information on cutoff points and percentages

country_names = set(total_populations['Country'])
for country in country_names:  # for each country
    cutoff = standard_rural_urban_distinction

    # prepare some variables for the specific country
    current_basins = country_basins[country_basins['CNTR_ID'] == country]
    populations_country = total_populations[total_populations['Country'] == country]
    current_population = populations_country['population'].item()
    filtered_population = populations_country['Filtered_pop'].item()
    unfiltered_population = populations_country['Unfiltered_pop'].item()

    rural_population = urban_population = -1  # some small number to enter the while loop
    cutoff_list = []  # list of previously used cutoff values to adjust step size.
    loop_count = 0  # counter that may not exceed a certain value to prevent infinite looping
    break_out = False  # set to True in case the code enters (what appears to be) an infinite loop

    # adjust the cutoff as long as we cannot write unfiltered_pop or filtered_pop as a percentage smaller than 100
    # of rural and urban pop.
    while rural_population < filtered_population + unfiltered_population:
        last_cutoff = cutoff  # used to adjust step size

        # calculate urban and rural populations given the cutoff
        rural_basins = current_basins[current_basins["pop_dens"] < cutoff]
        if len(rural_basins) == 0:
            rural_population = 0
        else:
            rural_population = rural_basins[["CNTR_ID", "population"]].groupby("CNTR_ID").sum().iloc[0].item()

        # adjust cutoff if unfiltered or filtered populations are not a fraction of urban and rural populations
        # respectively.
        cutoff += 5

        loop_count += 1
        if loop_count > 200:
            break_out = True
            break
    if break_out:  # if no satisfactory urban-rural distinction is found, simply assign percentages to population
        filtered_rural = filtered_population / current_population
        unfiltered_rural = unfiltered_population / current_population

        # writing to the basins
        basin_indicators = country_basins['CNTR_ID'] == country
        country_basins['filtered population'][basin_indicators] = \
            country_basins['population'][basin_indicators] * filtered_rural
        country_basins['unfiltered population'][basin_indicators] = \
            country_basins['population'][basin_indicators] * unfiltered_rural
    else:  # in case a rural_urban distinction is found.
        filtered_rural = filtered_population / rural_population
        unfiltered_rural = unfiltered_population / rural_population

        # writing to the basins
        basin_indicators = country_basins['CNTR_ID'] == country
        rural_indicators = basin_indicators & (country_basins['pop_dens'] < cutoff)

        country_basins['filtered population'][rural_indicators] = \
            country_basins['population'][rural_indicators] * filtered_rural
        country_basins['unfiltered population'][rural_indicators] = \
            country_basins['population'][rural_indicators] * unfiltered_rural

    # Writing some information on cutoffs for reference.
    df_row = [country, cutoff, filtered_rural, unfiltered_rural, break_out]
    countries_dataframe.append(df_row)

# unique adaptation for Serbia. assign urban populations to unfiltered and rural to filtered.
Serbia_basins = country_basins[country_basins["CNTR_ID"] == "RS"]
rural_urban_cutoff_RS = 103
rural_basins = Serbia_basins[Serbia_basins["pop_dens"] < rural_urban_cutoff_RS]
rural_populations = rural_basins[["CNTR_ID", "population"]].groupby("CNTR_ID").sum()
urban_basins = Serbia_basins[Serbia_basins["pop_dens"] > rural_urban_cutoff_RS]
urban_populations = urban_basins[["CNTR_ID", "population"]].groupby("CNTR_ID").sum()
total_populations_RS = Serbia_basins[["CNTR_ID", "population"]].groupby("CNTR_ID").sum()
total_populations_RS["filtered_pop_tot"] = numpy.array(total_populations["Filtered"]
                                                       [total_populations["Country"] == "RS"]) \
                                           * numpy.array(total_populations_RS["population"])

total_populations_RS["filtered_adap"] = total_populations_RS["filtered_pop_tot"] / rural_populations.iloc[0, 0]
total_populations_RS["unfiltered_pop_tot"] = numpy.array(total_populations["Unfiltered"]
                                                         [total_populations["Country"] == "RS"]) \
                                             * numpy.array(total_populations_RS["population"])
total_populations_RS["unfiltered_adap"] = total_populations_RS["unfiltered_pop_tot"] / urban_populations.iloc[0, 0]

# write the data for Serbia to the output
country_basins["filtered population"] = (country_basins["CNTR_ID"] != "RS") * country_basins["filtered population"]
country_basins["filtered population"] += ((country_basins["CNTR_ID"] == "RS") & (country_basins["pop_dens"]
                                                                                 < rural_urban_cutoff_RS)) * \
                                         country_basins["population"] * \
                                         total_populations_RS["filtered_adap"].iloc[0] / 100
country_basins["unfiltered population"] = (country_basins["CNTR_ID"] != "RS") * country_basins["unfiltered population"]
country_basins["unfiltered population"] += ((country_basins["CNTR_ID"] == "RS") & (country_basins["pop_dens"]
                                                                                   >= rural_urban_cutoff_RS)) * \
                                           country_basins["population"] * \
                                           total_populations_RS["unfiltered_adap"].iloc[0] / 100

# saving
country_basins.dropna(inplace=True)  # drops some empty basins
country_basins.to_file(country_basins_with_populations)
