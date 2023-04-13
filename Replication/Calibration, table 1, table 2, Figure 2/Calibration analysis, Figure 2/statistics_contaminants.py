import numpy
import pandas
import matplotlib.pyplot as plt
import os
import geopandas

from config.config import DATA_DIR


#set current directory as working directory
os.chdir(os.path.dirname(os.path.abspath(__file__)))


directory = os.path.join(DATA_DIR, 'Lumped 14')

dataframes = ['Lumped 14.shp']  # dataframes are in shapefile format

for dataframe in dataframes:
    df_location = os.path.join(directory, dataframe)

    current_dataframe = geopandas.read_file(df_location)
    plot_name_base = dataframe.split('.')
    plot_name_base = plot_name_base[0]
    ## make it a regression
    current_dataframe['ones'] = 1
    # first plot, all observations
    fig, ax = plt.subplots()
    ax.plot(numpy.log(current_dataframe['Prediction']+1), numpy.log(current_dataframe['Observatio']+1), 'o', markersize=4)
    ax.plot(numpy.log(current_dataframe['Observatio']+1), numpy.log(current_dataframe['Observatio']+1))
    plt.xlabel("Log Estimated Concentrations ng l$^{-1}$ ng$^{-1}$", fontdict={'size': 17, 'family': 'Times New Roman'})
    plt.ylabel("Log Observed Concentrations ng l$^{-1}$ ng$^{-1}$", fontdict={'size': 17, 'family': 'Times New Roman'})
    plt.title("All Occurences", fontdict={'size': 17, 'family': 'Times New Roman'})
    #SSR_model = numpy.sum(numpy.square(current_dataframe['Observatio'] - current_dataframe['Prediction']))
    SSR_model = numpy.sum(numpy.square(current_dataframe['Observatio'] - current_dataframe['Prediction']))
    SSR_mean = numpy.sum(numpy.square(numpy.mean(current_dataframe['Observatio']) - current_dataframe['Observatio']))
    R_2 = 1 - SSR_model/SSR_mean
    R_2 = str(R_2)
    R_2 = R_2[0:4]
    textstr = 'R² =' + R_2
    # these are matplotlib.patch.Patch properties
    props = dict(boxstyle='round', facecolor='wheat', alpha=0.5)

    # place a text box in upper left in axes coords
    ax.text(0.05, 0.95, textstr, transform=ax.transAxes, fontdict={'size': 17, 'family': 'Times New Roman'},
            verticalalignment='top', bbox=props)
    plt.tight_layout()
    plt.savefig(plot_name_base)
    plt.clf()

    # now the plots for the lower half of the observations and the upper half
    median = numpy.median(current_dataframe['discharge'])
    print("median is: " + str(median))

    upper_disch = current_dataframe[current_dataframe['discharge'] > median]
    fig, ax = plt.subplots()
    ax.plot(numpy.log(upper_disch['Prediction']+1), numpy.log(upper_disch['Observatio']+1), 'o', markersize=4)
    ax.plot(numpy.log(upper_disch['Prediction']+1), numpy.log(upper_disch['Prediction']+1))
    plt.xlabel("Log Estimated Concentrations ng l$^{-1}$ ng$^{-1}$", fontdict={'size': 17, 'family': 'Times New Roman'})
    plt.ylabel("Log Observed Concentrations ng l$^{-1}$ ng$^{-1}$", fontdict={'size': 17, 'family': 'Times New Roman'})
    plt.title("Above Median Discharge", fontdict={'size': 17, 'family': 'Times New Roman'})
    SSR_model = numpy.sum(numpy.square(upper_disch['Observatio'] - upper_disch['Prediction']))
    SSR_mean = numpy.sum(numpy.square(numpy.mean(current_dataframe['Observatio']) - upper_disch['Observatio']))
    R_2 = 1 - SSR_model/SSR_mean
    R_2 = str(R_2)
    R_2 = R_2[0:4]
    textstr = 'R² =' + R_2
    # these are matplotlib.patch.Patch properties
    props = dict(boxstyle='round', facecolor='wheat', alpha=0.5)

    # place a text box in upper left in axes coords
    ax.text(0.05, 0.95, textstr, transform=ax.transAxes, fontdict={'size': 17, 'family': 'Times New Roman'},
            verticalalignment='top', bbox=props)
    plt.tight_layout()
    plt.savefig(plot_name_base + '_high_disch')
    plt.clf()

    lower_disch = current_dataframe[current_dataframe['discharge'] <= median]
    fig, ax = plt.subplots()
    ax.plot(numpy.log(lower_disch['Prediction']+1), numpy.log(lower_disch['Observatio']+1), 'o', markersize=4)
    ax.plot(numpy.log(lower_disch['Prediction']+1), numpy.log(lower_disch['Prediction']+1))
    plt.xlabel("Log Estimated Concentrations ng l$^{-1}$ ng$^{-1}$", fontdict={'size': 17, 'family': 'Times New Roman'})
    plt.ylabel("Log Observed Concentrations ng l$^{-1}$ ng$^{-1}$", fontdict={'size': 17, 'family': 'Times New Roman'})
    plt.title("Below Median Discharge", fontdict={'size': 17, 'family': 'Times New Roman'})
    SSR_model = numpy.sum(numpy.square(lower_disch['Observatio'] - lower_disch['Prediction']))
    SSR_mean = numpy.sum(numpy.square(numpy.mean(current_dataframe['Observatio']) - lower_disch['Observatio']))
    R_2 = 1 - SSR_model/SSR_mean
    R_2 = str(R_2)
    R_2 = R_2[0:4]
    textstr = 'R² =' + R_2
    # these are matplotlib.patch.Patch properties
    props = dict(boxstyle='round', facecolor='wheat', alpha=0.5)

    # place a text box in upper left in axes coords
    ax.text(0.05, 0.95, textstr, transform=ax.transAxes, fontdict={'size': 17, 'family': 'Times New Roman'},
            verticalalignment='top', bbox=props)
    plt.tight_layout()
    plt.savefig(plot_name_base + '_low_disch')
    plt.clf()