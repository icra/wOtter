import pandas
import os

# Load data
directory = os.path.join(os.path.dirname(os.getcwd()), 'data')

# input

AGG_WWTP_location = os.path.join(directory, "AGG_WWTP_df.csv")
new_pollutant_location = os.path.join(directory, "Diclofenac consumption.csv")

AGG_WWTP = pandas.read_csv(AGG_WWTP_location)
new_pollutant = pandas.read_csv(new_pollutant_location)
new_pollutant = new_pollutant.rename(columns={'Country':'Country_copy'})
AGG_WWTP_new = pandas.merge(AGG_WWTP, new_pollutant, left_on="Country", right_on="Country_copy")
AGG_WWTP_new = AGG_WWTP_new.drop('Country_copy', 1)
AGG_WWTP_new.to_csv(AGG_WWTP_location)
