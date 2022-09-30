import pandas
import numpy
import os

directory = os.path.join(os.path.dirname(os.getcwd()), 'data')

# input
observation_df_location = os.path.join(directory, "Raw data/Wilkinson_z_normalized.csv")
moved_obs_location = os.path.join(directory, "Raw data/adjust observations/adjust coord.csv")
removed_obs_location = os.path.join(directory, "Raw data/adjust observations/remove_coord.csv")

observation_df = pandas.read_csv(observation_df_location)

# output
pollution_observed_adapted = os.path.join(directory, "Wilkinson_z_normalized_adapted.csv")

# location changes:
moved_obs = pandas.read_csv(moved_obs_location)
original_coord = moved_obs.iloc[:, 0]
moved_coord = moved_obs.iloc[:, 1]

for i in range(len(original_coord)):
    location = list(observation_df['Coord'] == original_coord[i])
    observation_df['Coord'][location] = moved_coord[i]

# deletions
removed_obs_df = pandas.read_csv(removed_obs_location)
check_values = removed_obs_df.iloc[:, 0]
for i in range(len(removed_obs_df)):
    location = list(observation_df['Coord'] == check_values[i])
    observation_df = observation_df.drop(observation_df.index[location])

observation_df.to_csv(pollution_observed_adapted)



