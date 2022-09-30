import pandas
import os

directory = os.path.join(os.path.dirname(os.getcwd()), 'data')

# input and output
data_frame_location = os.path.join(directory, "AGG_WWTP_df_no_treatment.csv")
output_location = os.path.join(directory, "AGG_WWTP_df.csv")

AGG_WWTP_df = pandas.read_csv(data_frame_location)

AGG_WWTP_df['Treatment_level'] = 0  # first treatment level; no treatment.
AGG_WWTP_df['Treatment_level'] += AGG_WWTP_df['Primary']  # 1 is primary
AGG_WWTP_df['Treatment_level'][AGG_WWTP_df['Secondary'] == 1] = 2  # 2 secondary

# define secondary + and secondary ++
AGG_WWTP_df['Biological treatment'] = AGG_WWTP_df['NRemoval'] + AGG_WWTP_df['PRemoval']

#AGG_WWTP_df['Treatment_level'][AGG_WWTP_df['Biological treatment'] == 1] = 3  # 3 is secondary plus
#AGG_WWTP_df['Treatment_level'][AGG_WWTP_df['Biological treatment'] == 2] = 4  # 4 is secondary plus plus
AGG_WWTP_df['Treatment_level'][AGG_WWTP_df['Ozonation'] == 1] = 3  # 5 is tertiary

fields = ['pixel_number', 'Treatment_level', 'Treat_a', 'Filt_a', 'Unfilt_a', 'pollution', 'dcpLatitud', 'dcpLongitu',
          'countryID']
AGG_WWTP_df = AGG_WWTP_df[fields]  # extracts only specified fields
AGG_WWTP_df.to_csv(output_location)  # saves output.
