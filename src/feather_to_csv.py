import feather
import pandas as pd
import os

print(os.getcwd())
wght_df = feather.read_dataframe('nhm_metdata_weights.feather')
print(wght_df.head())

wght_df.to_csv('hru_metdata_weights.csv')

wght_df2 = pd.read_csv('hru_metdata_weights.csv')
print(wght_df2.head())