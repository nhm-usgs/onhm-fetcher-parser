import geopandas as gpd
import pandas as pd
import numpy as np
import glob
import zipfile
import rasterio
import os
import xarray as xr
import json
from shapely.geometry import Polygon
import csv
import requests
from requests.exceptions import HTTPError

# Grab an example .nc file from Gridmet Thredds server
#university of Idaho url:http://thredds.northwestknowledge.net:8080/thredds/ncss/agg_met_tmmx_1979_CurrentYear_CONUS.nc?
# var=daily_maximum_temperature&disableLLSubset=on&disableProjSubset=on&horizStride=1&time_start=2014-01-01T00%3A00%3A00Z
# &time_end=2014-01-01T00%3A00%3A00Z&timeStride=1&accept=netcdf
# file saved as: agg_met_tmmx_1979_CurrentYear_CONUS_2014_01_01.nc

url2 = 'http://thredds.northwestknowledge.net:8080/thredds/ncss/agg_met_tmmx_1979_CurrentYear_CONUS.nc'
payload2 = {'var': 'daily_maximum_temperature',
            'disableLLSubset': 'on',
            'disableProjSubset': 'on',
            'horizStride': '1',
            'time_start': '2018-12-31T00:00:00Z',
            'time_end': '2018-12-31T00:00:00Z',
            'timeStride': '1',
            'accept': 'netcdf'}
# print(url)
try:
    #     myfile = requests.get(url, params=payload)
    myfile = requests.get(url2, params=payload2)
    myfile.raise_for_status()
except HTTPError as http_err:
    print(f'HTTP error occurred: {http_err}')  # Python 3.6
except Exception as err:
    print(f'Other error occurred: {err}')  # Python 3.6
else:
    print('Success!')
    #     print(myfile.headers)
    print(myfile.url)

with open(r'../DeschutesData/test3.nc', 'wb') as fh:
    fh.write(myfile.content)
    fh.close()

uofi_file = r'../DeschutesData/test3.nc'
ds = xr.open_dataset(uofi_file)
print(ds)
#Read NHM hru shapefiles into geopandas dataframe
print(os.getcwd())
from pathlib import Path
folder = Path(r'../DeschutesData') # assumes working directory is onhm-fetcher-parser
print(folder)
# shapefiles = folder.glob("*_0[1-2].shp")
shapefiles = folder.glob("*.shp")
gdf = pd.concat([
    gpd.read_file(shp)
    for shp in shapefiles
]).pipe(gpd.GeoDataFrame)
gdf.reset_index(drop=True, inplace=True)

# add new column uofi_tmax dataframe
gdf['tmax'] = 0.0

# create some variables from the UofI netcdf file for use later
print('\n The meta data is: \n', json.dumps(ds.attrs, indent=4))
lathandle = ds['lat']
lonhandle = ds['lon']
timehandle = ds['day']
datahandle = ds['daily_maximum_temperature']

# collect data to describe geotransform
lonmin = float(ds.attrs['geospatial_lon_min'])
latmax = float(ds.attrs['geospatial_lat_max'])
lonres = float(ds.attrs['geospatial_lon_resolution'])
latres = float(ds.attrs['geospatial_lon_resolution'])

# Print some information on the data

print('\n Data attributes, sizes, and coords \n')
# print('\n Data attributes are: \n', json.dumps(datahandle.attrs, indent=4))
print('\n Data sizes are: \n', datahandle.sizes)
print('\n Data coords are: \n', datahandle.coords)

ts = datahandle.sizes
print(type(ts))
print(ts['day'])
dayshape = ts['day']
Lonshape = ts['lon']
Latshape = ts['lat']
# dayshape,lonshape,latshape = datahandle.values.shape
print(dayshape, Lonshape, Latshape)

temp = ds.daily_maximum_temperature
temp1 = temp.isel(day=0)
temp1.plot()

lon, lat = np.meshgrid(lonhandle, lathandle)
df = pd.DataFrame({'temperature': temp1.values.flatten()})
res = 0.04166666/2.0
numcells = np.shape(lat)[0]*np.shape(lat)[1]
poly = []
index = np.zeros(numcells)
count = 0
# ncfcell = gpd.GeoDataFrame()
# ncfcell['geometry'] = None

for i in range(np.shape(lon)[0]):
    for j in range(np.shape(lon)[1]):
        lat_point_list = [lat[i,j]-res, lat[i,j]+res, lat[i,j]+res, lat[i,j]-res]
        lon_point_list = [lon[i,j]+res, lon[i,j]+res, lon[i,j]-res, lon[i,j]-res]
        poly.append(Polygon(zip(lon_point_list, lat_point_list)))
        index[count] = count
        count += 1
# ncfcells = gpd.GeoDataFrame(df, index=index, crs=ds['crs'], geometry=poly)
ncfcells = gpd.GeoDataFrame(df, index=index, geometry=poly)
# ncfcells.dropna(subset=['temperature'], inplace=True)
ncfcells.head()


spatial_index = ncfcells.sindex
tcount = 0
with open('../DeschutesData/Deshutes_weights_hru_v3.csv', 'w', newline='') as f:
    writer = csv.writer(f, delimiter=',', quotechar='"', quoting=csv.QUOTE_MINIMAL)
    for index, row in gdf.iterrows():
        count = 0
        if tcount == 0:
            writer.writerow(['grid_ids', 'HRU_ID', 'w'])
        possible_matches_index = list(spatial_index.intersection(row['geometry'].bounds))
        # if not(len(possible_matches_index) == 0):
        possible_matches = ncfcells.iloc[possible_matches_index]
        precise_matches = possible_matches[possible_matches.intersects(row['geometry'])]
        # if not(len(precise_matches) == 0):
        res_intersection = gpd.overlay(gdf.loc[[index]], precise_matches, how='intersection')
        for nindex, row in res_intersection.iterrows():
            tmpfloat = np.float(res_intersection.area.iloc[nindex]/gdf.loc[[index], 'geometry'].area)
            writer.writerow([np.int(precise_matches.index[count]), np.int(row['HRU_ID']), tmpfloat])
            count += 1
        tcount += 1
        if tcount%10 == 0:
            print(tcount, index)
        # else:
        #     print('no intersection: ', index, np.int(row['nhm_id']))
# f.close()