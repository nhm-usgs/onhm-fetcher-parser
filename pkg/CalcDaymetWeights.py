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
from datetime import datetime, timedelta
from urllib.parse import urlencode
from pathlib import Path

def distance(p1x, p1y, p2x, p2y):
    return np.sqrt(np.power((p2x-p1x),2) + np.power((p2y-p1y),2))

# Grab an example .nc file from Daymet Thredds server
#daymet url:https://thredds.daac.ornl.gov/thredds/ncss/ornldaac/1328/2018/daymet_v3_tmax_2018_na.nc4? \
#var=lat&var=lon&var=tmax&north=54&west=-126&east=-65&south=23&disableProjSubset=on&horizStride=1& \
#time_start=2018-01-01T12%3A00%3A00Z&time_end=2018-01-01T12%3A00%3A00Z&timeStride=1&accept=netcdf
# file saved as: daymet_v3_tmax_2018_na.nc4.nc

dm_file = Path(r'../Data/daymet_v3_tmax_2018_na.nc4.nc')
if dm_file.exists():
    print('netcdf file exists opening data in xarray')
    ds = xr.open_dataset(dm_file)
    print('finished reading netcdf file')
else:
    print('netcdf file does not exist - downloading ...')
    prcpurl = 'https://thredds.daac.ornl.gov/thredds/ncss/ornldaac/1328/2018/daymet_v3_tmax_2018_na.nc4'
    prcppayload = {
        'var': 'lat&var=lon&var=tmax',
        'north': '54',
        'west': '-126',
        'east': '-65',
        'south': '23',
        'disableProjSubset': 'on',
        'horizStride': '1',
        'time_start': '2018-01-01T12:00:00Z',
        'time_end': '2018-01-01T12:00:00Z',
        'timeStride': '1',
        'accept': 'netcdf'}
    try:
        s = requests.Session()
        # https://github.com/psf/requests/issues/1454
        qry = urlencode(prcppayload).replace('%26', '&')
        qry = qry.replace('%3D', '=')
        print(qry)
        tmaxfile = requests.get(prcpurl, params=qry)
        tmaxfile.raise_for_status()
    except HTTPError as http_err:
        print(f'HTTP error occured: {http_err}')
    except Exception as err:
        print(f'Other error occured: {err}')
    else:
        print('daymet data retrieved!')
        with open(r'../Data/daymet_v3_tmax_2018_na.nc4.nc', 'wb') as fh:
            fh.write(tmaxfile.content)
        fh.close
        dm_file = Path(r'../Data/daymet_v3_tmax_2018_na.nc4.nc')
        ds = xr.open_dataset(dm_file)
        print('finished downloading netcdf file')

print(ds)
#Read NHM hru shapefiles into geopandas dataframe
print(os.getcwd())

folder = Path(r'../Data') # assumes working directory is onhm-fetcher-parser
print(folder)
# shapefiles = folder.glob("*_0[1].shp")
shapefiles = folder.glob("*.shp")
print('reading hru shapefiles')
gdf = pd.concat([
    gpd.read_file(shp)
    for shp in shapefiles
]).pipe(gpd.GeoDataFrame)
gdf.reset_index(drop=True, inplace=True)
print('finished reading hru shapefiles')
# add new column uofi_tmax dataframe
gdf['tmax'] = 0.0

# create some variables from the UofI netcdf file for use later
# print('\n The meta data is: \n', json.dumps(ds.attrs, indent=4))
lathandle = ds['lat']
lonhandle = ds['lon']
timehandle = ds['time']
datahandle = ds['tmax']

# Print some information on the data

print('\n Data attributes, sizes, and coords \n')
# print('\n Data attributes are: \n', json.dumps(datahandle.attrs, indent=4))
print('\n Data sizes are: \n', datahandle.sizes)
print('\n Data coords are: \n', datahandle.coords)

# ts = datahandle.sizes
# print(type(ts))
# print(ts['day'])
# dayshape = ts['day']
# Lonshape = ts['lon']
# Latshape = ts['lat']
# # dayshape,lonshape,latshape = datahandle.values.shape
# print(dayshape, Lonshape, Latshape)

# temp = ds.tmax
# temp1 = temp.isel(day=0)
# temp1.plot()

lon = ds.lon.values
lat = ds.lat.values
df = pd.DataFrame({'temperature': ds.tmax.values.flatten()})
res = 0.04166666/2.0
numcells = np.shape(lat)[0]*np.shape(lat)[1]
poly = []
index = np.zeros(numcells)
count = 0
# ncfcell = gpd.GeoDataFrame()
# ncfcell['geometry'] = None

df = pd.DataFrame({'temperature': ds.tmax.values.flatten()})
res = 0.04166666 / 2.0
numcells = (np.shape(lat)[0] - 2) * (
            np.shape(lat)[1] - 2)  # -2 to ignore boundaries, daymet domain should well overlap conus
poly = []
index = np.zeros(numcells, dtype=int)
count = 0
# ncfcell = gpd.GeoDataFrame()
# ncfcell['geometry'] = None

cell_file = Path(r'../Data/daymet_cells_t.csv')
if cell_file.exists():
    print('daymet cells exist - reading file - may take a while')
    ncfcells = gpd.read_file(cell_file)
    print('finished reading daymet_cells_t.csv file')
else:
    print('daymet cells file does not exist so creating cells and file')
    for i in range(1, np.shape(lon)[0] - 1):
        if i % 10 == 0: print(i)
        for j in range(1, np.shape(lon)[1] - 1):
            tpoly_1_lon = [lon[i, j], lon[i, j - 1], lon[i + 1, j - 1], lon[i + 1, j]]
            tpoly_1_lat = [lat[i, j], lat[i, j - 1], lat[i + 1, j - 1], lat[i + 1, j]]
            tpoly_1 = Polygon(zip(tpoly_1_lon, tpoly_1_lat))
            p1 = tpoly_1.centroid

            tpoly_2_lon = [lon[i, j], lon[i + 1, j], lon[i + 1, j + 1], lon[i, j + 1]]
            tpoly_2_lat = [lat[i, j], lat[i + 1, j], lat[i + 1, j + 1], lat[i, j + 1]]
            tpoly_2 = Polygon(zip(tpoly_2_lon, tpoly_2_lat))
            p2 = tpoly_2.centroid

            tpoly_3_lon = [lon[i, j], lon[i, j + 1], lon[i - 1, j + 1], lon[i - 1, j]]
            tpoly_3_lat = [lat[i, j], lat[i, j + 1], lat[i - 1, j + 1], lat[i - 1, j]]
            tpoly_3 = Polygon(zip(tpoly_3_lon, tpoly_3_lat))
            p3 = tpoly_3.centroid

            tpoly_4_lon = [lon[i, j], lon[i - 1, j], lon[i - 1, j - 1], lon[i, j - 1]]
            tpoly_4_lat = [lat[i, j], lat[i - 1, j], lat[i - 1, j - 1], lat[i, j - 1]]
            tpoly_4 = Polygon(zip(tpoly_4_lon, tpoly_4_lat))
            p4 = tpoly_4.centroid

            lon_point_list = [p1.x, p2.x, p3.x, p4.x]
            lat_point_list = [p1.y, p2.y, p3.y, p4.y]

            poly.append(Polygon(zip(lon_point_list, lat_point_list)))
            index[count] = count
            count += 1

    ncfcells = gpd.GeoDataFrame(df, index=index, geometry=poly)
    # ncfcells.dropna(subset=['temperature'], inplace=True)
    print('finished creating daymet cells - now writting file - may take a while')
    # ncfcells.to_csv('../Data/daymet_cells_2.csv')
    print('finished writing daymet_cells.csv')

ncfcells.head()

print('Creating Spatial Index - This could take some time')
spatial_index = ncfcells.sindex
print('Finished Spatial Index')
tcount = 0
with open('tmp_weights.csv', 'w', newline='') as f:
    writer = csv.writer(f, delimiter=',', quotechar='"', quoting=csv.QUOTE_MINIMAL)
    for index, row in gdf.iterrows():
        count = 0
        if tcount == 0:
            # writer.writerow(['grid_ids', 'nhm_id', 'hru_id_nat', 'w'])
            writer.writerow(['grid_ids', 'hru_id_nat', 'w'])
        possible_matches_index = list(spatial_index.intersection(row['geometry'].bounds))
        if not(len(possible_matches_index) == 0):
            possible_matches = ncfcells.iloc[possible_matches_index]
            precise_matches = possible_matches[possible_matches.intersects(row['geometry'])]
            if not(len(precise_matches) == 0):
                res_intersection = gpd.overlay(gdf.loc[[index]], precise_matches, how='intersection')
                for nindex, row in res_intersection.iterrows():
                    tmpfloat = np.float(res_intersection.area.iloc[nindex]/gdf.loc[[index], 'geometry'].area)
                    # writer.writerow([np.int(precise_matches.index[count]), np.int(row['nhm_id']), np.int(row['hru_id_nat']), tmpfloat])
                    writer.writerow([np.int(precise_matches.index[count]), np.int(row['hru_id_nat']), tmpfloat])
                    count += 1
                tcount += 1
                if tcount%100 == 0:
                    print(tcount, index)
        else:
            print('no intersection', index, np.int(row['hru_id_nat']))
# f.close()