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
    return np.sqrt(np.power((p2x - p1x), 2) + np.power((p2y - p1y), 2))


# Grab an example .nc file from Daymet Thredds server
# daymet url:https://thredds.daac.ornl.gov/thredds/ncss/ornldaac/1328/2018/daymet_v3_tmax_2018_na.nc4? \
# var=lat&var=lon&var=tmax&north=54&west=-126&east=-65&south=23&disableProjSubset=on&horizStride=1& \
# time_start=2018-01-01T12%3A00%3A00Z&time_end=2018-01-01T12%3A00%3A00Z&timeStride=1&accept=netcdf
# file saved as: daymet_v3_tmax_2018_na.nc4.nc
input_f = Path('../Data_v1_1')
dm_file = Path('../Data_v1_1/daymet_for_wght.nc')
hruid = 'GFv11_id'
if dm_file.exists():
    print('netcdf file exists opening data in xarray', flush=True)
    # ds = xr.open_dataset(dm_file)
    print('finished reading netcdf file', flush=True)

else:
    print('netcdf file does not exist - downloading ...', flush=True)
    prcpurl = 'https://thredds.daac.ornl.gov/thredds/ncss/daymet-v3-agg/na.ncml'
    prcppayload = {
        #     'var': 'lat&var=lon&var=tmax',
        'var': 'lat&var=lon&var=prcp&var=srad&var=swe&var=tmax&var=tmin&var=vp',
        'north': '54',
        'west': '-126',
        'east': '-65',
        'south': '20',
        'disableProjSubset': 'on',
        'horizStride': '1',
        'time_start': '2018-12-31T12:00:00Z',
        'time_end': '2018-12-31T12:00:00Z',
        'timeStride': '1',
        'accept': 'netcdf'}
    try:
        s = requests.Session()
        # https://github.com/psf/requests/issues/1454
        qry = urlencode(prcppayload).replace('%26', '&')
        qry = qry.replace('%3D', '=')
        print(qry, flush=True)
        tmaxfile = requests.get(prcpurl, params=qry)
        tmaxfile.raise_for_status()
    except HTTPError as http_err:
        print(f'HTTP error occured: {http_err}', flush=True)
    except Exception as err:
        print(f'Other error occured: {err}', flush=True)
    else:
        print('daymet data retrieved!', flush=True)
        with open(dm_file, 'wb') as fh:
            fh.write(tmaxfile.content)
        fh.close()
        # dm_file = Path(r'../Data/daymet_v3_tmax_2018_na.nc4.nc')
        print('finished downloading netcdf file', flush=True)

ds = xr.open_dataset(dm_file)
print(ds, flush=True)
# Read NHM hru shapefiles into geopandas dataframe
print(os.getcwd(), flush=True)

# folder = Path(r'../Data') # assumes working directory is onhm-fetcher-parser
print(input_f, flush=True)
# shapefiles = folder.glob("*_0[1].shp")
shapefiles = input_f.glob("*.shp")
print('reading hru shapefiles', flush=True)

gdf = pd.concat([
    gpd.read_file(shp)
    for shp in shapefiles
]).pipe(gpd.GeoDataFrame)
gdf.reset_index(drop=True, inplace=True)
print('finished reading hru shapefiles', flush=True)

# add new column uofi_tmax dataframe
gdf['tmax'] = 0.0

print('\n The meta data is: \n', ds.attrs, flush=True)
lathandle = ds['lat']
lonhandle = ds['lon']
timehandle = ds['time']
datahandle = ds['tmax']
dhlat = ds['lat']
dhlon = ds['lon']
crshandle = ds['lambert_conformal_conic']
print('\n The crs meta data is \n', crshandle.attrs, flush=True)
print(datahandle, flush=True)

# collect data to describe geotransform
lonmin = float(ds.attrs['geospatial_lon_min'])
latmax = float(ds.attrs['geospatial_lat_max'])

# Print some information on the data

print('\n Data attributes, sizes, and coords \n', flush=True)
print('\n Data attributes are: \n', datahandle.attrs, flush=True)
print('\n Data sizes are: \n', datahandle.sizes, flush=True)
print('\n Data coords are: \n', datahandle.coords, flush=True)
print('\n Lat coords are: \n', dhlat.attrs, flush=True)

ts = datahandle.sizes
print(type(ts), flush=True)
print(ts['time'], flush=True)
dayshape = ts['time']
Lonshape = ts['x']
Latshape = ts['y']

print(dayshape, Lonshape, Latshape, flush=True)

lon = ds.lon.values
lat = ds.lat.values

# first create dataframe with temp
df = pd.DataFrame({'temperature': ds.tmax.values.flatten()})
res = 0.04166666 / 2.0
numcells2 = (np.shape(lat)[0] - 2) * (
        np.shape(lat)[1] - 2)  # -2 to ignore boundaries, daymet domain should well overlap conus
poly2 = []
index2 = np.zeros(numcells2)
count = 0
# ncfcell = gpd.GeoDataFrame()
# ncfcell['geometry'] = None

for i in range(1, np.shape(lon)[0] - 1):
    if i % 100 == 0: print(i, flush=True)
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

        poly2.append(Polygon(zip(lon_point_list, lat_point_list)))
        index2[count] = count
        count += 1
# ncfcells = gpd.GeoDataFrame(df, index=index, crs=ds['crs'], geometry=poly)
ncfcells2 = gpd.GeoDataFrame(df, index=index2, geometry=poly2)
ncfcells2.dropna(subset=['temperature'], inplace=True)
ncfcells2.head()
print('Creating Spatial Index - This could take some time', flush=True)
spatial_index2 = ncfcells2.sindex
print('Finished Spatial Index', flush=True)

tcount = 0
with open('tmp_daymet_weights_hru_v1_1_rnan.csv', 'w', newline='') as f:
    writer = csv.writer(f, delimiter=',', quotechar='"', quoting=csv.QUOTE_MINIMAL)
    for index, row in gdf.iterrows():
        count = 0
        if tcount == 0:
            writer.writerow(['grid_ids', hruid, 'w'])
        possible_matches_index = list(spatial_index2.intersection(row['geometry'].bounds))
        if not (len(possible_matches_index) == 0):
            possible_matches = ncfcells2.iloc[possible_matches_index]
            precise_matches = possible_matches[possible_matches.intersects(row['geometry'])]
            if not (len(precise_matches) == 0):
                res_intersection = gpd.overlay(gdf.loc[[index]], precise_matches, how='intersection')
                for nindex, row in res_intersection.iterrows():
                    ttt = gdf.loc[[index], 'geometry'].area
                    tmpfloat = np.float(res_intersection.area.iloc[nindex] / ttt[index])
                    writer.writerow([np.int(precise_matches.index[count]), np.int(row[hruid]), tmpfloat])
                    count += 1
                tcount += 1
                if tcount % 100 == 0:
                    print(tcount, index, flush=True)

        else:
            print('no intersection: ', index, np.int(row[hruid]), flush=True)

# f.close()
