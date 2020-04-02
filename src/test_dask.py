
import geopandas as gpd
import pandas as pd
import glob
import zipfile
import rasterio
import os
import xarray as xr
import numpy as np
from numpy.ma import masked
import dask.dataframe as dd
import netCDF4
def map_gridmet(df, wght_id, wghts, data):
    for index, row in df.iterrows():
        try:
            weight_id_rows = wghts.get_group(row[wght_id])
            tmp = np.nan_to_num(np_get_wval(data, weight_id_rows, row[wghts_id]) - 273.5)
            df.tmax.at[index] = tmp
            tmp2 = 0
        except:
            df.tmax.at[index] = netCDF4.default_fillvals['f8']

def np_get_wval(ndata, wghts, hru_id):
    """
    Returns weighted average of ndata with weights = grp
    1) mdata = the subset of values associated with the gridmet id's that are mapped to hru_id.
    2) Some of these values may have nans if the gridmet id is outside of conus so only return values
    that are inside of conus
    3) this means that hru's that are entirely outside of conus will return nans which will ultimately,
    outside of this function get assigned zero's.
    4) the value is assigned the weighted average
    :param ndata: float array of data values
    :param wghts: float array of weights
    :param hru_id hru id number
    :return: numpy weighted averaged - masked to deal with nans associated with
            ndata that is outside of the conus.
    """
    mdata = np.ma.masked_array(ndata[wghts['grid_ids'].values.astype(int)],
                               np.isnan(ndata[wghts['grid_ids'].values.astype(int)]))
#     if np.ma.is_masked(mdata):
#         print('returning masked value', hru_id)

    # mdata = np.ma.masked_where(ndata[wghts['grid_ids'].values.astype(int)] <= 0.0,
    #                            (ndata[wghts['grid_ids'].values.astype(int)]))
    tmp = np.ma.average(mdata, weights=wghts['w'])
    if tmp is masked:
#         print('returning masked value', hru_id)
        return netCDF4.default_fillvals['f8'] #np.nan

    else:
        return tmp

print(os.getcwd())

print(os.getcwd())
from pathlib import Path
# folder = Path(r'../Data') # assumes working directory is onhm-fetcher-parser
folder = Path(r'../Data_v1_1') # assumes working directory is onhm-fetcher-parser
print(folder)
# shapefiles = folder.glob("*_0[1-2].shp")
shapefiles = folder.glob("*2e*.shp")
gdf = pd.concat([
    gpd.read_file(shp)
    for shp in shapefiles
]).pipe(gpd.GeoDataFrame)
gdf.reset_index(drop=True, inplace=True)
# gdf.plot()
print(gdf)

import requests
from requests.exceptions import HTTPError
import json

# delete existing file if it exists
gmfile = Path('../Data_v1_1/test_gm3.nc')
exists = gmfile.exists()
if not exists:
    #     ds.close()
    # os.remove(gmfile)
    # print('removed existing file')

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

    with open(gmfile, 'wb') as fh:
        fh.write(myfile.content)
        fh.close()

ds = xr.open_dataset(gmfile)
print(ds)

# df = ds.to_dataframe()

print('\n The meta data is: \n', json.dumps(ds.attrs, indent=4))
lathandle = ds['lat']
lonhandle = ds['lon']
timehandle = ds['day']
# datahandle=ds['air_temperature'] # for non aggragated download
datahandle = ds['daily_maximum_temperature']  # for aggragated download


# collect data to describe geotransform
lonmin = float(ds.attrs['geospatial_lon_min'])
latmax = float(ds.attrs['geospatial_lat_max'])
lonres = float(ds.attrs['geospatial_lon_resolution'])
latres = float(ds.attrs['geospatial_lon_resolution'])

# Print some information on the data

print('\n Data attributes, sizes, and coords \n')
print('\n Data sizes are: \n', datahandle.sizes)
print('\n Data coords are: \n', datahandle.coords)

ts = datahandle.sizes
print(type(ts))
print(ts['day'])
dayshape = ts['day']
Lonshape = ts['lon']
Latshape = ts['lat']

print(dayshape, Lonshape, Latshape)

dfmap = pd.DataFrame(gdf.filter(['nhru_v11']))
print(type(dfmap))
print(dfmap)
nhm_id = dfmap.nhru_v11.values
nhm_id
wght_UofI = pd.read_csv('../Data_v1_1/tmp_Gridmet_weights_hru_v1_1e.csv')
wghts_id = wght_UofI.columns[1]
unique_hru_ids = wght_UofI.groupby(wghts_id)
ndata = datahandle.values[dayshape-1,:,:].flatten(order='K')
dfmap['tmax'] = 0.0
map_gridmet(dfmap, wghts_id, unique_hru_ids, ndata)

dfmap.tmax.plot()