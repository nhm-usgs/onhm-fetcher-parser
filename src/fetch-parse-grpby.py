
import geopandas as gpd
import pandas as pd
import netCDF4
import numpy as np
import glob
import zipfile
import rasterio
import os
import xarray as xr
import json
from rasterstats import zonal_stats
from rasterio.transform import from_origin

from nc_to_geotiff import *

def testnan(value, weight):
    if np.isnan(value):
        tvalue = 0.0
        tweight = 0.0
    else:
        tvalue = value
        tweight = weight
    return tvalue, tweight

def np_get_wval(grp, ndata):
    mdata = np.ma.masked_array(ndata[grp['grid_ids'].values.astype(int)],
                               np.isnan(ndata[grp['grid_ids'].values.astype(int)]))
    return np.ma.average(mdata, weights=grp['w'])

# Some constants
ktoc=273.15

#=========================================================
# Read all hru shapefiles into single geopandas dataframe
#=========================================================

print(os.getcwd())
from pathlib import Path
folder = Path(r'../Data') # assumes working directory is onhm-fetcher-parser
print(folder)
shapefiles = folder.glob("*.shp")
gdf = pd.concat([
    gpd.read_file(shp)
    for shp in shapefiles
]).pipe(gpd.GeoDataFrame)
gdf.reset_index(drop=True, inplace=True)

num_hru = len(gdf.index)
# print(gdf.head())
print('finished reading shapefiles')
# add new columns for temp and precip
gdf['tmax']=None
gdf['tmin']=None
gdf['ppt']=None


#=========================================================
#            GET CLIMATE DATA
#=========================================================
dirPath='http://thredds.northwestknowledge.net:8080'
tmaxfile='/thredds/dodsC/MET/tmmx/tmmx_2019.nc'
tminfile='/thredds/dodsC/MET/tmmn/tmmn_2019.nc'
pptfile='/thredds/dodsC/MET/pr/pr_2019.nc'

dstmax = xr.open_dataset(dirPath+tmaxfile)
dstmin = xr.open_dataset(dirPath+tminfile)
dstppt = xr.open_dataset(dirPath+pptfile)

print('finished opening netcdf files')
#=========================================================
# Get Shape/Lat/Lon/DataArray
# All the datahandle should be the same for each so will
# try to open one lat, lon, day, crs from tmax to start
#=========================================================

lat_h=dstmax['lat']
lon_h=dstmax['lon']
time_h=dstmax['day']
crs_h=dstmax['crs']

tmax_h=dstmax['air_temperature']
tmin_h=dstmin['air_temperature']
tppt_h=dstppt['precipitation_amount']

ts = tmax_h.sizes
dayshape = ts['day']
lonshape = ts['lon']
latshape = ts['lat']


#=========================================================
#       Read hru weights
#=========================================================

# wght_df_40 = pd.read_csv(r'./Data/hru_metdata_weights_40m.csv')
wght_UofI = pd.read_csv('../Data/hru_uofimetdata_weights.csv')
print('finished reading weight file')
unique_hru_ids = wght_UofI.groupby('hru_id_nat')
#=========================================================
#       Iterate over hrus and assign tmax, tmin, ppt
#       based on weight file
#=========================================================
tind = 0
zcount = 0

#=========================================================
#   flatten data arrays
#=========================================================
tmax_h_flt = np.nan_to_num(tmax_h.values[dayshape - 1,:,:].flatten(order='K'))
tmin_h_flt = np.nan_to_num(tmax_h.values[dayshape - 1,:,:].flatten(order='K'))
tppt_h_flt = np.nan_to_num(tmax_h.values[dayshape - 1,:,:].flatten(order='K'))

np_tmax = np.zeros(num_hru)
np_tmin = np.zeros(num_hru)
np_ppt = np.zeros(num_hru)

for index, row in gdf.iterrows():

    # weight_id_rows = wght_df_40.loc[wght_df_40['hru_id_nat'] == row['hru_id_nat']]
    weight_id_rows = unique_hru_ids.get_group(row['hru_id_nat'])
    np_tmax[index] = np_get_wval(weight_id_rows, tmax_h_flt) - 273.5
    np_tmin[index] = np_get_wval(weight_id_rows, tmin_h_flt) - 273.5
    np_ppt[index] = np_get_wval(weight_id_rows, tppt_h_flt)
#     ttmax = ttmin = tppt = 0.0
#     tmaxwght = tminwght = pptwght =0.0
#     tcount = 0
#
#     # based on metadata of the netcdf file the shape of the netcdf array is day,lat(y),lon(x)
#     for ind2, rw2 in weight_id_rows.iterrows():
# #           print(rw2['Y_ind'],rw2['X_ind'])
#         tmaxv, tmaxwt = testnan(tmax_h.values[dayshape - 1, int(rw2['Y_ind'])-1, int(rw2['X_ind'])-1], rw2['w'])
#         if tmaxwt > 0.0:
#             ttmax += tmaxwt * tmaxv
#             tmaxwght += tmaxwt
#             tcount += 1
#         tminv, tminwt = testnan(tmin_h.values[dayshape - 1, int(rw2['Y_ind'])-1, int(rw2['X_ind'])-1], rw2['w'])
#         if tminwt > 0.0:
#             ttmin += tminwt * tminv
#             tminwght += tminwt
#
#         pptv, pptwt = testnan(tppt_h.values[dayshape - 1, int(rw2['Y_ind']-1), int(rw2['X_ind'])-1], rw2['w'])
#         if pptwt > 0.0:
#             tppt += pptwt * pptv
#             pptwght += pptwt
    if index % 10000 == 0:
        print(index, row['hru_id_nat'])
#     if len(weight_id_rows) > 0 and (tmaxwght or tminwght or pptwght > 0.0):
#         if tmaxwght <= 0.0 or tminwght <= 0.0 or pptwght <= 0.0:
#             print(tind, row['hru_id_nat'], rw2['w'],
#                   ttmax, tmaxwght, ttmin,
#                   tminwght, tppt, pptwght)
#         gdf.loc[gdf.index[tind], 'tmax'] = ((ttmax / tmaxwght) - ktoc)
#         gdf.loc[gdf.index[tind], 'tmin'] = ((ttmin / tminwght) - ktoc)
#         gdf.loc[gdf.index[tind], 'ppt'] = tppt / pptwght
#     else:
#         zcount += 1
#         gdf.loc[gdf.index[tind], 'tmax'] = 0
#         gdf.loc[gdf.index[tind], 'tmin'] = 0
#         gdf.loc[gdf.index[tind], 'ppt'] = 0
    tind += 1



# try: ncfile.close() # just to be safe, make sure dataset is not already open.
# except: pass
print('zcount = ', zcount)
ncfile = netCDF4.Dataset('new.nc',mode='w',format='NETCDF4_CLASSIC')

# Global Attributes
ncfile.Conventions = 'CF-1.8'
ncfile.featureType = 'timeSeries'
ncfile.history = ''

sp_dim = len(gdf.index)
hruid_dim = ncfile.createDimension('hruid', sp_dim)     # hru_id
time_dim = ncfile.createDimension('time', None) # unlimited axis (can be appended to).
for dim in ncfile.dimensions.items():
    print(dim)

#Create Variables
time = ncfile.createVariable('time', np.int, ('time', ))
time.long_name = 'time'
time.standard_name = 'time'
time.units = 'days since '+'base_date'+' 00:00'+'time_zone'

hru = ncfile.createVariable('hruid', np.int, ('hruid', ))
hru.cf_role = 'timeseries_id'
hru.long_name = 'local model hru id'

lat = ncfile.createVariable('hru_lat', np.float32, ('hruid',))
lat.long_name = 'Latitude of HRU centroid'
lat.units = 'degrees_north'
lat.standard_name = 'hru_latitude'

lon = ncfile.createVariable('hru_lon', np.float32, ('hruid',))
lon.long_name = 'Longitude of HRU centroid'
lon.units = 'degrees_east'
lon.standard_name = 'hru_longitude'

prcp = ncfile.createVariable('prcp', np.float32, ('time', 'hruid'))
prcp.long_name = 'Daily precipitation rate'
prcp.units = 'mm/day'
prcp.standard_name = 'lwe_precipitation_rate'

tmax = ncfile.createVariable('tmax', np.float32, ('time', 'hruid'))
tmax.long_name = 'Maximum daily air temperature'
tmax.units = 'degree_Celsius'
tmax.standard_name = 'maximum_daily_air_temperature'

tmin = ncfile.createVariable('tmin', np.float32, ('time', 'hruid'))
tmin.long_name = 'Minimum daily air temperature'
tmin.units = 'degree_Celsius'
tmin.standard_name = 'minimum_daily_air_temperature'

# fill variables with available data
def getXY(pt):
    return (pt.x, pt.y)
centroidseries = gdf['geometry'].centroid
tlon, tlat = [list(t) for t  in zip(*map(getXY, centroidseries))]
# print(lon, lat)
lon[:] = tlon
lat[:] = tlat
hru[:] = gdf['hru_id_nat'].values
# print(hruid)
# tmax[0,:] = gdf['tmax'].values
# tmin[0,:] = gdf['tmin'].values
# prcp[0,:] = gdf['ppt'].values

tmax[0,:] = np_tmax
tmin[0,:] = np_tmin
prcp[0,:] = np_ppt

ncfile.close(); print("dataset is closed")