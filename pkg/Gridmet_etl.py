import argparse
from pathlib import Path
import xarray as xr
import pandas as pd
import geopandas as gpd
import datetime as dt
import numpy as np
from numpy.ma import masked
import netCDF4
import os


def finalize(odir, year, gdf, numdays, start_date, wght_id,
             mprcp, mtmax, mtmin, mrhmax, mrhmin, mws):
    print(os.getcwd(), flush=True)
    os.chdir(odir)
    print(os.getcwd(), flush=True)
    ncfile = netCDF4.Dataset('gm_' + 'climate_' + str(year) + '.nc',
                             mode='w', format='NETCDF4_CLASSIC')

    # Global Attributes
    ncfile.Conventions = 'CF-1.8'
    ncfile.featureType = 'timeSeries'
    ncfile.history = ''

    sp_dim = len(gdf.index)

    hruid_dim = ncfile.createDimension('hruid', sp_dim)  # hru_id
    time_dim = ncfile.createDimension('time', numdays)  # unlimited axis (can be appended to).
    for dim in ncfile.dimensions.items():
        print(dim, flush=True)

    # Create Variables
    time = ncfile.createVariable('time', 'f4', ('time',))
    time.long_name = 'time'
    time.standard_name = 'time'
    time.units = 'days since ' + start_date.strftime("%Y-%m-%d %H:%M:%S")
    time.calendar = 'standard'

    hru = ncfile.createVariable('hruid', 'i', ('hruid',))
    hru.cf_role = 'timeseries_id'
    hru.long_name = 'local model hru id'


    lat = ncfile.createVariable('hru_lat', np.dtype(np.float32).char, ('hruid',))
    lat.long_name = 'Latitude of HRU centroid'
    lat.units = 'degrees_north'
    lat.standard_name = 'hru_latitude'

    lon = ncfile.createVariable('hru_lon', np.dtype(np.float32).char, ('hruid',))
    lon.long_name = 'Longitude of HRU centroid'
    lon.units = 'degrees_east'
    lon.standard_name = 'hru_longitude'

    prcp = ncfile.createVariable('prcp', np.dtype(np.float32).char, ('time', 'hruid'))
    prcp.long_name = 'daily total precipitation'
    prcp.units = 'mm/day'
    prcp.standard_name = 'daily_total_precipitation'
    prcp.fill_value = netCDF4.default_fillvals['f8']

    tmax = ncfile.createVariable('tmax', np.dtype(np.float32).char, ('time', 'hruid'))
    tmax.long_name = 'Maximum daily air temperature'
    tmax.units = 'degree_Celsius'
    tmax.standard_name = 'maximum_daily_air_temperature'
    tmax.fill_value = netCDF4.default_fillvals['f8']

    tmin = ncfile.createVariable('tmin', np.dtype(np.float32).char, ('time', 'hruid'))
    tmin.long_name = 'Minimum daily air temperature'
    tmin.units = 'degree_Celsius'
    tmin.standard_name = 'minimum_daily_air_temperature'
    tmin.fill_value = netCDF4.default_fillvals['f8']

    rhmax = ncfile.createVariable('rhmax', np.dtype(np.float32).char, ('time', 'hruid'))
    rhmax.long_name = 'daylight average incident shortwave radiation'
    rhmax.units = 'W/m2'
    rhmax.standard_name = 'daylight_average_incident_shortwave_radiation'
    rhmax.fill_value = netCDF4.default_fillvals['f8']

    rhmin = ncfile.createVariable('rhmin', np.dtype(np.float32).char, ('time', 'hruid'))
    rhmin.long_name = 'snow water equivalent'
    rhmin.units = 'kg/m2'
    rhmin.standard_name = 'snow_water_equivalent'
    rhmin.fill_value = netCDF4.default_fillvals['f8']

    ws = ncfile.createVariable('ws', np.dtype(np.float32).char, ('time', 'hruid'))
    ws.long_name = 'daily average vapor pressure'
    ws.units = 'Pa'
    ws.standard_name = 'daily_average_vapor_pressure'
    ws.fill_value = netCDF4.default_fillvals['f8']

    # fill variables with available data
    def getxy(pt):
        return pt.x, pt.y

    cs = gdf.geometry.apply(lambda x: x.centroid)
    # tlon, tlat = [list(t) for t in zip(*map(getxy, cs))]
    # print(lon, lat)
    time[:] = np.arange(0, numdays, dtype=np.float)
    lon[:] = cs.x.values
    lat[:] = cs.y.values
    hru[:] = np.asarray(gdf.index)
    # print(hruid, flush=True)
    # tmax[0,:] = gdf['tmax'].values
    # tmin[0,:] = gdf['tmin'].values
    # prcp[0,:] = gdf['ppt'].values

    tmax[:, :] = mtmax[:, :]
    tmin[:, :] = mtmin[:, :]
    prcp[:, :] = mprcp[:, :]
    rhmax[:, :] = mrhmax[:, :]
    rhmin[:, :] = mrhmin[:, :]
    ws[:, :] = mws[:, :]

    ncfile.close()
    print("dataset is closed", flush=True)

def get_gpd_from_shapefile(idir):
    shapefiles = idir.glob('*.shp')
    gdf = pd.concat([
        gpd.read_file(shp)
        for shp in shapefiles
    ]).pipe(gpd.GeoDataFrame)
    gdf.reset_index(drop=True, inplace=True)
    return gdf

def get_dm_files(idir, year):
    fprcp = None
    ftmax = None
    ftmin = None

    wsrcp = 'ppt'
    vtmax = 'tmax'
    vtmin = 'tmin'
    vrhmin = 'rhmin'
    vrhmax = 'rhmax'
    vws = 'ws'

    str = f'{year}_gm_{wsrcp}*.nc'
    fprcp = list(idir.glob(f'{year}_gm_{wsrcp}*.nc'))
    ftmax = list(idir.glob(f'{year}_gm_{vtmax}*.nc'))
    ftmin = list(idir.glob(f'{year}_gm_{vtmin}*.nc'))
    frhmax = list(idir.glob(f'{year}_gm_{vrhmax}*.nc'))
    frhmin = list(idir.glob(f'{year}_gm_{vrhmin}*.nc'))
    fws = list(idir.glob(f'{year}_gm_{vws}*.nc'))


    return xr.open_dataset(fprcp[0]), \
           xr.open_dataset(ftmax[0]), \
           xr.open_dataset(ftmin[0]), \
           xr.open_dataset(frhmax[0]), \
           xr.open_dataset(frhmin[0]), \
           xr.open_dataset(fws[0])

def main():
    idir = None
    odir = None
    dmyear = None
    wght_file = None
    date = None
    ndata = None
    # default_double_value = 9.96920996838687e+36
    my_parser = argparse.ArgumentParser(prog='Gridmet_etl',
                                        description='Convert gridmet to prms cbh in netcdf format')

    my_parser.add_argument('-i', '--inputdir', type=str,
                           help='Input directory containing daymet files',
                           default=None, required=True)
    my_parser.add_argument('-o', '--outputdir', type=str,
                           help='Directory to output netcdf prms cbh files',
                           default=None, required=True)
    my_parser.add_argument('-y', '--year', type=int,
                           help='Year of daymet file to permorm etl',
                           default=None, required=True)
    my_parser.add_argument('-w', '--weightsfile', type=str,
                           help='path/weight.csv - path to weight file', metavar='weight_file',
                           default=None, required=True)

    args = my_parser.parse_args()

    if args.inputdir is not None:
        idir = Path(args.inputdir)
        if not idir.exists():
            print(f'Input directory: {idir} - does not exist - exiting', flush=True)
            return
    if args.outputdir is not None:
        odir = Path(args.outputdir)
        if not odir.exists():
            print(f'Output directory: {odir} - does not exist - exiting', flush=True)
            return
    if args.year is not None:
        dmyear = args.year
    if args.weightsfile is not None:
        wght_file = Path(args.weightsfile)
        if not wght_file.exists():
            print(f'Weights file: {wght_file} - does not exist -exiting', flush=True)
    date = dt.datetime(year=int(dmyear), month=1, day=1, hour=0)
    start_date = date
    dprcp, dtmax, dtmin, drhmax, drhmin, dws = get_dm_files(idir, dmyear)

    wght_dm = pd.read_csv(wght_file)
    wght_id = wght_dm.columns[1]
    unique_hru_ids = wght_dm.groupby(wght_id)
    print(f'Using weight id: {wght_id}', flush=True)

    gdf = get_gpd_from_shapefile(idir)
    gdf1 = gdf.sort_values(wght_id).dissolve(by=wght_id)

    numdays = dprcp.sizes['day']
    # numdays = 1
    mprcp = np.zeros((numdays, len(gdf1.index)))
    mtmax = np.zeros((numdays, len(gdf1.index)))
    mtmin = np.zeros((numdays, len(gdf1.index)))
    mrhmax = np.zeros((numdays, len(gdf1.index)))
    mrhmin= np.zeros((numdays, len(gdf1.index)))
    mws = np.zeros((numdays, len(gdf1.index)))

    lon = dprcp.lon.values
    lat = dprcp.lat.values

    def getaverage(data, wghts):
        try:
            v_ave = np.average(data, weights=wghts)
        except ZeroDivisionError:
            v_ave = netCDF4.default_fillvals['f8']
        return v_ave

    tindex = np.asarray(gdf1.index)

    for day in np.arange(numdays):
        print(date, flush=True)
        # if day > 0: break
        d_year = np.zeros((7, len(tindex)))
        ndata = np.zeros((7, (np.shape(lon)[0]) * (np.shape(lat)[0])))
        ndata[0, :] = dprcp.precipitation_amount.values[day,:,:].flatten()
        ndata[1, :] = dtmax.daily_maximum_temperature.values[day,:,:].flatten()
        ndata[2, :] = dtmin.daily_minimum_temperature.values[day,:,:].flatten()
        ndata[3, :] = drhmax.daily_maximum_relative_humidity.values[day,:,:].flatten()
        ndata[4, :] = drhmin.daily_minimum_relative_humidity.values[day,:,:].flatten()
        ndata[5, :] = dws.daily_mean_wind_speed.values[day,:,:].flatten()

        # for index, row in gdf.iterrows():
        for i in np.arange(len(tindex)):
            try:
                weight_id_rows = unique_hru_ids.get_group(tindex[i])
                tw = weight_id_rows.w.values
                tgid = weight_id_rows.grid_ids.values
                d_year[0, i] = getaverage(ndata[0, tgid], tw)
                d_year[1, i] = getaverage(ndata[1, tgid]-273.5, tw)
                d_year[2, i] = getaverage(ndata[2, tgid]-273.5, tw)
                d_year[3, i] = getaverage(ndata[3, tgid], tw)
                d_year[4, i] = getaverage(ndata[4, tgid], tw)
                d_year[5, i] = getaverage(ndata[5, tgid], tw)
            except KeyError:
                d_year[0, i] = netCDF4.default_fillvals['f8']
                d_year[1, i] = netCDF4.default_fillvals['f8']
                d_year[2, i] = netCDF4.default_fillvals['f8']
                d_year[3, i] = netCDF4.default_fillvals['f8']
                d_year[4, i] = netCDF4.default_fillvals['f8']
                d_year[5, i] = netCDF4.default_fillvals['f8']

        date += dt.timedelta(days=1)
        # if i == 0:
        #     break

        # np.nan_to_num traps values that had been partially filled by average calc above
        mprcp[day, :] = np.nan_to_num(d_year[0, :], nan=netCDF4.default_fillvals['f8'])
        mtmax[day, :] = np.nan_to_num(d_year[1, :], nan=netCDF4.default_fillvals['f8'])
        mtmin[day, :] = np.nan_to_num(d_year[2, :], nan=netCDF4.default_fillvals['f8'])
        mrhmax[day, :] = np.nan_to_num(d_year[3, :], nan=netCDF4.default_fillvals['f8'])
        mrhmin[day, :] = np.nan_to_num(d_year[4, :], nan=netCDF4.default_fillvals['f8'])
        mws[day, :] = np.nan_to_num(d_year[5, :], nan=netCDF4.default_fillvals['f8'])

    finalize(odir, dmyear, gdf1, numdays, start_date, wght_id, mprcp, mtmax, mtmin, mrhmax, mrhmin, mws)

if __name__ == "__main__":
    main()