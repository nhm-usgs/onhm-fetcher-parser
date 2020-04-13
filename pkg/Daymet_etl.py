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
import sys

def finalize(odir, year, gdf, numdays, start_date, wght_id,
             mprcp, mtmax, mtmin, msrad, mswe, mvp, mdayl):
    print(os.getcwd(), flush=True)
    os.chdir(odir)
    print(os.getcwd(), flush=True)
    ncfile = netCDF4.Dataset('dm_' + 'climate_' + str(year) + '.nc',
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

    srad = ncfile.createVariable('srad', np.dtype(np.float32).char, ('time', 'hruid'))
    srad.long_name = 'daylight average incident shortwave radiation'
    srad.units = 'W/m2'
    srad.standard_name = 'daylight_average_incident_shortwave_radiation'
    srad.fill_value = netCDF4.default_fillvals['f8']

    swe = ncfile.createVariable('swe', np.dtype(np.float32).char, ('time', 'hruid'))
    swe.long_name = 'snow water equivalent'
    swe.units = 'kg/m2'
    swe.standard_name = 'snow_water_equivalent'
    swe.fill_value = netCDF4.default_fillvals['f8']

    vp = ncfile.createVariable('vp', np.dtype(np.float32).char, ('time', 'hruid'))
    vp.long_name = 'daily average vapor pressure'
    vp.units = 'Pa'
    vp.standard_name = 'daily_average_vapor_pressure'
    vp.fill_value = netCDF4.default_fillvals['f8']

    dayl = ncfile.createVariable('dayl', np.dtype(np.float32).char, ('time', 'hruid'))
    dayl.long_name = 'daylength'
    dayl.units = 's'
    dayl.standard_name = 'daylength'
    dayl.fill_value = netCDF4.default_fillvals['f8']

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
    srad[:, :] = msrad[:, :]
    swe[:, :] = mswe[:, :]
    vp[:, :] = mvp[:, :]
    dayl[:, :] = mdayl[:, :]

    ncfile.close()
    print("dataset is closed", flush=True)


def get_dm_xy_bounds():
    xmin = -2754.25
    xmax = 3252.75
    ymin = -2011.0
    ymax = 1687.0
    return xmin, xmax, ymin, ymax

def np_get_wval(ndata, tgid, wghts, hru_id):
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
    mdata = np.ma.masked_array(ndata, np.isnan(ndata))

    # mdata = np.ma.masked_where(ndata[wghts['grid_ids'].values.astype(int)] <= 0.0,
    #                            (ndata[wghts['grid_ids'].values.astype(int)]))
    tmp = np.ma.average(mdata, weights=wghts)
    return tmp
    # if tmp is masked:
    #     # print(f'returning masked value: {hru_id}', ndata)
    #     return netCDF4.default_fillvals['f8']

    # else:
    #     return tmp

def get_gpd_from_shapefile(idir):
    shapefiles = idir.glob('*.shp')
    gdf = pd.concat([
        gpd.read_file(shp)
        for shp in shapefiles
    ]).pipe(gpd.GeoDataFrame)
    gdf.reset_index(drop=True, inplace=True)
    return gdf

def get_dm_files(idir, year, xmin, xmax, ymin, ymax):
    fprcp = None
    ftmax = None
    ftmin = None

    vprcp = 'prcp'
    vtmax = 'tmax'
    vtmin = 'tmin'
    vsrad = 'srad'
    vswe = 'swe'
    vvp = 'vp'
    vdayl = 'dayl'

    fprcp = idir / f'daymet_v3_{vprcp}_{year}_na.nc4'
    ftmax = idir / f'daymet_v3_{vtmax}_{year}_na.nc4'
    ftmin = idir / f'daymet_v3_{vtmin}_{year}_na.nc4'
    fsrad = idir / f'daymet_v3_{vsrad}_{year}_na.nc4'
    fswe = idir / f'daymet_v3_{vswe}_{year}_na.nc4'
    fvp = idir / f'daymet_v3_{vvp}_{year}_na.nc4'
    fdayl = idir / f'daymet_v3_{vdayl}_{year}_na.nc4'

    return xr.open_dataset(fprcp).sel(
    y=slice(ymax*1000.0, ymin*1000.0),
    x=slice(xmin*1000.0,xmax*1000.0)), \
           xr.open_dataset(ftmax).sel(
    y=slice(ymax*1000.0, ymin*1000.0),
    x=slice(xmin*1000.0,xmax*1000.0)), \
           xr.open_dataset(ftmin).sel(
    y=slice(ymax*1000.0, ymin*1000.0),
    x=slice(xmin*1000.0,xmax*1000.0)), \
           xr.open_dataset(fsrad).sel(
    y=slice(ymax*1000.0, ymin*1000.0),
    x=slice(xmin*1000.0,xmax*1000.0)), \
           xr.open_dataset(fswe).sel(
    y=slice(ymax*1000.0, ymin*1000.0),
    x=slice(xmin*1000.0,xmax*1000.0)), \
           xr.open_dataset(fvp).sel(
    y=slice(ymax*1000.0, ymin*1000.0),
    x=slice(xmin*1000.0,xmax*1000.0)), \
           xr.open_dataset(fdayl).sel(
    y=slice(ymax*1000.0, ymin*1000.0),
    x=slice(xmin*1000.0,xmax*1000.0))

def main():
    idir = None
    odir = None
    dmyear = None
    wght_file = None
    date = None
    ndata = None
    # default_double_value = 9.96920996838687e+36
    my_parser = argparse.ArgumentParser(prog='Daymet_etl',
                                        description='Convert daymet to prms cbh in netcdf format')

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
    date = dt.datetime(year=int(dmyear), month=1, day=1, hour=12)
    start_date = date
    xmin, xmax, ymin, ymax = get_dm_xy_bounds()
    dprcp, dtmax, dtmin, dsrad, dswe, dvp, ddayl = get_dm_files(idir, dmyear, xmin, xmax, ymin, ymax)

    wght_dm = pd.read_csv(wght_file)
    wght_id = wght_dm.columns[1]
    unique_hru_ids = wght_dm.groupby(wght_id)
    print(f'Using weight id: {wght_id}', flush=True)

    gdf = get_gpd_from_shapefile(idir)
    gdf1 = gdf.sort_values(wght_id).dissolve(by=wght_id)

    numdays = dprcp.sizes['time']
    # numdays = 1
    mprcp = np.zeros((numdays, len(gdf1.index)))
    mtmax = np.zeros((numdays, len(gdf1.index)))
    mtmin = np.zeros((numdays, len(gdf1.index)))
    msrad = np.zeros((numdays, len(gdf1.index)))
    mswe= np.zeros((numdays, len(gdf1.index)))
    mvp = np.zeros((numdays, len(gdf1.index)))
    mdayl = np.zeros((numdays, len(gdf1.index)))

    lon = dprcp.lon.values
    lat = dprcp.lat.values

    def getaverage(data, wghts, tindex):
        try:
            v_ave = np.average(data, weights=wghts)
        except ZeroDivisionError:
            # print(f'ZeroDivisionError index: {tindex}')
            v_ave = netCDF4.default_fillvals['f8']
        return v_ave

    tindex = np.asarray(gdf1.index)

    for day in np.arange(numdays):
        print(date, flush=True)
        # if day > 0: break
        d_year = np.zeros((7, len(tindex)))
        ndata = np.zeros((7, (np.shape(lon)[1] - 2) * (np.shape(lon)[0] - 2)))
        ndata[0, :] = dprcp.prcp.sel(time=date).values[1:np.shape(lon)[0] - 1, 1:np.shape(lon)[1] - 1].flatten()
        ndata[1, :] = dtmax.tmax.sel(time=date).values[1:np.shape(lon)[0] - 1, 1:np.shape(lon)[1] - 1].flatten()
        ndata[2, :] = dtmin.tmin.sel(time=date).values[1:np.shape(lon)[0] - 1, 1:np.shape(lon)[1] - 1].flatten()
        ndata[3, :] = dsrad.srad.sel(time=date).values[1:np.shape(lon)[0] - 1, 1:np.shape(lon)[1] - 1].flatten()
        ndata[4, :] = dswe.swe.sel(time=date).values[1:np.shape(lon)[0] - 1, 1:np.shape(lon)[1] - 1].flatten()
        ndata[5, :] = dvp.vp.sel(time=date).values[1:np.shape(lon)[0] - 1, 1:np.shape(lon)[1] - 1].flatten()
        ndata[6, :] = ddayl.dayl.sel(time=date).values[1:np.shape(lon)[0] - 1, 1:np.shape(lon)[1] - 1].flatten()

        # for index, row in gdf.iterrows():
        fill = np.zeros(7)
        for i in np.arange(len(tindex)):
            try:
                weight_id_rows = unique_hru_ids.get_group(tindex[i])
                tw = weight_id_rows.w.values
                tgid = weight_id_rows.grid_ids.values
                
                if np.isnan(getaverage(ndata[0, tgid], tw, tindex[i])):
                    # np.set_printoptions(threshold=sys.maxsize)
                    # print(len(tw), tw)
                    # print(len(tgid), tgid)
                    # print(len(ndata[1,tgid]), ndata[1, tgid])
                    # # print(getaverage(ndata[1, tgid], tw))

                    d_year[0, i] = np_get_wval(ndata[0, tgid], tgid, tw, tindex[i])
                    d_year[1, i] = np_get_wval(ndata[1, tgid], tgid, tw, tindex[i])
                    d_year[2, i] = np_get_wval(ndata[2, tgid], tgid, tw, tindex[i])
                    d_year[3, i] = np_get_wval(ndata[3, tgid], tgid, tw, tindex[i])
                    d_year[4, i] = np_get_wval(ndata[4, tgid], tgid, tw, tindex[i])
                    d_year[5, i] = np_get_wval(ndata[5, tgid], tgid, tw, tindex[i])
                    d_year[6, i] = np_get_wval(ndata[6, tgid], tgid, tw, tindex[i])
                else:
                    d_year[0, i] = getaverage(ndata[0, tgid], tw, tindex[i])
                    d_year[1, i] = getaverage(ndata[1, tgid], tw, tindex[i])
                    d_year[2, i] = getaverage(ndata[2, tgid], tw, tindex[i])
                    d_year[3, i] = getaverage(ndata[3, tgid], tw, tindex[i])
                    d_year[4, i] = getaverage(ndata[4, tgid], tw, tindex[i])
                    d_year[5, i] = getaverage(ndata[5, tgid], tw, tindex[i])
                    d_year[6, i] = getaverage(ndata[6, tgid], tw, tindex[i])

            except KeyError:
                print(f'KeyError, {tindex[i]}')
                d_year[0, i] = netCDF4.default_fillvals['f8']
                d_year[1, i] = netCDF4.default_fillvals['f8']
                d_year[2, i] = netCDF4.default_fillvals['f8']
                d_year[3, i] = netCDF4.default_fillvals['f8']
                d_year[4, i] = netCDF4.default_fillvals['f8']
                d_year[5, i] = netCDF4.default_fillvals['f8']
                d_year[6, i] = netCDF4.default_fillvals['f8']

        date += dt.timedelta(days=1)
        
        mprcp[day, :] = d_year[0, :]
        mtmax[day, :] = d_year[1, :]
        mtmin[day, :] = d_year[2, :]
        msrad[day, :] = d_year[3, :]
        mswe[day, :] = d_year[4, :]
        mvp[day, :] = d_year[5, :]
        mdayl[day, :] = d_year[6, :]

    finalize(odir, dmyear, gdf1, numdays, start_date, wght_id, mprcp, mtmax, mtmin, msrad, mswe, mvp, mdayl)

if __name__ == "__main__":
    main()
