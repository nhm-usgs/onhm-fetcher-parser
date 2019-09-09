import rasterio
from rasterio.transform import from_origin
import numpy as np
import datetime as dt
from numpy.ma import masked


def gridmet_nc_to_geotiff(ds, time_index, path, filename, dsname):

    # collect data to describe geotransform
    lonmin = float(ds.attrs['geospatial_lon_min'])
    latmax = float(ds.attrs['geospatial_lat_max'])
    lonres = float(ds.attrs['geospatial_lon_resolution'])
    latres = float(ds.attrs['geospatial_lon_resolution'])

    # get data shape
    ts = ds.sizes
    dayshape = ts['day']
    lonshape = ts['lon']
    latshape = ts['lat']

    # get transform
    transform = from_origin(lonmin, latmax, lonres, latres)

    # create geodiff from netcdf data
    file = path + filename
    new_dataset = rasterio.open(file, 'w',
                                driver='GTiff',
                                height=lonshape, width=latshape,
                                count=1, dtype=str(ds.dtype),
                                crs={'init': 'epsg:4326'},
                                transform=transform)
    vals = ds[dsname].values[dayshape - time_index, :, :]
    ud = np.flipud(vals)
    new_dataset.write(ud, 1)
    new_dataset.close()


def np_get_wval(ndata, wghts, hru_id):
    """
    Returns weighted average of ndata with weights = grp
    :param ndata: float array of data values
    :param wghts: float array of weights
    :return: numpy weighted averaged - masked to deal with nans associated with
            ndata that is outside of the conus.
    """
    mdata = np.ma.masked_array(ndata[wghts['grid_ids'].values.astype(int)],
                               np.isnan(ndata[wghts['grid_ids'].values.astype(int)]))

    # mdata = np.ma.masked_where(ndata[wghts['grid_ids'].values.astype(int)] <= 0.0,
    #                            (ndata[wghts['grid_ids'].values.astype(int)]))
    tmp = np.ma.average(mdata, weights=wghts['w'])
    if tmp is masked:
        print('returning masked value', hru_id, mdata, wghts['w'])
        return 0.0

    else:
        return tmp


def get_gm_url(type, dataset, numdays=None, startdate=None, enddate=None,  ctype='GridMetSS'):
    """
    This helper function returns a url and payload to be used with requests
    to get climate data.  Returned values can be used in a request for example:
    myfile = requests.get(url, params=payload)

    :param numdays: proceeding number of days to retrieve
    :param dataset: datset to retrieve can be:
        'tmax', 'tmin', 'ppt'
    :param ctype: Type of url to retrieve:
        'GridMet':
    :return: URL for retrieving GridMet subset data and payload of options
    """
    sformat = "%Y-%m-%d"
    if type == 'days':
        dt1 = dt.timedelta(days=1) # because Gridmet data release today is yesterdays data
        dt2 = dt.timedelta(days=numdays)

        end = dt.datetime.now() - dt1
        start = dt.datetime.now() - dt2
        str_start = start.strftime(sformat) + "T00:00:00Z"
        str_end = end.strftime(sformat) + "T00:00:00Z"
        str_start_cf = start.strftime(sformat) + " 00:00:00"
    elif type == 'date':
        str_start = startdate.strftime(sformat) + "T00:00:00Z"
        str_end = enddate.strftime(sformat) + "T00:00:00Z"
        str_start_cf = startdate.strftime(sformat) + " 00:00:00"

    dsvar = None
    url = None
    if ctype == 'GridMetSS': #extract data using NetcdfSubset service
        if dataset == 'tmax':
            dsvar = 'daily_maximum_temperature'
            url = 'http://thredds.northwestknowledge.net:8080/thredds/ncss/agg_met_tmmx_1979_CurrentYear_CONUS.nc'
        elif dataset == 'tmin':
            dsvar = 'daily_minimum_temperature'
            url = 'http://thredds.northwestknowledge.net:8080/thredds/ncss/agg_met_tmmn_1979_CurrentYear_CONUS.nc'
        elif dataset == 'ppt':
            dsvar = 'precipitation_amount'
            url = 'http://thredds.northwestknowledge.net:8080/thredds/ncss/agg_met_pr_1979_CurrentYear_CONUS.nc'
        elif dataset == 'rhmax':
            dsvar = 'daily_maximum_relative_humidity'
            url = 'http://thredds.northwestknowledge.net:8080/thredds/ncss/grid/agg_met_rmax_1979_CurrentYear_CONUS.nc'
        elif dataset == 'rhmin':
            dsvar = 'daily_minimum_relative_humidity'
            url = 'http://thredds.northwestknowledge.net:8080/thredds/ncss/grid/agg_met_rmin_1979_CurrentYear_CONUS.nc'
        elif dataset == 'ws':
            dsvar = 'daily_mean_wind_speed'
            url = 'http://thredds.northwestknowledge.net:8080/thredds/ncss/grid/agg_met_vs_1979_CurrentYear_CONUS.nc'
        payload = {
            'var': dsvar,
            'north': '49.4000',
            'west': '-124.7666',
            'east': '-67.0583',
            'south': '25.0666',
            'disableLLSubset': 'on',
            'disableProjSubset': 'on',
            'horizStride': '1',
            'time_start': str_start,
            'time_end': str_end,
            'timeStride': '1',
            'accept': 'netcdf'}
        return str_start_cf, url, payload
