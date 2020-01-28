import geopandas as gpd
import pandas as pd
import netCDF4
import numpy as np
import glob
# import rasterio
import os
import sys
import xarray as xr
#import json
#from rasterstats import zonal_stats

# from rasterio.transform import from_origin
from helper import np_get_wval, get_gm_url, getxml
import requests
from requests.exceptions import HTTPError
from datetime import datetime, timedelta
from pathlib import Path

class FpoNHM:
    """ Class for fetching climate data and parsing into netcdf
        input files for use with the USGS operational National Hydrologic
        Model (oNHM).  Workflow:
            1) Initialize(): fetch climate data
            2) Run(): map/interpolate onto hru
            3) Finalize(): write netcdf input files
        Mapping options:
            1) weighted average based on intersection area of hru
                with netcdf file cells.
            2) rasterstats - zonal averaging

    """
    def __init__(self, climsource='GridMetSS'):
        """
        Initialize class

        :param  numdays: number of days past to retrieve
        :param climsource: Constant for now but may have multiple
            choice for data sources in the future.  Currently default is
            GridMet:  http://www.climatologylab.org/gridmet.html
        """
        self.climsource = climsource
        if climsource == 'GridMetSS':
            self.gmss_vars = {
                'tmax': 'daily_maximum_temperature',
                'tmin': 'daily_minimum_temperature',
                'ppt': 'precipitation_amount',
                'rhmax': 'daily_maximum_relative_humidity',
                'rhmin': 'daily_minimum_relative_humidity',
                'ws': 'daily_mean_wind_speed'}
        # type of retrieval (days) retrieve by previous number of days - used in operational mode
        # or (date) used to retrieve specific period of time
        self.type = None

        self.numdays = None

        #prefix for file names - default is ''.
        self.fileprefix = None

        # xarray containers for tempurature max, temperature min and precipitation
        self.dstmax = None
        self.dstmin = None
        self.dsppt = None
        self.dsrhmax = None
        self.dsrhmin = None
        self.dsws = None

        # Geopandas dataframe that will hold hru id and geometry
        self.gdf = None

        # input and output path directories
        self.iptpath = None
        self.optpath = None

        # weights file
        self.wghts_file = None

        # start and end dates of using type == date
        self.start_date = None
        self.end_date = None

        # handles to netcdf climate data
        # Coordinates
        self.lat_h = None
        self.lon_h = None
        self.time_h = None
        # Geotransform
        self.crs_h = None
        # Climate data
        self.tmax_h = None
        self.tmin_h = None
        self.tppt_h = None
        self.rhmax_h = None
        self.rhmin_h = None
        self.ws_h = None
        # Dimensions
        self.dayshape = None
        self.lonshape = None
        self.latshape = None

        # num HRUs
        self.num_hru = None

        # grouby hru_id_nat on wieghts file
        self.unique_hru_ids = None

        # numpy arrays to store mapped climate data
        self.np_tmax = None
        self.np_tmin = None
        self.np_ppt = None
        self.np_rhmax = None
        self.np_rhmin = None
        self.np_ws = None

        # logical use_date
        self.use_date = False

        # Starting date based on numdays
        self.str_start = None

    def initialize(self, iptpath, optpath, weights_file, type=None, days=None,
                   start_date=None, end_date=None, fileprefix=''):
        """
        Initialize the fp_ohm class:
            1) initialize geopandas dataframe of concatenated hru_shapefiles
            2) initialize climate data using xarray

        :param iptpath: directory containing hru shapefiles and weight file,
                        geotiffs if using rasterstats
        :param optpath: directory to save netcdf input files
        :return: success or failure
        """
        self.iptpath = Path(iptpath)
        if self.iptpath.exists():
            print('input path exits')
        else:
            print('input path does not exist')

        self.optpath = Path(optpath)
        if self.iptpath.exists():
            print('output path exits')
        else:
            print('output path does not exist')

        self.wghts_file = Path(weights_file)
        if self.iptpath.exists():
            print('weights file exits')
        else:
            print('weights file not exist')
        self.wghts_id = None
        self.type = type
        self.numdays = days
        self.start_date = start_date
        self.end_date = end_date
        self.fileprefix = fileprefix
        # #check tmax for available dates
        # now = datetime.today().date()
        # yesterday = now - timedelta(days=1)
        # tdata = getxml()['gridDataset']['TimeSpan']['end']
        # tenddate = datetime.strptime(tdata[:10],'%Y-%m-%d').date()
        # if self.type == 'date':
        #     if self.end_date.date() <= tenddate:
        #         print('Gridmet is available for time period specified - proceeding with pull')
        #     else:
        #         print(f'Requested end date {self.end_date.date()} is greater than available Gridmet end date {tenddate}')
        #         print('process exiting')
        #         sys.exit(1)
        # # else:
        # #     print(f'Requested end date {self.end_date.date()} is greater than yesterdays (Gridment updated to yesterday) date {yesterday}')
        # #     print('process exiting')
        # #     sys.exit(1)
        # else:
        #     if yesterday <= tenddate:
        #         print('Gridmet is updated - proceeding with pull')
        #     else:
        #         print(f'Requested end date {yesterday} is greater than available Gridmet end date {tenddate}')
        #         print('process exiting')
        #         sys.exit(1)

        print(os.getcwd())
        os.chdir(self.iptpath)
        print(os.getcwd())
        if self.type == 'date':
            print(f'start_date: {self.start_date} and end_date: {self.end_date}')
        else:
            print(f'number of days: {self.numdays}')
        # glob.glob produces different results on Win and Linux. Adding sorted makes result consistent
        filenames = sorted(glob.glob('*.shp'))
        self.gdf = pd.concat([gpd.read_file(f) for f in filenames], sort=True).pipe(gpd.GeoDataFrame)
        self.gdf.reset_index(drop=True, inplace=True)
        print(filenames)
        print(self.gdf.head())

        self.num_hru = len(self.gdf.index)
        tmaxfile = None
        tminfile = None
        pptfile = None
        rhminfile = None
        rhmaxfile = None
        wsfile = None
        str_start = None
        if self.type == 'date':
            self.numdays = ((self.end_date - self.start_date).days + 1)
        # Download netcdf subsetted data
        #get_gm_url(type, dataset, numdays=None, startdate=None, enddate=None,  ctype='GridMetSS'):
        try:
            #Maximum Temperature
            self.str_start, tmxurl, tmxparams = get_gm_url(self.type, 'tmax', self.numdays,
                                                           self.start_date, self.end_date)
            tmaxfile = requests.get(tmxurl, params=tmxparams)
            tmaxfile.raise_for_status()
            #Minimum Temperature
            self.str_start, tmnurl, tmnparams = get_gm_url(self.type, 'tmin', self.numdays,
                                                        self.start_date, self.end_date)
            tminfile = requests.get(tmnurl, params=tmnparams)
            tminfile.raise_for_status()
            #Precipitation
            self.str_start, ppturl, pptparams = get_gm_url(self.type, 'ppt', self.numdays,
                                                           self.start_date, self.end_date)
            pptfile = requests.get(ppturl, params=pptparams)
            pptfile.raise_for_status()
            #Maximum Relative Humidity
            self.str_start, rhmaxurl, rhmaxparams = get_gm_url(self.type, 'rhmax', self.numdays,
                                                           self.start_date, self.end_date)
            rhmaxfile = requests.get(rhmaxurl, params=rhmaxparams)
            rhmaxfile.raise_for_status()
            #Minimum Relative Humidity
            self.str_start, rhminurl, rhminparams = get_gm_url(self.type, 'rhmin', self.numdays,
                                                           self.start_date, self.end_date)
            rhminfile = requests.get(rhminurl, params=rhminparams)
            rhminfile.raise_for_status()
            #Mean daily Wind Speed
            self.str_start, wsurl, wsparams = get_gm_url(self.type, 'ws', self.numdays,
                                                           self.start_date, self.end_date)
            wsfile = requests.get(wsurl, params=wsparams)
            wsfile.raise_for_status()

        except HTTPError as http_err:
            print('HTTP error occured: {http_err}')
            if self.numdays == 1:
                sys.exit("numdays == 1: Gridmet not updated")
            else:
                sys.exit("GridMet not available or a bad request")
        except Exception as err:
            print('Other error occured: {err}')
        else:
            print('Gridmet data retrieved!')

        # write downloaded data to local netcdf files and open as xarray
        ncfile = (self.fileprefix + 'tmax_' + (datetime.now().strftime('%Y_%m_%d')) + '.nc',
                  self.fileprefix + 'tmin_' + str(datetime.now().strftime('%Y_%m_%d')) + '.nc',
                  self.fileprefix + 'ppt_' + str(datetime.now().strftime('%Y_%m_%d')) + '.nc',
                  self.fileprefix + 'rhmax_' + str(datetime.now().strftime('%Y_%m_%d')) + '.nc',
                  self.fileprefix + 'rhmin_' + str(datetime.now().strftime('%Y_%m_%d')) + '.nc',
                  self.fileprefix + 'ws_' + str(datetime.now().strftime('%Y_%m_%d')) + '.nc',)

        for index, tfile in enumerate(ncfile):
            with open(tfile, 'wb') as fh:
                if index == 0:
                    fh.write(tmaxfile.content)
                elif index == 1:
                    fh.write(tminfile.content)
                elif index == 2:
                    fh.write(pptfile.content)
                elif index == 3:
                    fh.write(rhmaxfile.content)
                elif index == 4:
                    fh.write(rhminfile.content)
                elif index == 5:
                    fh.write(wsfile.content)

            fh.close()
            if index == 0:
                self.dstmax = xr.open_dataset(tfile)
            elif index == 1:
                self.dstmin = xr.open_dataset(tfile)
            elif index == 2:
                self.dsppt = xr.open_dataset(tfile)
            elif index == 3:
                self.dsrhmax = xr.open_dataset(tfile)
            elif index == 4:
                self.dsrhmin = xr.open_dataset(tfile)
            elif index == 5:
                self.dsws = xr.open_dataset(tfile)

        # =========================================================
        # Get handles to shape/Lat/Lon/DataArray
        # All the datahandles including shape, lat, lon should be
        # the same for each netcdf file.  In the future this may change
        # but for now will open and assume these data handles are the
        # same for each of the climate netcdf files, so grab them from dstmax.
        # =========================================================

        self.lat_h = self.dstmax['lat']
        self.lon_h = self.dstmax['lon']
        self.time_h = self.dstmax['day']
        # self.crs_h = self.dstmax['crs']

        if self.climsource == 'GridMetSS':
            self.tmax_h = self.dstmax[self.gmss_vars['tmax']]
            self.tmin_h = self.dstmin[self.gmss_vars['tmin']]
            self.tppt_h = self.dsppt[self.gmss_vars['ppt']]
            self.rhmax_h = self.dsrhmax[self.gmss_vars['rhmax']]
            self.rhmin_h = self.dsrhmin[self.gmss_vars['rhmin']]
            self.ws_h = self.dsws[self.gmss_vars['ws']]
        else:
            print('Error: climate source data not specified')

        ts = self.tmax_h.sizes
        self.dayshape = ts['day']
        self.lonshape = ts['lon']
        self.latshape = ts['lat']

        # if self.type == 'days':
        print(f'Gridmet returned days = {self.dayshape} and expected number of days {self.numdays}')
        if self.dayshape == self.numdays:
            return True
        else:
            print('returned and expected days not equal')
            return False
        # else:
        #     if self.dayshape == self.numdays:
        #         return True
        #     else:
        #         return False


    def run_weights(self):

        # =========================================================
        #       Read hru weights
        # =========================================================

        wght_uofi = pd.read_csv(self.wghts_file)
        self.wghts_id = wght_uofi.columns[1]
        self.unique_hru_ids = wght_uofi.groupby(self.wghts_id)
        print('finished reading weight file')

        # intialize numpy arrays to store climate vars
        self.np_tmax = np.zeros((self.numdays, self.num_hru))
        self.np_tmin = np.zeros((self.numdays, self.num_hru))
        self.np_ppt = np.zeros((self.numdays, self.num_hru))
        self.np_rhmax = np.zeros((self.numdays, self.num_hru))
        self.np_rhmin = np.zeros((self.numdays, self.num_hru))
        self.np_ws = np.zeros((self.numdays, self.num_hru))

        for day in np.arange(self.numdays):
            print(day)
            tmax = np.zeros(self.num_hru)
            tmin = np.zeros(self.num_hru)
            ppt = np.zeros(self.num_hru)
            rhmax = np.zeros(self.num_hru)
            rhmin = np.zeros(self.num_hru)
            ws = np.zeros(self.num_hru)
            # =========================================================
            #   flatten data arrays
            # =========================================================
            # tmax_h_flt = np.nan_to_num(self.tmax_h.values[day, :, :].flatten(order='K'))
            # tmin_h_flt = np.nan_to_num(self.tmin_h.values[day, :, :].flatten(order='K'))
            # tppt_h_flt = np.nan_to_num(self.tppt_h.values[day, :, :].flatten(order='K'))
            # trhmax_h_flt = np.nan_to_num(self.rhmax_h.values[day, :, :].flatten(order='K'))
            # trhmin_h_flt = np.nan_to_num(self.rhmin_h.values[day, :, :].flatten(order='K'))
            # tws_h_flt = np.nan_to_num(self.ws_h.values[day, :, :].flatten(order='K'))

            tmax_h_flt = self.tmax_h.values[day, :, :].flatten(order='K')
            tmin_h_flt = self.tmin_h.values[day, :, :].flatten(order='K')
            tppt_h_flt = self.tppt_h.values[day, :, :].flatten(order='K')
            trhmax_h_flt = self.rhmax_h.values[day, :, :].flatten(order='K')
            trhmin_h_flt = self.rhmin_h.values[day, :, :].flatten(order='K')
            tws_h_flt = self.ws_h.values[day, :, :].flatten(order='K')

            for index, row in self.gdf.iterrows():
                # weight_id_rows = wght_df_40.loc[wght_df_40[self.wghts_id] == row[self.wghts_id]]
                weight_id_rows = self.unique_hru_ids.get_group(row[self.wghts_id])
                tmax[index] = np.nan_to_num(np_get_wval(tmax_h_flt, weight_id_rows, index+1) - 273.5)
                tmin[index] = np.nan_to_num(np_get_wval(tmin_h_flt, weight_id_rows, index+1) - 273.5)
                ppt[index] = np.nan_to_num(np_get_wval(tppt_h_flt, weight_id_rows, index+1))
                rhmax[index] = np.nan_to_num(np_get_wval(trhmax_h_flt, weight_id_rows, index+1))
                rhmin[index] = np.nan_to_num(np_get_wval(trhmin_h_flt, weight_id_rows, index+1))
                ws[index] = np.nan_to_num(np_get_wval(tws_h_flt, weight_id_rows, index+1))

                if index % 10000 == 0:
                    print(index, row[self.wghts_id])

            self.np_tmax[day, :] = tmax
            self.np_tmin[day, :] = tmin
            self.np_ppt[day, :] = ppt
            self.np_rhmax[day, :] = rhmax
            self.np_rhmin[day, :] = rhmin
            self.np_ws[day, :] = ws

        # close xarray datasets
        self.dstmax.close()
        self.dstmin.close()
        self.dsppt.close()
        self.dsrhmax.close()
        self.dsrhmin.close()
        self.dsws.close()

    def run_rasterstat(self):
        tmp = 0

    def finalize(self):
        print(os.getcwd())
        os.chdir(self.optpath)
        print(os.getcwd())
        ncfile = netCDF4.Dataset(self.fileprefix + 'climate_' + str(datetime.now().strftime('%Y_%m_%d')) + '.nc',
                                 mode='w', format='NETCDF4_CLASSIC')

        # Global Attributes
        ncfile.Conventions = 'CF-1.8'
        ncfile.featureType = 'timeSeries'
        ncfile.history = ''

        sp_dim = len(self.gdf.index)

        hruid_dim = ncfile.createDimension('hruid', sp_dim)  # hru_id
        time_dim = ncfile.createDimension('time', self.numdays)  # unlimited axis (can be appended to).
        for dim in ncfile.dimensions.items():
            print(dim)

        # Create Variables
        time = ncfile.createVariable('time', 'i', ('time',))
        time.long_name = 'time'
        time.standard_name = 'time'
        time.units = 'days since ' + self.str_start

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
        prcp.long_name = 'Daily precipitation rate'
        prcp.units = 'mm/day'
        prcp.standard_name = 'lwe_precipitation_rate'

        tmax = ncfile.createVariable('tmax', np.dtype(np.float32).char, ('time', 'hruid'))
        tmax.long_name = 'Maximum daily air temperature'
        tmax.units = 'degree_Celsius'
        tmax.standard_name = 'maximum_daily_air_temperature'

        tmin = ncfile.createVariable('tmin', np.dtype(np.float32).char, ('time', 'hruid'))
        tmin.long_name = 'Minimum daily air temperature'
        tmin.units = 'degree_Celsius'
        tmin.standard_name = 'minimum_daily_air_temperature'

        rhmax = ncfile.createVariable('rhmax', np.dtype(np.float32).char, ('time', 'hruid'))
        rhmax.long_name = 'Maximum daily relative humidity'
        rhmax.units = 'percent'
        rhmax.standard_name = 'daily_maximum_relative_humidity'

        rhmin = ncfile.createVariable('rhmin', np.dtype(np.float32).char, ('time', 'hruid'))
        rhmin.long_name = 'Minimum daily relative humidity'
        rhmin.units = 'percent'
        rhmin.standard_name = 'daily_minimum_relative_humidity'

        ws = ncfile.createVariable('ws', np.dtype(np.float32).char, ('time', 'hruid'))
        ws.long_name = 'Mean daily wind speed'
        ws.units = 'm/s'
        ws.standard_name = 'daily_mean_wind_speed'


        # fill variables with available data
        def getxy(pt):
            return pt.x, pt.y

        centroidseries = self.gdf.geometry.centroid
        tlon, tlat = [list(t) for t in zip(*map(getxy, centroidseries))]
        # print(lon, lat)
        time[:] = np.arange(0, self.numdays)
        lon[:] = tlon
        lat[:] = tlat
        hru[:] = self.gdf[self.wghts_id].values
        # print(hruid)
        # tmax[0,:] = gdf['tmax'].values
        # tmin[0,:] = gdf['tmin'].values
        # prcp[0,:] = gdf['ppt'].values

        tmax[:, :] = self.np_tmax[:, :]
        tmin[:, :] = self.np_tmin[:, :]
        prcp[:, :] = self.np_ppt[:, :]
        rhmax[:, :] = self.np_rhmax[:, :]
        rhmin[:, :] = self.np_rhmin[:, :]
        ws[:, :] = self.np_ws[:, :]

        ncfile.close()
        print("dataset is closed")

    def setNumdays(self, num_d):
        self.numdays = num_d
