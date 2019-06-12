import geopandas as gpd
import pandas as pd
import netCDF4
import numpy as np
import glob
# import rasterio
import os
import xarray as xr
import json
from rasterstats import zonal_stats

# from rasterio.transform import from_origin
from helper import np_get_wval, get_gm_url
import requests
from requests.exceptions import HTTPError
from datetime import datetime as datetime
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
    def __init__(self, numdays, climsource='GridMetSS'):
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
                'ppt': 'precipitation_amount'}
        self.numdays = numdays
        # xarray containers for tempurature max, temperature min and precipitation
        self.dstmax = None
        self.dstmin = None
        self.dsppt = None

        # Geopandas dataframe that will hold hru id and geometry
        self.gdf = None

        # input and output path directories
        self.iptpath = None
        self.optpath = None

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

        # Starting date based on numdays
        self.str_start = None

    def initialize(self, iptpath, optpath):
        """
        Initialize the fp_ohm class:
            1) initialize geopandas dataframe of concatenated hru_shapefiles
            2) initialize climate data using xarray

        :param iptpath: directory containing hru shapefiles and weight file,
                        geotiffs if using rasterstats
        :param optpath: directory to save netcdf input files
        :return: success or failure
        """
        self.iptpath = iptpath
        self.optpath = optpath
        print(os.getcwd())
        os.chdir(self.iptpath)
        print(os.getcwd())
        filenames = glob.glob('*.shp')
        self.gdf = pd.concat([gpd.read_file(f) for f in filenames]).pipe(gpd.GeoDataFrame)
        self.gdf.reset_index(drop=True, inplace=True)
        print(filenames)

        self.num_hru = len(self.gdf.index)
        tmaxfile = None
        tminfile = None
        pptfile = None
        str_start = None
        # Download netcdf subsetted data
        try:
            self.str_start, tmxurl, tmxparams = get_gm_url(self.numdays, 'tmax')
            tmaxfile = requests.get(tmxurl, params=tmxparams)
            tmaxfile.raise_for_status()

            self.str_start, tmnurl, tmnparams = get_gm_url(self.numdays, 'tmin')
            tminfile = requests.get(tmnurl, params=tmnparams)
            tminfile.raise_for_status()

            self.str_start, ppturl, pptparams = get_gm_url(self.numdays, 'ppt')
            pptfile = requests.get(ppturl, params=pptparams)
            pptfile.raise_for_status()

        except HTTPError as http_err:
            print('HTTP error occured: {http_err}')
        except Exception as err:
            print('Other error occured: {err}')
        else:
            print('Success!')

        # write downloaded data to local netcdf files and open as xarray
        ncfile = ('tmax_' + str(datetime.now().strftime('%Y_%m_%d')) + '.nc',
                  'tmin_' + str(datetime.now().strftime('%Y_%m_%d')) + '.nc',
                  'ppt_' + str(datetime.now().strftime('%Y_%m_%d')) + '.nc')
        for index, tfile in enumerate(ncfile):
            with open(tfile, 'wb') as fh:
                if index == 0:
                    fh.write(tmaxfile.content)
                elif index == 1:
                    fh.write(tminfile.content)
                elif index == 2:
                    fh.write(pptfile.content)
            fh.close()
            if index == 0:
                self.dstmax = xr.open_dataset(tfile)
            elif index == 1:
                self.dstmin = xr.open_dataset(tfile)
            elif index == 2:
                self.dsppt = xr.open_dataset(tfile)

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
        else:
            print('Error: climate source data not specified')

        ts = self.tmax_h.sizes
        self.dayshape = ts['day']
        self.lonshape = ts['lon']
        self.latshape = ts['lat']

        if self.dayshape == self.numdays:
            return True
        else:
            return False

    def run_weights(self):

        # =========================================================
        #       Read hru weights
        # =========================================================

        wght_uofi = pd.read_csv(self.iptpath / Path('hru_uofimetdata_weights.csv'))
        self.unique_hru_ids = wght_uofi.groupby('hru_id_nat')
        print('finished reading weight file')

        # intialize numpy arrays to store climate vars
        self.np_tmax = np.zeros((self.numdays, self.num_hru))
        self.np_tmin = np.zeros((self.numdays, self.num_hru))
        self.np_ppt = np.zeros((self.numdays, self.num_hru))

        for day in np.arange(self.numdays):
            print(day)
            tmax = np.zeros(self.num_hru)
            tmin = np.zeros(self.num_hru)
            ppt = np.zeros(self.num_hru)
            # =========================================================
            #   flatten data arrays
            # =========================================================
            tmax_h_flt = np.nan_to_num(self.tmax_h.values[day, :, :].flatten(order='K'))
            tmin_h_flt = np.nan_to_num(self.tmin_h.values[day, :, :].flatten(order='K'))
            tppt_h_flt = np.nan_to_num(self.tppt_h.values[day, :, :].flatten(order='K'))

            for index, row in self.gdf.iterrows():
                # weight_id_rows = wght_df_40.loc[wght_df_40['hru_id_nat'] == row['hru_id_nat']]
                weight_id_rows = self.unique_hru_ids.get_group(row['hru_id_nat'])
                tmax[index] = np_get_wval(tmax_h_flt, weight_id_rows) - 273.5
                tmin[index] = np_get_wval(tmin_h_flt, weight_id_rows) - 273.5
                ppt[index] = np_get_wval(tppt_h_flt, weight_id_rows)
                if index % 10000 == 0:
                    print(index, row['hru_id_nat'])

            self.np_tmax[day, :] = tmax
            self.np_tmin[day, :] = tmin
            self.np_ppt[day, :] = ppt

        # close xarray datasets
        self.dstmax.close()
        self.dstmin.close()
        self.dsppt.close()

    def run_rasterstat(self):
        tmp = 0

    def finalize(self):
        print(os.getcwd())
        os.chdir(self.optpath)
        print(os.getcwd())
        ncfile = netCDF4.Dataset('climate_' + str(datetime.now().strftime('%Y_%m_%d')) + '.nc',
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

        # fill variables with available data
        def getxy(pt):
            return pt.x, pt.y

        centroidseries = self.gdf['geometry'].centroid
        tlon, tlat = [list(t) for t in zip(*map(getxy, centroidseries))]
        # print(lon, lat)
        time[:] = np.arange(0, self.numdays)
        lon[:] = tlon
        lat[:] = tlat
        hru[:] = self.gdf['hru_id_nat'].values
        # print(hruid)
        # tmax[0,:] = gdf['tmax'].values
        # tmin[0,:] = gdf['tmin'].values
        # prcp[0,:] = gdf['ppt'].values

        tmax[:, :] = self.np_tmax[:, :]
        tmin[:, :] = self.np_tmin[:, :]
        prcp[:, :] = self.np_ppt[:, :]

        ncfile.close()
        print("dataset is closed")

    def setNumdays(self, num_d):
        self.numdays = num_d
