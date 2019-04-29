import rasterio
from rasterio.transform import from_origin
import numpy as np

def gridmet_nc_to_geotiff(ds, time_index, path, filename, dsname ):

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

    #get transform
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