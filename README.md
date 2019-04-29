# onhm-fetcher-parser

Experiments in fetching and parsing  climate data for oNHM

setup environment using conda

	* conda create --name gridmet_3 python=3.6 netcdf4 numpy pandas xarray dask bottleneck jupyter notebook=5.7.4 matplotlib cartopy scipy
	* conda activate gridmet_3
	* conda install -c pyviz geoviews
	* conda install rasterio
	* conda install -c conda-forge rasterstats
	
* Add Data folder to root directory and get data from:
* https://drive.google.com/open?id=1WMu_sAceVZv7BrovudtOG_whTaIjZUgT
* unzip Data.7z into Data folder - contains shapefiles of hru by region
