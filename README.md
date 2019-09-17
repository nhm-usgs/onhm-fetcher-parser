# onhm-fetcher-parser

Experiments in fetching and parsing  climate data for oNHM

setup environment using conda

	* conda env create -f pgk/environment.yml
	* conda activate ofp_env

	
* Add Data folder to root directory and get data from:
* ftp://ftpext.usgs.gov/pub/cr/co/denver/BRR-CR/pub/rmcd/Data_hru_shp_v2.tar.gz
* unzip Data_hru_shp_v2.tar.gz into Data folder - contains shapefiles of hru by region

To run onhm-fetcher-parser from a conda cmd prompt using the ofp_env environment 

	* python pkg/climate_etl.py -t date -p 2015-01-01 2015-12-31 -f 2015_ -i Data -o Output -w Data/weights.csv
