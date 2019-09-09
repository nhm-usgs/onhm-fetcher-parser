#!/bin/bash
python ../pkg/climate_etl.py -t date -p 2015-01-01 2015-12-31 -f 2015_ -i ../Data -o ../Output -w ../Data/weights.csv &
python ../pkg/climate_etl.py -t date -p 2016-01-01 2016-12-31 -f 2016_ -i ../Data -o ../Output -w ../Data/weights.csv &
python ../pkg/climate_etl.py -t date -p 2017-01-01 2017-12-31 -f 2017_ -i ../Data -o ../Output -w ../Data/weights.csv  &
python ../pkg/climate_etl.py -t date -p 2018-01-01 2018-12-31 -f 2018_ -i ../Data -o ../Output -w ../Data/weights.csv &
python ../pkg/climate_etl.py -t date -p 2019-01-01 2019-07-31 -f 2019_ -i ../Data -o ../Output -w ../Data/weights.csv 
