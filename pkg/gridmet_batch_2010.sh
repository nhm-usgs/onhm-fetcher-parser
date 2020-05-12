#! /bin/sh
#SBATCH -J gridmet_etl_2010
#SBATCH -t 24:00:00
#SBATCH -o %j-dm_out.out
#SBATCH -p workq
#SBATCH -A wbeep
#SBATCH -N 1
##SBATCH -c 10
#SBATCH --mail-type=ALL
#SBATCH --mail-user=rmcd@usgs.gov

export PATH="$PATH:$HOME/miniconda3/bin"

source activate ofp_env

python Gridmet_etl.py -i ../../../gridmet_raw/ -o ../../../gm_gf_v11 -y 2010 -w ../tmp_Gridmet_weights_hru_v1_1e.csv > 2010_gm.out &
python Gridmet_etl.py -i ../../../gridmet_raw/ -o ../../../gm_gf_v11 -y 2011 -w ../tmp_Gridmet_weights_hru_v1_1e.csv > 2011_gm.out &
python Gridmet_etl.py -i ../../../gridmet_raw/ -o ../../../gm_gf_v11 -y 2012 -w ../tmp_Gridmet_weights_hru_v1_1e.csv > 2012_gm.out &
python Gridmet_etl.py -i ../../../gridmet_raw/ -o ../../../gm_gf_v11 -y 2013 -w ../tmp_Gridmet_weights_hru_v1_1e.csv > 2013_gm.out &
python Gridmet_etl.py -i ../../../gridmet_raw/ -o ../../../gm_gf_v11 -y 2014 -w ../tmp_Gridmet_weights_hru_v1_1e.csv > 2014_gm.out &
# python Gridmet_etl.py -i ../../../gridmet_raw/ -o ../../../gm_gf_v11 -y 2015 -w ../tmp_Gridmet_weights_hru_v1_1e.csv > 2015_gm.out &
# python Gridmet_etl.py -i ../../../gridmet_raw/ -o ../../../gm_gf_v11 -y 2016 -w ../tmp_Gridmet_weights_hru_v1_1e.csv > 2016_gm.out &
# python Gridmet_etl.py -i ../../../gridmet_raw/ -o ../../../gm_gf_v11 -y 2017 -w ../tmp_Gridmet_weights_hru_v1_1e.csv > 2017_gm.out &
# python Gridmet_etl.py -i ../../../gridmet_raw/ -o ../../../gm_gf_v11 -y 2018 -w ../tmp_Gridmet_weights_hru_v1_1e.csv > 2018_gm.out &
# python Gridmet_etl.py -i ../../../gridmet_raw/ -o ../../../gm_gf_v11 -y 2019 -w ../tmp_Gridmet_weights_hru_v1_1e.csv > 2019_gm.out &
wait

