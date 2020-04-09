#! /bin/sh
#SBATCH -J daymet_etl
#SBATCH -t 24:00:00
#SBATCH -o %j-dm_out.out
#SBATCH -p workq
#SBATCH -A wbeep
#SBATCH -N 1
#SBATCH -c 10
#SBATCH --mail-type=ALL
#SBATCH --mail-user=rmcd@usgs.gov

export PATH="$PATH:$HOME/miniconda3/bin"

source activate ofp_env

python Daymet_etl.py -i ../../../daymet_v3_raw/north_america/ -o ../../../dm_gf_v11 -y 2010 -w ../daymet_weights_hru_v1_1e.csv > 2010.out &
python Daymet_etl.py -i ../../../daymet_v3_raw/north_america/ -o ../../../dm_gf_v11 -y 2011 -w ../daymet_weights_hru_v1_1e.csv > 2011.out &
python Daymet_etl.py -i ../../../daymet_v3_raw/north_america/ -o ../../../dm_gf_v11 -y 2012 -w ../daymet_weights_hru_v1_1e.csv > 2012.out &
python Daymet_etl.py -i ../../../daymet_v3_raw/north_america/ -o ../../../dm_gf_v11 -y 2013 -w ../daymet_weights_hru_v1_1e.csv > 2013.out &
python Daymet_etl.py -i ../../../daymet_v3_raw/north_america/ -o ../../../dm_gf_v11 -y 2014 -w ../daymet_weights_hru_v1_1e.csv > 2014.out &
python Daymet_etl.py -i ../../../daymet_v3_raw/north_america/ -o ../../../dm_gf_v11 -y 2015 -w ../daymet_weights_hru_v1_1e.csv > 2015.out &
python Daymet_etl.py -i ../../../daymet_v3_raw/north_america/ -o ../../../dm_gf_v11 -y 2016 -w ../daymet_weights_hru_v1_1e.csv > 2016.out &
python Daymet_etl.py -i ../../../daymet_v3_raw/north_america/ -o ../../../dm_gf_v11 -y 2017 -w ../daymet_weights_hru_v1_1e.csv > 2017.out &
python Daymet_etl.py -i ../../../daymet_v3_raw/north_america/ -o ../../../dm_gf_v11 -y 2018 -w ../daymet_weights_hru_v1_1e.csv > 2018.out &
#python Daymet_etl.py -i ../../../daymet_v3_raw/north_america/ -o ../../../dm_gf_v11 -y 2019 -w ../daymet_weights_hru_v1_1e.csv > 2019.out &
wait

