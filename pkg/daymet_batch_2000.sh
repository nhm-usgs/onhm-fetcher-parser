#! /bin/sh
#SBATCH -J daymet_etl
#SBATCH -t 1-08:00
#SBATCH -o %j-dm_out.out
#SBATCH -p workq
#SBATCH -A wbeep
#SBATCH -N 1
#SBATCH -c 10
#SBATCH --mail-type=ALL
#SBATCH --mail-user=rmcd@usgs.gov

export PATH="$PATH:$HOME/miniconda3/bin"

source activate ofp_env

python Daymet_etl.py -i ../../../daymet_v3_raw/north_america/ -o ../../../dm_gf_v11 -y 2000 -w ../daymet_weights_hru_v1_1e.csv > 2000.out &
python Daymet_etl.py -i ../../../daymet_v3_raw/north_america/ -o ../../../dm_gf_v11 -y 2001 -w ../daymet_weights_hru_v1_1e.csv > 2001.out &
python Daymet_etl.py -i ../../../daymet_v3_raw/north_america/ -o ../../../dm_gf_v11 -y 2002 -w ../daymet_weights_hru_v1_1e.csv > 2002.out &
python Daymet_etl.py -i ../../../daymet_v3_raw/north_america/ -o ../../../dm_gf_v11 -y 2003 -w ../daymet_weights_hru_v1_1e.csv > 2003.out &
python Daymet_etl.py -i ../../../daymet_v3_raw/north_america/ -o ../../../dm_gf_v11 -y 2004 -w ../daymet_weights_hru_v1_1e.csv > 2004.out &
python Daymet_etl.py -i ../../../daymet_v3_raw/north_america/ -o ../../../dm_gf_v11 -y 2005 -w ../daymet_weights_hru_v1_1e.csv > 2005.out &
python Daymet_etl.py -i ../../../daymet_v3_raw/north_america/ -o ../../../dm_gf_v11 -y 2006 -w ../daymet_weights_hru_v1_1e.csv > 2006.out &
python Daymet_etl.py -i ../../../daymet_v3_raw/north_america/ -o ../../../dm_gf_v11 -y 2007 -w ../daymet_weights_hru_v1_1e.csv > 2007.out &
python Daymet_etl.py -i ../../../daymet_v3_raw/north_america/ -o ../../../dm_gf_v11 -y 2008 -w ../daymet_weights_hru_v1_1e.csv > 2008.out &
python Daymet_etl.py -i ../../../daymet_v3_raw/north_america/ -o ../../../dm_gf_v11 -y 2009 -w ../daymet_weights_hru_v1_1e.csv > 2009.out &
wait

