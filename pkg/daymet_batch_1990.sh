#! /bin/sh
#SBATCH -J daymet_etl_1990
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

python Daymet_etl.py -i ../../../daymet_v3_raw/north_america/ -o ../../../dm_gf_v11 -y 1990 -w ../daymet_weights_hru_v1_1e.csv > 1990.out &
python Daymet_etl.py -i ../../../daymet_v3_raw/north_america/ -o ../../../dm_gf_v11 -y 1991 -w ../daymet_weights_hru_v1_1e.csv > 1991.out &
python Daymet_etl.py -i ../../../daymet_v3_raw/north_america/ -o ../../../dm_gf_v11 -y 1992 -w ../daymet_weights_hru_v1_1e.csv > 1992.out &
python Daymet_etl.py -i ../../../daymet_v3_raw/north_america/ -o ../../../dm_gf_v11 -y 1993 -w ../daymet_weights_hru_v1_1e.csv > 1993.out &
python Daymet_etl.py -i ../../../daymet_v3_raw/north_america/ -o ../../../dm_gf_v11 -y 1994 -w ../daymet_weights_hru_v1_1e.csv > 1994.out &
python Daymet_etl.py -i ../../../daymet_v3_raw/north_america/ -o ../../../dm_gf_v11 -y 1995 -w ../daymet_weights_hru_v1_1e.csv > 1995.out &
python Daymet_etl.py -i ../../../daymet_v3_raw/north_america/ -o ../../../dm_gf_v11 -y 1996 -w ../daymet_weights_hru_v1_1e.csv > 1996.out &
python Daymet_etl.py -i ../../../daymet_v3_raw/north_america/ -o ../../../dm_gf_v11 -y 1997 -w ../daymet_weights_hru_v1_1e.csv > 1997.out &
python Daymet_etl.py -i ../../../daymet_v3_raw/north_america/ -o ../../../dm_gf_v11 -y 1998 -w ../daymet_weights_hru_v1_1e.csv > 1998.out &
python Daymet_etl.py -i ../../../daymet_v3_raw/north_america/ -o ../../../dm_gf_v11 -y 1999 -w ../daymet_weights_hru_v1_1e.csv > 1999.out &
wait

