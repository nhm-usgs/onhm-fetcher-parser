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

python Daymet_etl.py -i ../../../daymet_v3_raw/north_america/ -o ../../../dm_gf_v11 -y 1980 -w ../daymet_weights_hru_v1_1e.csv > 1980.out &
python Daymet_etl.py -i ../../../daymet_v3_raw/north_america/ -o ../../../dm_gf_v11 -y 1981 -w ../daymet_weights_hru_v1_1e.csv > 1981.out &
python Daymet_etl.py -i ../../../daymet_v3_raw/north_america/ -o ../../../dm_gf_v11 -y 1982 -w ../daymet_weights_hru_v1_1e.csv > 1982.out &
python Daymet_etl.py -i ../../../daymet_v3_raw/north_america/ -o ../../../dm_gf_v11 -y 1983 -w ../daymet_weights_hru_v1_1e.csv > 1983.out &
python Daymet_etl.py -i ../../../daymet_v3_raw/north_america/ -o ../../../dm_gf_v11 -y 1984 -w ../daymet_weights_hru_v1_1e.csv > 1984.out &
python Daymet_etl.py -i ../../../daymet_v3_raw/north_america/ -o ../../../dm_gf_v11 -y 1985 -w ../daymet_weights_hru_v1_1e.csv > 1985.out &
python Daymet_etl.py -i ../../../daymet_v3_raw/north_america/ -o ../../../dm_gf_v11 -y 1986 -w ../daymet_weights_hru_v1_1e.csv > 1986.out &
python Daymet_etl.py -i ../../../daymet_v3_raw/north_america/ -o ../../../dm_gf_v11 -y 1987 -w ../daymet_weights_hru_v1_1e.csv > 1987.out &
python Daymet_etl.py -i ../../../daymet_v3_raw/north_america/ -o ../../../dm_gf_v11 -y 1988 -w ../daymet_weights_hru_v1_1e.csv > 1988.out &
python Daymet_etl.py -i ../../../daymet_v3_raw/north_america/ -o ../../../dm_gf_v11 -y 1989 -w ../daymet_weights_hru_v1_1e.csv > 1989.out &
wait

