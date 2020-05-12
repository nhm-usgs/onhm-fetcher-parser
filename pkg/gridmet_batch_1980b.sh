#! /bin/sh
#SBATCH -J gridmet_etl_1980
#SBATCH -t 1-08:00
#SBATCH -o %j-dm_out.out
#SBATCH -p workq
#SBATCH -A wbeep
#SBATCH -N 1
##SBATCH -c 10
#SBATCH --mail-type=ALL
#SBATCH --mail-user=rmcd@usgs.gov

export PATH="$PATH:$HOME/miniconda3/bin"

source activate ofp_env

# python Gridmet_etl.py -i ../../../gridmet_raw/ -o ../../../gm_gf_v11 -y 1979 -w ../tmp_Gridmet_weights_hru_v1_1e.csv > 1979_gm.out &
# python Gridmet_etl.py -i ../../../gridmet_raw/ -o ../../../gm_gf_v11 -y 1980 -w ../tmp_Gridmet_weights_hru_v1_1e.csv > 1980_gm.out &
# python Gridmet_etl.py -i ../../../gridmet_raw/ -o ../../../gm_gf_v11 -y 1981 -w ../tmp_Gridmet_weights_hru_v1_1e.csv > 1981_gm.out &
# python Gridmet_etl.py -i ../../../gridmet_raw/ -o ../../../gm_gf_v11 -y 1982 -w ../tmp_Gridmet_weights_hru_v1_1e.csv > 1982_gm.out &
# python Gridmet_etl.py -i ../../../gridmet_raw/ -o ../../../gm_gf_v11 -y 1983 -w ../tmp_Gridmet_weights_hru_v1_1e.csv > 1983_gm.out &
# python Gridmet_etl.py -i ../../../gridmet_raw/ -o ../../../gm_gf_v11 -y 1984 -w ../tmp_Gridmet_weights_hru_v1_1e.csv > 1984_gm.out &
python Gridmet_etl.py -i ../../../gridmet_raw/ -o ../../../gm_gf_v11 -y 1985 -w ../tmp_Gridmet_weights_hru_v1_1e.csv > 1985_gm.out &
python Gridmet_etl.py -i ../../../gridmet_raw/ -o ../../../gm_gf_v11 -y 1986 -w ../tmp_Gridmet_weights_hru_v1_1e.csv > 1986_gm.out &
python Gridmet_etl.py -i ../../../gridmet_raw/ -o ../../../gm_gf_v11 -y 1987 -w ../tmp_Gridmet_weights_hru_v1_1e.csv > 1987_gm.out &
python Gridmet_etl.py -i ../../../gridmet_raw/ -o ../../../gm_gf_v11 -y 1988 -w ../tmp_Gridmet_weights_hru_v1_1e.csv > 1988_gm.out &
python Gridmet_etl.py -i ../../../gridmet_raw/ -o ../../../gm_gf_v11 -y 1989 -w ../tmp_Gridmet_weights_hru_v1_1e.csv > 1989_gm.out &
wait

