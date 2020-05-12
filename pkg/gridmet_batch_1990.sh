#! /bin/sh
#SBATCH -J gridmet_etl_1990
#SBATCH -t 1-09:00
#SBATCH -o %j-dm_out.out
#SBATCH -p workq
#SBATCH -A wbeep
#SBATCH -N 1
##SBATCH -c 10
#SBATCH --mail-type=ALL
#SBATCH --mail-user=rmcd@usgs.gov

export PATH="$PATH:$HOME/miniconda3/bin"

source activate ofp_env

python Gridmet_etl.py -i ../../../gridmet_raw/ -o ../../../gm_gf_v11 -y 1990 -w ../tmp_Gridmet_weights_hru_v1_1e.csv > 1990_gm.out &
python Gridmet_etl.py -i ../../../gridmet_raw/ -o ../../../gm_gf_v11 -y 1991 -w ../tmp_Gridmet_weights_hru_v1_1e.csv > 1991_gm.out &
python Gridmet_etl.py -i ../../../gridmet_raw/ -o ../../../gm_gf_v11 -y 1992 -w ../tmp_Gridmet_weights_hru_v1_1e.csv > 1992_gm.out &
python Gridmet_etl.py -i ../../../gridmet_raw/ -o ../../../gm_gf_v11 -y 1993 -w ../tmp_Gridmet_weights_hru_v1_1e.csv > 1993_gm.out &
python Gridmet_etl.py -i ../../../gridmet_raw/ -o ../../../gm_gf_v11 -y 1994 -w ../tmp_Gridmet_weights_hru_v1_1e.csv > 1994_gm.out &
# python Gridmet_etl.py -i ../../../gridmet_raw/ -o ../../../gm_gf_v11 -y 1995 -w ../tmp_Gridmet_weights_hru_v1_1e.csv > 1995_gm.out &
# python Gridmet_etl.py -i ../../../gridmet_raw/ -o ../../../gm_gf_v11 -y 1996 -w ../tmp_Gridmet_weights_hru_v1_1e.csv > 1996_gm.out &
# python Gridmet_etl.py -i ../../../gridmet_raw/ -o ../../../gm_gf_v11 -y 1997 -w ../tmp_Gridmet_weights_hru_v1_1e.csv > 1997_gm.out &
# python Gridmet_etl.py -i ../../../gridmet_raw/ -o ../../../gm_gf_v11 -y 1998 -w ../tmp_Gridmet_weights_hru_v1_1e.csv > 1998_gm.out &
# python Gridmet_etl.py -i ../../../gridmet_raw/ -o ../../../gm_gf_v11 -y 1999 -w ../tmp_Gridmet_weights_hru_v1_1e.csv > 1999_gm.out &
wait

