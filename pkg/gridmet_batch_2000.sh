#! /bin/sh
#SBATCH -J gridmet_etl_2000
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

python Gridmet_etl.py -i ../../../gridmet_raw/ -o ../../../gm_gf_v11 -y 2000 -w ../tmp_Gridmet_weights_hru_v1_1e.csv > 2000_gm.out &
python Gridmet_etl.py -i ../../../gridmet_raw/ -o ../../../gm_gf_v11 -y 2001 -w ../tmp_Gridmet_weights_hru_v1_1e.csv > 2001_gm.out &
python Gridmet_etl.py -i ../../../gridmet_raw/ -o ../../../gm_gf_v11 -y 2002 -w ../tmp_Gridmet_weights_hru_v1_1e.csv > 2002_gm.out &
python Gridmet_etl.py -i ../../../gridmet_raw/ -o ../../../gm_gf_v11 -y 2003 -w ../tmp_Gridmet_weights_hru_v1_1e.csv > 2003_gm.out &
python Gridmet_etl.py -i ../../../gridmet_raw/ -o ../../../gm_gf_v11 -y 2004 -w ../tmp_Gridmet_weights_hru_v1_1e.csv > 2004_gm.out &
# python Gridmet_etl.py -i ../../../gridmet_raw/ -o ../../../gm_gf_v11 -y 2005 -w ../tmp_Gridmet_weights_hru_v1_1e.csv > 2005_gm.out &
# python Gridmet_etl.py -i ../../../gridmet_raw/ -o ../../../gm_gf_v11 -y 2006 -w ../tmp_Gridmet_weights_hru_v1_1e.csv > 2006_gm.out &
# python Gridmet_etl.py -i ../../../gridmet_raw/ -o ../../../gm_gf_v11 -y 2007 -w ../tmp_Gridmet_weights_hru_v1_1e.csv > 2007_gm.out &
# python Gridmet_etl.py -i ../../../gridmet_raw/ -o ../../../gm_gf_v11 -y 2008 -w ../tmp_Gridmet_weights_hru_v1_1e.csv > 2008_gm.out &
# python Gridmet_etl.py -i ../../../gridmet_raw/ -o ../../../gm_gf_v11 -y 2009 -w ../tmp_Gridmet_weights_hru_v1_1e.csv > 2009_gm.out &
wait

