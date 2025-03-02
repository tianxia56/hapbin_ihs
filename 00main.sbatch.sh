#!/bin/bash
#SBATCH --partition=ycga
#SBATCH --time=2-00:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=20000

#bash submit_making_bfile.sh

#bash bfile.qc.sh

#Rscript convert_AA_bfile.R

#python convert_rmap.py

#bash sub_make_tped.sh

#python tped_to_hap.py

#bash sub_run_hapbin_ihs.sh 

#python three_norms.py

#python compare_outputs.py

