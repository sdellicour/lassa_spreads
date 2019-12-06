#!/bin/bash
#SBATCH --job-name=LASV_L_2_NGA_skygrid
#SBATCH --time=744:00:00
#SBATCH --partition=Lgpu
#SBATCH --gres=gpu:1
#SBATCH --mem-per-cpu=5120

module load beagle-lib/2.1.2-goolf-1.7.20

java -jar beast_1104_prerelease_051118.jar -beagle_gpu -beagle_double -beagle_order 1 -overwrite LASV_L_2_NGA_skygrid.xml
