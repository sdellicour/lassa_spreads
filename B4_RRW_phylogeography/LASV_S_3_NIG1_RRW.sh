#!/bin/bash

#PBS -m abe
#PBS -M sdellico@ulb.ac.be
#PBS -q gpu
#PBS -l walltime=120:00:00
#PBS -l nodes=1:ppn=1:gpus=1:gpgpu
#PBS -l feature=kepler
#PBS -o LASV_S_3_NIG1_RRW.out
#PBS -e LASV_S_3_NIG1_RRW.err

# module load Java/1.8.0_92
module load beagle-lib/3.0.2-foss-2018b-CUDA-9.2.88

# cd $PBS_O_WORKDIR

java -jar beast_1104_prerelease_051118.jar -beagle_gpu -beagle_double -beagle_order 1 -overwrite LASV_S_3_NIG1_RRW.xml
