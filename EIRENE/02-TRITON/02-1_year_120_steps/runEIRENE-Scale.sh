#!/bin/bash

#PBS -V
#PBS -q fill
#PBS -l nodes=c0101:ppn=64

cd $PBS_O_WORKDIR

hostname

module unload mpi
module load openmpi/2.1.6
#module load scale/6.3.b15
module load scale/ondrejch
export DATA=/opt/scale6.3_data

#export scale_input_dir=/home/ondrejch/projects/myMSR/SCALE
export HDF5_USE_FILE_LOCKING=FALSE

scalerte -m -N 64 EIRENE.inp

