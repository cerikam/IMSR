#!/bin/bash

#PBS -V
#PBS -q fill
#PBS -l nodes=1:ppn=64

cd $PBS_O_WORKDIR

hostname

module unload mpi
module load openmpi/2.1.6
module load scale/6.3.1-mpi
export DATA=/opt/scale6.3_data

export HDF5_USE_FILE_LOCKING=FALSE

scalerte -m -N 32 EIRENE.inp

