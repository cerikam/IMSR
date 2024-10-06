#!/bin/bash

#PBS -V
#PBS -q amd
#PBS -l nodes=1:ppn=64

cd $PBS_O_WORKDIR

hostname

module load mpi
module load scale/6.3.2-mpi
export DATA=/opt/scale6.3_data

export HDF5_USE_FILE_LOCKING=FALSE

#scalerte -m -N $PBS_NUM_PPN ThEIRENE.inp
scalerte -m -N 32 ThEIRENE.inp

