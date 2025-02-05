#!/bin/bash
#PBS -V
#PBS -q xeon
#PBS -l nodes=1:ppn=8

hostname
cd ${PBS_O_WORKDIR}
module load openmc

openmc > done.out

