#!/bin/bash
#PBS -V
#PBS -N job_SIRENE
#PBS -q fill
#PBS -l nodes=1:ppn=64

hostname
rm -f done.dat
cd ${PBS_O_WORKDIR}
module load mpi
module load serpent

sss2 -omp $PBS_NUM_PPN SIRENE > myout.out
awk 'BEGIN{ORS="\t"} /ANA_KEFF/ || /CONVERSION/ {print $7" "$8;}' SIRENE_res.m > done.out

