#!/bin/bash
#SBATCH -J FHI-AIMS
#SBATCH --partition=serial
#SBATCH --time=00:19:59
##SBATCH -N 1
#SBATCH -n 8
#SBATCH --mem-per-cpu=1900
#SBATCH --cpus-per-task=1
#SBATCH -o AIMS_output
#SBATCH -e AIMS_error
#
#
#
aims_dir=/homeappl/home/hyvonen1/fhi-aims.160328_3/bin/aims.160328_3.mpi.x
#
#
#
srun $aims_dir >& aims.out 
