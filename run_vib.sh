#!/bin/bash
#SBATCH -J AIMS_vib
#SBATCH --partition=serial
#SBATCH --time=01:59:59
##SBATCH -N 1
#SBATCH -n 16
#SBATCH --mem-per-cpu=1900
#SBATCH --cpus-per-task=1
#SBATCH -o AIMS_output
#SBATCH -e AIMS_error
#
export KMP_AFFINITY=compact
export OMP_NUM_THREADS=1
export OMP_STACKSIZE="1G"
#
#
#
#
#

aims_dir=/homeappl/home/hyvonen1/fhi-aims.160328_3/bin/aims.160328_3.mpi.x

delta=0.0025
mkdir delta_$delta
cp control.in delta_$delta/
cp geometry.in delta_$delta/
cp get_vibrations.py delta_$delta/
cd delta_$delta

module load python-env

python get_vibrations.py -r $aims_dir run 0 >& vib_$delta.out
python get_vibrations.py run 1 >& vib_post_$delta.out
