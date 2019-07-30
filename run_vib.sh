#!/bin/bash
#SBATCH -J AIMS_vib
#SBATCH --partition=serial
#SBATCH --time=05:59:59
##SBATCH -N 1
#SBATCH -n 8
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

# suffix when running aims e.g. 'mpirun -np 4'
suffix=srun

# aims binary
aims_bin=/homeappl/home/hyvonen1/fhi-aims.160328_3/bin/aims.160328_3.mpi.x

# modules that must be loaded
module load python-env



delta=0.0025
mkdir delta_$delta
cp control.in delta_$delta/
cp geometry.in delta_$delta/
cp get_vibrations.py delta_$delta/
cd delta_$delta

# -s suffix when running aims
# -r path to aims binary
python get_vibrations.py -s $suffix -r $aims_bin run 0 >& vib_$delta.out
python get_vibrations.py run 1 >& vib_post_$delta.out
