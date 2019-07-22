#!/bin/bash
#SBATCH -J AIMS_vib
#SBATCH --partition=serial
#SBATCH --time=00:59:59
##SBATCH -N 1
#SBATCH -n 4
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
module load python-env


delta=0.0025
mkdir delta_$delta
cp control.in delta_$delta/
cp geometry.in delta_$delta/
cp get_vibrations_occ.py delta_$delta/
cd delta_$delta

# -s suffix when running aims
# -r path to aims binary
python get_vibrations_occ.py -s srun -r $aims_dir run 0 >& vib_$delta.out
python get_vibrations_occ.py run 1 >& vib_post_$delta.out
