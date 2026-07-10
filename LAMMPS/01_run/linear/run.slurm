#!/bin/bash
#SBATCH --job-name=tree_n20_seed101
#SBATCH --output=job.out
#SBATCH --error=job.err
#SBATCH --ntasks=1
#SBATCH --time=08:00:00
#SBATCH --partition=phd_student
#SBATCH --qos=phd_student

source /apps/compilers/intel/oneapi/setvars.sh
module load lammps-cpu
export I_MPI_FABRICS=shm

SEED=101
NSTEPS=2000000

/apps/codes/lammps/cpu/bin/lmp_oneapi -in in.poly -var seed ${SEED} -var NSTEPS ${NSTEPS}
