#!/bin/bash
#SBATCH -N 1
#SBATCH -C knl
#SBATCH -q interactive 
#SBATCH -J test
#SBATCH --mail-user=fernando_tta@berkeley.edu
#SBATCH --mail-type=ALL
#SBATCH -t 02:30:00

#OpenMP settings:
export OMP_NUM_THREADS=272
export OMP_PLACES=threads
export OMP_PROC_BIND=spread


#run the application:
srun -n 1 -c 272 --cpu_bind=cores myapp.x
#change to salloc for interactive SHELL queue
