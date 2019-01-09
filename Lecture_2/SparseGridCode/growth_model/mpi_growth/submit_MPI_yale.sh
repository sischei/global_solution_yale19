#!/bin/bash

#SBATCH --job-name=mpi
#SBATCH --output=mpi_job.txt
#SBATCH --ntasks=4
#SBATCH --time=00:10:00

#make MPI available
module load Langs/Python
module load MPI/OpenMPI
module load Libs/NUMPY
module load Libs/MPI4PY

mpirun python main.py
