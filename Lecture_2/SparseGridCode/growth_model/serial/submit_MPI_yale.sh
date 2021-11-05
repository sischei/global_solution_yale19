#!/bin/bash

#SBATCH --job-name=SparseGrid
#SBATCH --output=SG_run.txt
#SBATCH --ntasks=1
#SBATCH --time=00:10:00

#make MPI available
module load Langs/Python
module load MPI/OpenMPI
module load Libs/NUMPY
module load Libs/MPI4PY

python main.py
