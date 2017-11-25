#!/bin/bash

#SBATCH -N 3
#SBATCH -c 12
#SBATCH --ntasks-per-node 1

srun messung

echo fertig > job.out
