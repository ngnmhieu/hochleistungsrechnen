#!/bin/bash

#SBATCH -N 1
#SBATCH -c 12
#SBATCH --ntasks-per-node 1

srun messung2

echo fertig > job.out