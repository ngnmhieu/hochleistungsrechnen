#!/bin/bash

#SBATCH -N 3
#SBATCH -c 12
#SBATCH -p west
#SBATCH --ntasks-per-node 1

srun messung-seq.sh

echo "fertig" > job.out
