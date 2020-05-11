#!/bin/bash -l

#SBATCH -n 2
#SBATCH -t 48:00:00
#SBATCH -J running_subsample
#SBATCH -e subs.e
#SBATCH -o subs.o

# Load modules
module load miniconda/3
module load gromacs/2019

# Your commands
python subsample.py
