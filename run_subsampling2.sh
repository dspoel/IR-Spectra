#!/bin/bash -l

#SBATCH -n 2
#SBATCH -t 48:00:00
#SBATCH -J running_subsample
#SBATCH -e subs2.e
#SBATCH -o subs2.o

# Load modules
module load miniconda/3
module load gromacs/2019

# Your commands
python subsample2.py
