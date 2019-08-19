#!/bin/bash -l

#SBATCH -n 2
#SBATCH -t 5:00:00
#SBATCH -J running_all_spectra
#SBATCH --mail-type=ALL
#SBATCH --mail-user alfred.andersson.9942@student.uu.se

module load miniconda/3

python /home/alfred/IR-Spectra/scripts/GenerateSpectra/get_exp_data.py
