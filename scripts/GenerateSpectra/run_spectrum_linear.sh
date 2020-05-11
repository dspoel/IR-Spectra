#!/bin/bash -l

#SBATCH -n 2
#SBATCH -t 4:00:00
#SBATCH -J running_all_spectra
#SBATCH -e all_spectra.e
#SBATCH -o all_spectra.o

# Load modules
module load miniconda/3
module load gromacs/2019

# Your commands
python /home/henning/projects/IR-Spectra/scripts/GenerateSpectra/mk_spectra_linear.py \
-e /home/henning/JCAMP-DX_3 \
-ffd /home/henning/wd/THERMO \
-ffs CGenFF \
-o /home/henning/projects/IR-Spectra/results_linear \
--pdf

