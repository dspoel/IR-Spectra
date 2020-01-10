#!/bin/bash -l

#SBATCH -n 2
#SBATCH -t 4:00:00
#SBATCH -J running_all_spectra
#SBATCH -e all_spectra_eic.e
#SBATCH -o all_spectra_eic.o

# Load modules
module load miniconda/3
module load gromacs/2019

# Your commands
python /home/henning/projects/IR-Spectra/scripts/GenerateSpectra/cmp_spectra.py \
-e /home/henning/JCAMP-DX \
-md /home/henning/projects/IR-Spectra/spectra_from_MD \
-pool 12 \
-ffd /home/henning/wd/THERMO \
-ffs CGenFF \
-o /home/henning/projects/IR-Spectra/results_comp_MD \
-max 4000 \
-g 24 \
--png

