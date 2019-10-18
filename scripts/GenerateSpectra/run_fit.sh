#!/bin/bash -l

#SBATCH -n 2
#SBATCH -t 6:00:00
#SBATCH -J running_all_spectra
#SBATCH -e fit_spectra.e 
#SBATCH -o fit_spectra.o

# Load modules
module load miniconda/3
module load gromacs/2019

# Your commands
python /home/henning/projects/IR-Spectra/scripts/GenerateSpectra/calc_scale_factor.py \
-e /home/henning/JCAMP-DX \
-qmd /home/henning/Liquids/MOLECULES \
-qms G4 OEP \
-ffd /home/henning/wd/THERMO \
-ffs CGenFF GAFF-BCC GAFF-ESP OPLS \
-o /home/henning/projects/IR-Spectra/results_with_opls \
--png
