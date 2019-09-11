#!/bin/bash -l

#SBATCH -n 2
#SBATCH -t 4:00:00
#SBATCH -J running_all_spectra
#SBATCH --mail-type=ALL
#SBATCH -e all_spectra.e
#SBATCH -o all_spectra.o
#SBATCH --mail-user henning.henschel@icm.uu.se

# Load modules
module load miniconda/3
module load gromacs/2019

# Your commands
python /home/henning/projects/IR-Spectra/scripts/GenerateSpectra/mk_spectra.py \
-e /home/henning/JCAMP-DX \
-qmd /home/henning/Liquids/MOLECULES \
-qms G4 OEP \
-ffd /home/henning/wd/THERMO \
-ffs CGenFF GAFF-BCC GAFF-ESP \
-o /home/henning/projects/IR-Spectra/scaled_results \
--png
