#!/bin/bash -l

#SBATCH -n 2
#SBATCH -t 4:00:00
#SBATCH -J running_all_spectra
#SBATCH --mail-type=ALL
#SBATCH --mail-user alfred.andersson.9942@student.uu.se

# Load modules
module load miniconda/3
module load gromacs/2019

# Your commands
python /home/alfred/IR-Spectra/scripts/GenerateSpectra/mk_spectra.py \
-e /home/alfred/JCAMP-DX/ABSORBANCE \
-qmd /home/spoel/Liquids/MOLECULES \
-qms G4 OEP \
-ffd /home/spoel/wd/THERMO \
-ffs CGenFF GAFF-BCC GAFF-ESP \
-o /home/alfred/IR-Spectra/results \
--png
