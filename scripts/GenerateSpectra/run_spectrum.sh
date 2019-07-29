#!/bin/bash -l

#SBATCH -n 2
#SBATCH -t 4:00:00
#SBATCH -J running_all_spectra
#SBATCH --mail-type=ALL
#SBATCH --mail-user alfred.andersson.9942@student.uu.se

# Load modules
module load python/3.6.2
module load miniconda/3
module load gromacs/2019

# Your commands
python3.6 mk_spectra.py \
-g4 /home/spoel/Liquids/MOLECULES/G4 \
-ffd /home/spoel/wd/THERMO \
-ffs GAFF-ESP CGenFF GAFF-BCC \
-o /home/alfred/IR-Spectra/scripts/GenerateSpectra/results \
--png
