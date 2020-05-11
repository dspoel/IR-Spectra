#!/bin/bash -l

#SBATCH -n 2
#SBATCH -t 10:00:00
#SBATCH -J running_all_spectra
#SBATCH -e all_spectra_eic.e
#SBATCH -o all_spectra_eic.o

# Load modules
module load miniconda/3
module load gromacs/2019

# Your commands
python /home/henning/projects/IR-Spectra/scripts/GenerateSpectra/mk_spectra.py \
-e /home/henning/JCAMP-DX_4 \
-qmd /home/henning/Liquids/MOLECULES \
-qms G4 OEP \
-ffd /home/henning/wd/THERMO \
-ffs CGenFF GAFF-BCC GAFF-ESP OPLS \
-o /home/henning/projects/IR-Spectra/results_full \
-min 550 \
-max 3846 \
--pdf

python /home/henning/projects/IR-Spectra/scripts/GenerateSpectra/mk_spectra.py \
-e /home/henning/JCAMP-DX_4 \
-qmd /home/henning/Liquids/MOLECULES \
-qms G4 OEP \
-ffd /home/henning/wd/THERMO \
-ffs CGenFF GAFF-BCC GAFF-ESP OPLS \
-o /home/henning/projects/IR-Spectra/results_medium \
-min 550 \
-max 2000 \
--pdf

python /home/henning/projects/IR-Spectra/scripts/GenerateSpectra/mk_spectra.py \
-e /home/henning/JCAMP-DX_4 \
-qmd /home/henning/Liquids/MOLECULES \
-qms G4 OEP \
-ffd /home/henning/wd/THERMO \
-ffs CGenFF GAFF-BCC GAFF-ESP OPLS \
-o /home/henning/projects/IR-Spectra/results_short \
-min 550 \
-max 1650 \
--pdf

python /home/henning/projects/IR-Spectra/scripts/GenerateSpectra/mk_spectra.py \
-e /home/henning/JCAMP-DX_4 \
-qmd /home/henning/Liquids/MOLECULES \
-qms G4 OEP \
-ffd /home/henning/wd/THERMO \
-ffs CGenFF GAFF-BCC GAFF-ESP OPLS \
-o /home/henning/projects/IR-Spectra/results_inverse \
-min 2000 \
-max 3846 \
--pdf

