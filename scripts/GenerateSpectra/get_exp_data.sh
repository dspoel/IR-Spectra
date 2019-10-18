#!/bin/bash -l

#SBATCH -n 2
#SBATCH -t 5:00:00
#SBATCH -J get_spectra
#SBATCH -o get_spectra.o
#SBATCH -e get_spectra.e
#SBATCH --mail-type=ALL
#SBATCH --mail-user henning.henschel@icm.uu.se

module load miniconda/3
export PYTHONPATH="/home/henning/Liquids/PYTHON/"
export LIQUIDS="/home/henning/Liquids/"

python3 /home/henning/projects/IR-Spectra/scripts/GenerateSpectra/get_exp_data.py
