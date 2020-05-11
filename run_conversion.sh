#!/bin/bash -l

#SBATCH -n 2
#SBATCH -t 40:00:00
#SBATCH -J running_conversion
#SBATCH -e conversion.e
#SBATCH -o conversion.o

# Load modules
module load miniconda/3
module load gromacs/2019

# Your commands
for dir in results_full results_medium results_short results_inverse
do

cd $dir
cd CSV

python /home/henning/projects/IR-Spectra/convert_crosscorrtables.py

cd ../..

done
