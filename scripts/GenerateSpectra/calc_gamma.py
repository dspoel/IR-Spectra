#!/usr/bin/env python3

import argparse, os, sys
from spectrum_functions import *
from pathlib import Path
from scipy.optimize import minimize
from scipy.optimize import minimize_scalar
import numpy as np

pwd = os.getcwd()
parser = argparse.ArgumentParser(description='Generate the IR-spectrum for a molecule using GROMACS output files.')
parser.add_argument('-e'         , type=str  , required=True           , help='<Required> Experimental data directory')
parser.add_argument('-qmd'       , type=str  , required=True           , help='<Required> QM directory')
parser.add_argument('-qms'	 , type=str  , required=True, nargs='+', help='<Required> List the names of the quantum mechanics directories')
parser.add_argument('-ffd'       , type=str  , required=True           , help='<Required> Input GROMACS directory with force field directories')
parser.add_argument('-ffs'       , type=str  , required=True, nargs='+', help='<Required> List the names of the force field directories found in the input GROMACS directory' )
parser.add_argument('-o'         , type=str  , default=pwd             , help='Output spectrum directory. DEFAULT: working directory')
parser.add_argument('-min'       , type=int  , default=0               , help='Lowest frequency to visualize in spectrum. DEFAULT: 0')
parser.add_argument('-max'       , type=int  , default=4000            , help='Highest frequency to visualize in spectrum. DEFAULT: 4000')
parser.add_argument('-n'         , type=int  , default=4               , help='Number of points on the frequency axis. DEFAULT: 1001')
parser.add_argument('-g'         , type=float, default=24              , help='Set gamma for Cauchy distribution. DEFAULT: 24')
parser.add_argument('--linear'   ,                                       help='Molecule is linear'                    , action='store_true')
parser.add_argument('--png'      ,                                       help='Generate the spectrum as a PNG'        , action='store_true')
parser.add_argument('--pdf'      ,                                       help='Generate the spectrum as a PDF'        , action='store_true')
parser.add_argument('--svg'      ,                                       help='Generate the spectrum as a SVG'        , action='store_true')


def produce_measure(gamma):
	total_score = 0
	for molecule in molecules:
		spectra = []
		exp_spectrum, start, stop, deltax = read_exp_data(exp_dir, molecule)
		npoints = int(np.round(((stop - start) / deltax) + 1))
		spectra.append(exp_spectrum)	
		if method in qms:
			spectra.append(generate_spectrum(qm_dir, method, molecule, None, start, stop, npoints, gamma, scaling_factor))
		if method in ffs:
			spectra.append(generate_spectrum(ff_dir, method, molecule, eigfreq_count[molecule], start, stop, npoints, gamma, scaling_factor))
		pearson_score   = pearsonr(spectra[0][1], spectra[1][1])[0]
		#spearman_score  = spearmanr(spectra[0][1], spectra[1][1])[0]
		total_score += pearson_score
		#total_score += spearman_score
	return -1 * total_score / len(molecules)

if __name__ == "__main__":
  
	args = parser.parse_args()
	
	exp_dir      = args.e
	qm_dir       = args.qmd
	qms          = args.qms
	ff_dir       = args.ffd
	ffs          = args.ffs
	output_dir   = args.o
	linear       = args.linear
	start        = args.min
	stop         = args.max
	npoints      = args.n
	gamma        = args.g
	generate_png = args.png
	generate_pdf = args.pdf
	generate_svg = args.svg
	
	molecules = find_molecules(exp_dir, qm_dir, qms, ff_dir, ffs)
	
	print('\nThe following number of molecules were found in all listed directories and will be processed:', len(molecules))

	types = qms + ffs	

	exp_range_eigen_count = []
	all_scaling_factors = {}
	all_measures = {}
	unscaled_measures = {}
	#molecules = molecules[0:5]
	eigfreq_count = {}
	for molecule in molecules:
		spectrum = generate_spectrum(qm_dir, qms[0], molecule, None, start, stop, npoints, gamma, 1.0)
		eigfreq_count[molecule] = len(spectrum[2])

	for method in types:
		all_scaling_factors[method] = 1.0
	all_scaling_factors['G4'] = 0.965
	all_scaling_factors['OEP'] = 0.968

	for method in types: 
		scaling_factor = all_scaling_factors[method]
		gamma = minimize_scalar(produce_measure, bracket = (1.0, 50.0)).x
		all_measures[method] = -1 * produce_measure(gamma)
		unscaled_measures[method] = -1 * produce_measure(24.0)
		print(method, gamma, all_measures[method])
		print(method, "24.0" , unscaled_measures[method])
