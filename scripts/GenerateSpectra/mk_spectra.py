#!/usr/bin/env python3

import argparse, os, sys
from spectrum_functions import *
from pathlib import Path

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
parser.add_argument('-n'         , type=int  , default=1001            , help='Number of points on the frequency axis. DEFAULT: 1001')
parser.add_argument('-g'         , type=float, default=24              , help='Set gamma for Cauchy distribution. DEFAULT: 24')
parser.add_argument('--cross'    ,                                       help='Run cross comparison (requires all experimental spectra to have the same resolution)'    , action='store_true')
parser.add_argument('--linear'   ,                                       help='Molecule is linear'                    , action='store_true')
parser.add_argument('--png'      ,                                       help='Generate the spectrum as a PNG'        , action='store_true')
parser.add_argument('--pdf'      ,                                       help='Generate the spectrum as a PDF'        , action='store_true')
parser.add_argument('--svg'      ,                                       help='Generate the spectrum as a SVG'        , action='store_true')

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
	do_cross     = args.cross
	generate_png = args.png
	generate_pdf = args.pdf
	generate_svg = args.svg
	
	molecules = find_molecules(exp_dir, qm_dir, qms, ff_dir, ffs)
	#molecules = ['cis-decalin', '13-oxazole', 'ethyl-sulfate', 'diethyl-oxalate']
	#molecules = molecules[0:5]

	print('\nThe following number of molecules were found in all listed directories and will be processed:', len(molecules))

	csv_dir = output_dir + "/CSV"
	if Path(csv_dir).is_dir() and os.listdir(csv_dir):
		print("the directory " + csv_dir + " already exists and it has contents. That directory will be emptied")
		os.system("rm -r " + csv_dir + "/*")
	elif Path(csv_dir).is_dir():
		print("the directory " + csv_dir + " already exists and but it has no contents. That directory will not be emptied")
	else:
		print("creating directory " + csv_dir)
		os.system("mkdir " + csv_dir)
	stats_dir = csv_dir + "/SINGLE" 
	os.system("mkdir " + stats_dir)
	
	types = qms + ffs	

	for type in types:
		statistics_file = stats_dir + "/" + type + '_statistics.csv'
		with open(statistics_file, 'w') as csvfile:
			writer = csv.writer(csvfile, delimiter=',')
			writer.writerow(['molecule', 'cos', 'pearson', 'spearman'])
		check_or_die(statistics_file, True)	

	all_spectra = {} 

	for molecule in molecules:
		print('\nNOW PROCESSING:', molecule)
		all_spectra[molecule] = save_spectrum(exp_dir, qm_dir, qms, ff_dir, ffs, molecule, output_dir,
                                                     start, stop, npoints, gamma, generate_png, generate_pdf, generate_svg)

	if do_cross:
		cross_compare(molecules, all_spectra, types, output_dir)
		exp_intracorr(molecules, all_spectra, output_dir)
