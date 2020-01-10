#!/usr/bin/env python3

import argparse, os, sys
from spectrum_functions import *
from pathlib import Path

pwd = os.getcwd()
parser = argparse.ArgumentParser(description='Generate the IR-spectrum for a molecule using GROMACS output files.')
parser.add_argument('-e'         , type=str  , required=True           , help='<Required> Experimental data directory')
parser.add_argument('-md'        , type=str  , required=True           , help='<Required> MD spectra directory')
parser.add_argument('-pool'      , type=int  , default=12              , help='Downsampling factor for MD spectra')
#parser.add_argument('-qmd'       , type=str  , required=True           , help='<Required> QM directory')
#parser.add_argument('-qms'	 , type=str  , required=True, nargs='+', help='<Required> List the names of the quantum mechanics directories')
parser.add_argument('-ffd'       , type=str  , required=True           , help='<Required> Input GROMACS directory with force field directories')
parser.add_argument('-ffs'       , type=str  , required=True, nargs='+', help='<Required> List the names of the force field directories found in the input GROMACS directory' )
parser.add_argument('-o'         , type=str  , default=pwd             , help='Output spectrum directory. DEFAULT: working directory')
parser.add_argument('-min'       , type=int  , default=500             , help='Lowest frequency to visualize in spectrum. DEFAULT: 0')
parser.add_argument('-max'       , type=int  , default=4000            , help='Highest frequency to visualize in spectrum. DEFAULT: 4000')
parser.add_argument('-n'         , type=int  , default=876             , help='Number of points on the frequency axis. DEFAULT: 1001')
parser.add_argument('-g'         , type=float, default=24              , help='Set gamma for Cauchy distribution. DEFAULT: 24')
#parser.add_argument('--linear'   ,                                       help='Molecule is linear'                    , action='store_true')
parser.add_argument('--png'      ,                                       help='Generate the spectrum as a PNG'        , action='store_true')
parser.add_argument('--pdf'      ,                                       help='Generate the spectrum as a PDF'        , action='store_true')
parser.add_argument('--svg'      ,                                       help='Generate the spectrum as a SVG'        , action='store_true')

if __name__ == "__main__":
  
	args = parser.parse_args()
	
	exp_dir      = args.e
	md_dir       = args.md
	n_pool       = args.pool
#	qm_dir       = args.qmd
#	qms          = args.qms
	ff_dir       = args.ffd
	ffs          = args.ffs
	output_dir   = args.o
#	linear       = args.linear
	start        = args.min
	stop         = args.max
	npoints      = args.n
	gamma        = args.g
	generate_png = args.png
	generate_pdf = args.pdf
	generate_svg = args.svg
	
	#molecules = find_molecules(exp_dir, qm_dir, qms, ff_dir, ffs)
	molecules = ['butane-23-dione', 'hexane-2-thiol', 'propane-13-diol', 'quinoline']
	argons = [0, 5, 20]
	#molecules = molecules[0:5]

	#print('\nThe following number of molecules were found in all listed directories and will be processed:', len(molecules))

	csv_dir = output_dir + "/CSV"
	if Path(csv_dir).is_dir() and os.listdir(csv_dir):
		print("the directory " + csv_dir + " already exists and it has contents. That directory will be emptied")
		os.system("rm -r " + csv_dir + "/*")
	elif Path(csv_dir).is_dir():
		print("the directory " + csv_dir + " already exists and but it has no contents. That directory will not be emptied")
	else:
		print("creating directory " + csv_dir)
		os.system("mkdir " + csv_dir)
	#stats_dir = csv_dir + "/SINGLE" 
	#os.system("mkdir " + stats_dir)
	
	types = argons + ffs	

	combis = ['molecule']
	oldtypes = []
	for type1 in types:
		oldtypes.append(type1)
		for type2 in types:
			if not type2 in oldtypes:
				combis.append(str(type1) + ' vs. ' + str(type2))
	statistics_file = csv_dir + '/cmp_statistics.csv'
	with open(statistics_file, 'w') as csvfile:
		writer = csv.writer(csvfile, delimiter='|')
		writer.writerow(combis)
		check_or_die(statistics_file, True)	
#
	all_spectra = {} 

	for molecule in molecules:
		print('\nNOW PROCESSING:', molecule)
		all_spectra[molecule] = cmp_spectra(exp_dir, md_dir, n_pool, ff_dir, ffs, molecule, argons, output_dir,
                                                     start, stop, npoints, gamma, generate_png, generate_pdf, generate_svg)

