import argparse, os
from spectrum_functions import *

pwd = os.getcwd()

parser = argparse.ArgumentParser(description='Generate the IR-spectrum for a molecule using GROMACS output files.')
parser.add_argument('-i1'      , type=str  , default=pwd         , help='First input molecule directory with conf.gro and topol.top. DEFAULT: working directory')
parser.add_argument('-i2'      , type=str  , default=None        , help='Second input molecule directory. DEFAULT: None')
parser.add_argument('-mnm1'    , type=str  , default="Molecule 1", help='Name of first molecule. DEFAULT: Molecule 1')
parser.add_argument('-mnm2'    , type=str  , default="Molecule 2", help='Name of second molecule. DEFAULT: Molecule 2')
parser.add_argument('-mdp'     , type=str  , default=pwd         , help='Input mdp directory with cg.mdp and nm.mdp. DEFAULT: working directory')
parser.add_argument('-o'       , type=str  , default="spectrum"  , help='Output spectrum name. DEFAULT: spectrum')
parser.add_argument('-od'      , type=str  , default=pwd         , help='Output spectrum directory. DEFAULT: working directory')
parser.add_argument('-vmin'    , type=int  , default=0           , help='Lowest frequency to visualize in spectrum. DEFAULT: 0')
parser.add_argument('-vmax'    , type=int  , default=4000        , help='Highest frequency to visualize in spectrum. DEFAULT: 4000')
parser.add_argument('-s'       , type=float, default=4           , help='Step size on the frequency axis. DEFAULT: 4')
parser.add_argument('-g'       , type=float, default=24          , help='Set gamma for Cauchy distribution. DEFAULT: 24')
parser.add_argument('--log1'   , 	                           help='Read first molecule from log-file'     , action='store_true')
parser.add_argument('--log2'   , 	                           help='Read second molecule from log-file'    , action='store_true')
parser.add_argument('--no_nm'  ,                                   help='Do not run nm'                         , action='store_true')
parser.add_argument('--linear1',                                   help='Molecule is linear'                    , action='store_true')
parser.add_argument('--linear2',                                   help='Molecule is linear'                    , action='store_true')
parser.add_argument('--header' ,                                   help='CSV file is generated without a header', action='store_false')
parser.add_argument('--csv'    ,                                   help='Generate the spectrum as a CSV'        , action='store_true')
parser.add_argument('--png'    ,                                   help='Generate the spectrum as a PNG'        , action='store_true')
parser.add_argument('--pdf'    ,                                   help='Generate the spectrum as a PDF'        , action='store_true')
parser.add_argument('--svg'    ,                                   help='Generate the spectrum as a SVG'        , action='store_true')

args = parser.parse_args()

molecule_path1 = args.i1
molecule_path2 = args.i2
molecule_name1 = args.mnm1
molecule_name2 = args.mnm2
mdpdir         = args.mdp
outname        = args.o
outdir         = args.od
no_nm          = args.no_nm
start          = args.vmin
stop           = args.vmax
step_size      = args.s
gamma          = args.g
log1           = args.log1
log2           = args.log2
linear1        = args.linear1
linear2        = args.linear2
header         = args.header
generate_csv   = args.csv
generate_png   = args.png
generate_pdf   = args.pdf
generate_svg   = args.svg

if not no_nm:
	run_one_nm(True, 1, 1, mdpdir, "nm_output")

save_spectrum(molecule_path1, molecule_name1, linear1, log1, molecule_path2, molecule_name2, linear2, log2, start, stop, step_size, gamma, outdir, outname, header, generate_csv, generate_png, generate_pdf, generate_svg)

