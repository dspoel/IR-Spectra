import argparse, os
from spectrum_functions import *

pwd = os.getcwd() + '/'

parser = argparse.ArgumentParser(description='Generate the IR-spectrum for a molecule using GROMACS output files.')
parser.add_argument('-i'   , '--input'  , type=str  , default=pwd       , help='input molecule directory')
parser.add_argument('-o'   , '--outnm'  , type=str  , default="spectrum", help='output spectrum name')
parser.add_argument('-odir', '--outdir' , type=str  , default=pwd       , help='output spectrum directory')
parser.add_argument('-l'   , '--linear' , type=bool , default=False     , help='states whether molecule is linear or not')
parser.add_argument('-vmin', '--start'  , type=int  , default=0         , help='lowest frequency to visualize in spectrum')
parser.add_argument('-vmax', '--stop'   , type=int  , default=4000      , help='highest frequency to visualize in spectrum')
parser.add_argument('-ss'  , '--step_sz', type=float, default=4         , help='step size (resoulution) on the frequency axis')
parser.add_argument('-g'   , '--gamma'  , type=float, default=24        , help='input molecule directory')
parser.add_argument('-hdr' , '--header' , type=bool , default=True      , help='csv file is generated with a header')
parser.add_argument('-csv' , '--csv'    , type=bool , default=False     , help='if true, a csv file will be generated')
parser.add_argument('-png' , '--png'    , type=bool , default=False     , help='if true, a png file will be generated')
parser.add_argument('-pdf' , '--pdf'    , type=bool , default=False     , help='if true, a pdf file will be generated')
parser.add_argument('-svg' , '--svg'    , type=bool , default=False     , help='if true, a svg file will be generated')

args = parser.parse_args()

molecule_path = args.input
linear        = args.linear
start         = args.start
stop          = args.stop
step_size     = args.step_sz
gamma         = args.gamma
outdir        = args.outdir
outname       = args.outnm
header        = args.header
generate_csv  = args.csv
generate_png  = args.png
generate_pdf  = args.pdf
generate_svg  = args.svg

save_spectrum(molecule_path, linear, start, stop, step_size, gamma, outdir, outname, header, generate_csv, generate_png, generate_pdf, generate_svg)

