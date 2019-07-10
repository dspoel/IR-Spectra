import os, math, re, csv
import numpy as np
import matplotlib.pyplot as plt
from numpy import linalg as la
from distutils.spawn import find_executable
from pathlib import Path
from spectrum_classes import *

def find_gmx():
	"""Find where GROMACS is installed"""
	gmx = None
	for mpi in [ "_mpi", "" ]:
		for double in [ "_d", ""]:
			gmx = find_executable("gmx" + mpi + double)
			if gmx:
				return gmx
	if not gmx:
		sys.exit("GROMACS is not installed!!!") 

def extract_eigenfrequencies(molecule_path):
	"""Extract eigenfrequencies from GROMACS eigenfrequency file"""
	eigenfrequencies = []
	for line in open(molecule_path + "eigenfreq.xvg", "r").readlines():
		if not line == line.lstrip():
			words = line.strip().split()
			eigenfrequencies.append(float(words[1]))
	return eigenfrequencies

def extract_eigenvectors(molecule_path):
	"""Extract eigenvectors from GROMACS eigenvector file"""
	eigenvectors = []       
	eigenvector  = np.empty([0,3])
	os.system(find_gmx() + " dump -f " + molecule_path + "eigenvec.trr > " + molecule_path + "eigenvec.txt")
	for line in open(molecule_path + "eigenvec.txt", "r").readlines():
		if re.match(r"\s+x\[", line):
			values = re.split("[{},]",line)
			eigenvector = np.append(eigenvector, [[float(values[1]),float(values[2]),float(values[3])]], axis=0)
		elif len(eigenvector) > 0:
			eigenvectors.append(eigenvector)
			eigenvector = np.empty([0,3])
	eigenvectors.append(eigenvector)
	os.remove(molecule_path + "eigenvec.txt")
	return eigenvectors[1:]

def extract_atomic_properties(molecule_path):
	"""Extract atomic properties from GROMACS topology file"""
	squared_masses = []
	charge_mass_factors = []
	for line in open(molecule_path + "topol.top", "r").readlines():
		words = line.strip().split()
		if len(words) == 11:
			squared_masses.append(math.sqrt(float(words[7])))
			charge_mass_factors.append(float(words[6])/math.sqrt(float(words[7])))            
	return squared_masses, charge_mass_factors

def generate_atoms(molecule_path):
	"""Initialize Atom objects and return them as a list"""
	squared_masses, charge_mass_ratios = extract_atomic_properties(molecule_path)
	atoms = []
	for i in range(len(squared_masses)):
		atom = Atom(squared_masses[i], charge_mass_ratios[i])
		atoms.append(atom)
	return atoms

def generate_normal_modes(molecule_path):
	"""Initialize NormalMode objects and return them as a list"""
	eigenfrequencies = extract_eigenfrequencies(molecule_path)
	eigenvectors     = extract_eigenvectors(molecule_path)
	normal_modes     = []
	for i in range(len(eigenfrequencies)):
		normal_mode = NormalMode(eigenfrequencies[i],eigenvectors[i])
		normal_modes.append(normal_mode)
	return normal_modes

def generate_molecule(molecule_path, linear):
	"""Initialize Molecule object"""
	atoms        = generate_atoms(molecule_path)
	normal_modes = generate_normal_modes(molecule_path)
	molecule     = Molecule(linear, atoms, normal_modes)
	for normal_mode in molecule.normal_modes:
		normal_mode.calculate_intensity(molecule.atoms)
	return molecule

def generate_cauchy_distribution(frequencies, eigenfrequency, gamma, intensity):
	"""Generate a Cauchy distribution"""
	cauchy = np.zeros(len(frequencies))
	for i in range(len(frequencies)):
		cauchy[i] = intensity*(1/math.pi)*(((1/2)*gamma)/((frequencies[i]-eigenfrequency)**2+((1/2)*gamma)**2))
	return cauchy

def generate_spectrum(molecule_path, linear, start, stop, step_size, gamma):
	frequencies = np.linspace(start, stop, int((stop-start)/step_size)+1)
	spectrum_all = np.zeros(len(frequencies))
	molecule = generate_molecule(molecule_path, linear)
	normal_modes = molecule.normal_modes
	for normal_mode in normal_modes:
		spectrum = generate_cauchy_distribution(frequencies, normal_mode.eigenfrequency, gamma, normal_mode.intensity)
		spectrum_all += spectrum
	return frequencies, spectrum_all

def save_spectrum_as_csv(frequencies, intensities, outdir, outname, header):
	"""Write the spectrum of all normal modes of a molecule to a CSV-file"""
	outformat = ".csv"
	output    = outdir + outname + outformat
	with open(output, 'w') as csvfile:
		writer = csv.writer(csvfile, delimiter=',')
		if header:
			writer.writerow(['frequency', 'intensity'])
		for i, frequency in enumerate(frequencies):
			writer.writerow([frequency, intensities[i]])

def save_spectrum_as_figure(frequencies, intensities, outdir, outname, outformat):
	"""Write the spectrum of all normal modes of a molecule as a PNG, PDF or SVG"""
	output = outdir + outname + '.' + outformat
	plt.figure(figsize=(18, 8))
	plt.plot(frequencies, intensities, label='Spectrum for all non-vibrational normal modes')
	plt.legend(loc='upper right')
	plt.xlabel('Frequency, $cm^{-1}$')
	plt.ylabel('IR intensity')		
	plt.yticks([])
	plt.rcParams.update({'font.size': 18})
	plt.savefig(output, format=outformat)

def save_spectrum(molecule_path, linear, start, stop, step_size, gamma, outdir, outname, header, csv, png, pdf, svg):
	possible_formats         = ["csv", "png", "pdf", "svg"]
	desired_formats          = [ csv,   png,   pdf,   svg ]
	frequencies, intensities = generate_spectrum(molecule_path, linear, start, stop, step_size, gamma)
	for i, desired_format in enumerate(desired_formats):
		if desired_format:
			selected_format = possible_formats[i]
			if Path(outdir + outname + '.' + selected_format).is_file():
				raise Exception('desired output file name "' + outname + "." + selected_format + '" is already in use. Please (re)move the old file or pick another output file name')
			elif not Path(outdir).is_dir():
				raise Exception('the specified output folder  "' + outdir + '" does not exist. Please create it or pick another output folder that exists')
			elif selected_format in [ "png", "pdf", "svg" ]:
				save_spectrum_as_figure(frequencies, intensities, outdir, outname, selected_format)
			elif selected_format in [ "csv" ]:
				save_spectrum_as_csv(frequencies, intensities, outdir, outname, header)
