import os, math, re, csv, gzip
import numpy as np
import matplotlib.pyplot as plt
from numpy import linalg as la
from distutils.spawn import find_executable
from pathlib import Path
from spectrum_classes import *

def check_or_die(filenm, die):
	if not os.path.isfile(filenm):
		print("ERROR: generating " + filenm)
		if die:
			exit(1)

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

def run_one_nm(die, sigma, scale_factor, mdpdir, output):
	gmx = find_gmx()
	conf = "conf.gro"
	for minimize in [ "cg" ]:
		tpr  = minimize + ".tpr"
		os.system((gmx + " grompp -f " + mdpdir + "/%s.mdp -o %s -v -maxwarn 1 -c conf.gro") % (minimize, tpr))
		check_or_die(tpr, die)
		conf = ("after_%s.g96" % minimize )
		os.system((gmx + " mdrun -s %s -c %s -v") % (tpr, conf))
		check_or_die(conf, die)
	tpr = "nm.tpr"
	os.system((gmx + " grompp -c %s -f " + mdpdir + "/nm.mpd -o %s -v") % (conf, tpr))
	check_or_die(tpr, die)
	mtz = "nm.mtx"
	os.system(gmx + " mdrun -s %s -v -mtx %s" % ( tpr, mtz ))
	check_or_die(mtz, die )
	os.system(("%s nmeig -last 1000 -s %s -f %s -sigma %d > %s") % (gmx, tpr, mtz, sigma, output))

def extract_eigenfrequencies(molecule_path):
	"""Extract eigenfrequencies from GROMACS eigenfrequency file"""
	eigenfrequencies = []
	for line in open(molecule_path + "/eigenfreq.xvg", "r").readlines():
		if not line == line.lstrip():
			words = line.strip().split()
			eigenfrequencies.append(float(words[1]))
	return eigenfrequencies

def extract_eigenvectors(molecule_path):
	"""Extract eigenvectors from GROMACS eigenvector file"""
	eigenvectors = []       
	eigenvector  = np.empty([0,3])
	os.system(find_gmx() + " dump -f " + molecule_path + "/eigenvec.trr > " + molecule_path + "/eigenvec.txt")
	for line in open(molecule_path + "/eigenvec.txt", "r").readlines():
		if re.match(r"\s+x\[", line):
			values = re.split("[{},]",line)
			eigenvector = np.append(eigenvector, [[float(values[1]),float(values[2]),float(values[3])]], axis=0)
		elif len(eigenvector) > 0:
			eigenvectors.append(eigenvector)
			eigenvector = np.empty([0,3])
	eigenvectors.append(eigenvector)
	os.remove(molecule_path + "/eigenvec.txt")
	return eigenvectors[1:]

def extract_atomic_properties(molecule_path):
	"""Extract atomic properties from GROMACS topology file"""
	squared_masses = []
	charge_mass_factors = []
	for line in open(molecule_path + "/topol.top", "r").readlines():
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
	for normal_mode in molecule.normal_modes():
		normal_mode.calculate_intensity(molecule.atoms())
	return molecule

def generate_cauchy_distribution(frequencies, eigenfrequency, gamma, intensity):
	"""Generate a Cauchy distribution"""
	cauchy = np.zeros(len(frequencies))
	for i in range(len(frequencies)):
		cauchy[i] = intensity*(1/math.pi)*(((1/2)*gamma)/((frequencies[i]-eigenfrequency)**2+((1/2)*gamma)**2))
	return cauchy

def generate_spectrum(molecule_path, molecule_name, linear, log, start, stop, step_size, gamma):
	if log:
		return generate_spectrum_from_log(molecule_path, molecule_name, start, stop, step_size, gamma)
	else:
		frequencies         = np.linspace(start, stop, int((stop-start)/step_size)+1)
		intensity_all_modes = np.zeros(len(frequencies))
		molecule            = generate_molecule(molecule_path, linear)
		normal_modes        = molecule.normal_modes()
		for normal_mode in normal_modes:
			intensity_one_mode   = generate_cauchy_distribution(frequencies, normal_mode.eigenfrequency(), gamma, normal_mode.intensity())
			intensity_all_modes += intensity_one_mode
		return [frequencies, intensity_all_modes, molecule_name]

def save_spectrum_as_csv(spectrum, outdir, outname, header):
	"""Write the spectrum of all normal modes of a molecule to a CSV-file"""
	outformat = ".csv"
	output    = outdir + "/" + outname + outformat
	with open(output, 'w') as csvfile:
		writer = csv.writer(csvfile, delimiter=',')
		if header:
			writer.writerow(['frequency', 'intensity'])
		for i, frequency in enumerate(spectrum[0]):
			writer.writerow([frequency, spectrum[1][i]])
	check_or_die(output, False)

def save_spectrum_as_figure(spectra, outdir, outname, outformat):
	"""Write the spectrum of all normal modes of a molecule as a PNG, PDF or SVG"""
	output = outdir + "/" + outname + '.' + outformat
	plt.figure(figsize=(18, 8))
	for spectrum in spectra:
		plt.plot(spectrum[0], spectrum[1]/max(spectrum[1]), label='Spectrum for ' + spectrum[2])
	plt.legend(loc='upper right')
	plt.xlabel('Frequency, $cm^{-1}$')
	plt.ylabel('IR intensity')		
	plt.yticks([])
	plt.rcParams.update({'font.size': 18})
	plt.savefig(output, format=outformat)
	check_or_die(output, False)

def save_spectrum(molecule_path1, molecule_name1, linear1, log1, 
                  molecule_path2, molecule_name2, linear2, log2, 
                  start, stop, step_size, gamma, outdir, outname, 
                  header, csv, png, pdf, svg):

	possible_formats         = ["csv", "png", "pdf", "svg"]
	desired_formats          = [ csv,   png,   pdf,   svg ]
	spectra = [generate_spectrum(molecule_path1, molecule_name1, linear1, log1, start, stop, step_size, gamma)]
	if molecule_path2:
		spectra.append(generate_spectrum(molecule_path2, molecule_name2, linear2, log2, start, stop, step_size, gamma))
	for i, desired_format in enumerate(desired_formats):
		if desired_format:
			selected_format = possible_formats[i]
			if Path(outdir + "/" + outname + '.' + selected_format).is_file():
				raise Exception('desired output file name "' + outname + "." + selected_format + '" is already in use. Please (re)move the old file or pick another output file name')
			elif not Path(outdir).is_dir():
				raise Exception('the specified output folder  "' + outdir + '" does not exist. Please create it or pick another output folder that exists')
			elif selected_format in [ "png", "pdf", "svg" ]:
				save_spectrum_as_figure(spectra, outdir, outname, selected_format)
			elif selected_format in [ "csv" ]:
				save_spectrum_as_csv(spectra, outdir, outname, header)

def generate_spectrum_from_log(log, molecule_name, start, stop, step_size, gamma):
	if os.path.exists(log):
		eigenfrequencies = []
		intensities      = []
		if log.endswith(".gz"):
			lines = gzip.open(log, 'rt').readlines()
		else:
			lines = open(log, 'r').readlines()
		for line in lines:
			if "Frequencies" in line:
				words  = re.findall(r"[-+]?\d*\.\d+|\d+", line)
				for word in words:
					eigenfrequencies.append(float(word))
			elif "IR Inten" in line:
				words  = re.findall(r"[-+]?\d*\.\d+|\d+", line)
				for word in words:
					intensities.append(float(word))
	else:
        	sys.exit("The QM log file does not exist!!!")
	if eigenfrequencies and intensities:
		frequencies         = np.linspace(start, stop, int((stop-start)/step_size)+1)
		intensity_all_modes = np.zeros(len(frequencies))
		for i, eigenfrequency in enumerate(eigenfrequencies):
			intensity_one_mode   = generate_cauchy_distribution(frequencies, eigenfrequency, gamma, intensities[i])
			intensity_all_modes += intensity_one_mode
		return [frequencies, intensity_all_modes, molecule_name]
	else:
		sys.exit("There are no frequency and intensty in the QM log file!!!")
