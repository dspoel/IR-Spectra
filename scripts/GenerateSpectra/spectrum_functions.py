import os, sys, math, re, csv, glob, gzip, itertools
import numpy as np
import matplotlib.pyplot as plt
plt.switch_backend('agg')
plt.rcParams.update({'font.size': 11})
from numpy import linalg as la
from numpy import trapz
from distutils.spawn import find_executable
from pathlib import Path
from spectrum_classes import *

def intersection(lst1, lst2): 
	return set(lst1).intersection(lst2)

def find_molecules(exp_dir, qm_dir, qms, ff_dir, ffs):
	exp_molecules    = [os.path.splitext(os.path.basename(f.path))[0] for f in os.scandir(exp_dir) if f.is_file()]
	print("Number of molecules in %s data directory: %d" % ("experimental", len(exp_molecules)))
	common_molecules = exp_molecules
	for qm in qms:
		qm_molecules     = [f.path.split("/")[-1] for f in os.scandir(qm_dir + "/" + qm) if f.is_dir()]
		print("Number of molecules in %s data directory: %d" % (qm, len(qm_molecules)))
		common_molecules = intersection(common_molecules, qm_molecules)
	for ff in ffs:
		ff_molecules     = [f.path.split("/")[-1] for f in os.scandir(ff_dir + "/" + ff) if f.is_dir()]
		print("Number of molecules in %s data directory: %d" % (ff, len(ff_molecules)))
		common_molecules = intersection(common_molecules, ff_molecules)
	return common_molecules

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
	full_path = molecule_path + "/eigenfreq.xvg"
	print('reading eigenfrequencies from:', full_path)
	for line in open(full_path, "r").readlines():
		if not line == line.lstrip():
			words = line.strip().split()
			eigenfrequencies.append(float(words[1]))
	return eigenfrequencies

def extract_eigenvectors(molecule_path):
	"""Extract eigenvectors from GROMACS eigenvector file"""
	eigenvectors = []
	eigenvector  = np.empty([0,3])
	full_path    = molecule_path + "/eigenvec.trr"
	temp_txt     = "eigenvec.txt"
	print('reading eigenvectors from:', full_path )
	os.system(find_gmx() + " dump -f " + full_path + " -quiet > " + temp_txt)
	for line in open(temp_txt, "r").readlines():
		if re.match(r"\s+x\[", line):
			values = re.split("[{},]",line)
			eigenvector = np.append(eigenvector, [[float(values[1]),float(values[2]),float(values[3])]], axis=0)
		elif len(eigenvector) > 0:
			eigenvectors.append(eigenvector)
			eigenvector = np.empty([0,3])
	eigenvectors.append(eigenvector)
	os.remove(temp_txt)
	return eigenvectors[1:]

def read_topol(topol):
	squared_masses = []
	charges        = []
	for line in open(topol, "r").readlines():
		words = line.strip().split()
		if len(words) == 11:
			squared_masses.append([math.sqrt(float(words[7]))])
			charges.append([float(words[6])]*3)
	return np.array(squared_masses), np.array(charges)

def read_mu(mu):
	charges = []
	row     = []
	counter = 1
	for line in open(mu, "r").readlines():
		values = line.split()
		row.append(float(values[counter]))
		if counter == 3:
			charges.append(row)
			row     = []
			counter = 1
		else:
			counter += 1
	return np.array(charges)

def extract_atomic_properties(molecule_path):
	"""Extract atomic properties from GROMACS topology file"""
	topol     = molecule_path + "/topol.top"
	mu        = molecule_path + "/mu.txt"
	mu_exists = os.path.isfile(mu)
	msg       = 'reading atomic properties from: %s' % (topol)
	if mu_exists:
		msg = "%s and %s" % (msg, mu)
	print(msg)
	squared_masses, charges = read_topol(topol)
	if mu_exists:
		charges = read_mu(mu)
		squared_masses = squared_masses[:charges.shape[0],:]            
	charge_mass_factors = np.zeros(charges.shape)
	for i in range(charges.shape[0]):
		for j in range(3):
			charge_mass_factors[i,j] = charges[i,j]/squared_masses[i]
	return squared_masses, charge_mass_factors

def generate_atoms(molecule_path):
	"""Initialize Atom objects and return them as a list"""
	squared_masses, charge_mass_ratios = extract_atomic_properties(molecule_path)
	atoms = []
	for i in range(squared_masses.shape[0]):
		atom = Atom(squared_masses[i], charge_mass_ratios[i,:])
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

def generate_molecule(molecule_path, eigfreq_count):
	"""Initialize Molecule object"""
	atoms        = generate_atoms(molecule_path)
	normal_modes = generate_normal_modes(molecule_path)
	molecule     = Molecule(eigfreq_count, atoms, normal_modes)
	for normal_mode in molecule.normal_modes():
		normal_mode.calculate_intensity(molecule.atoms())
	return molecule

def generate_cauchy_distribution(frequencies, eigenfrequency, gamma, intensity):
	"""Generate a Cauchy distribution"""
	cauchy = np.zeros(len(frequencies))
	for i in range(len(frequencies)):
		cauchy[i] = intensity*(1/math.pi)*(((1/2)*gamma)/((frequencies[i]-eigenfrequency)**2+((1/2)*gamma)**2))
	return cauchy

def generate_spectrum(input_dir, origin, molecule, eigfreq_count, start, stop, npoints, gamma):
	print("\n<" + origin + ">")
	if origin in ["G4", "OEP"]:
		path = "%s/%s/%s" % (input_dir, origin, molecule)
		return generate_spectrum_from_log(path, origin, start, stop, npoints, gamma)
	elif origin in ["CGenFF", "GAFF-ESP", "GAFF-BCC"]:
		frequencies         = np.linspace(start, stop, npoints)
		intensity_all_modes = np.zeros(len(frequencies))
		path = input_dir + "/" + origin + "/" + molecule
		molecule            = generate_molecule(path, eigfreq_count)
		normal_modes        = molecule.normal_modes()
		eigenfrequencies    = []
		for normal_mode in normal_modes:
			intensity_one_mode   = generate_cauchy_distribution(frequencies, normal_mode.eigenfrequency(), gamma, normal_mode.intensity())
			intensity_all_modes += intensity_one_mode
			eigenfrequencies.append(normal_mode.eigenfrequency())
		return [frequencies, intensity_all_modes, np.array(eigenfrequencies), origin]
	else:
		raise Exception("%s is not a supported format. Please add it to code!" % (origin))

def normalize_spectra(spectra):
	for i, spectrum in enumerate(spectra):
		spectra[i][1] = spectrum[1]/np.trapz(spectrum[1],spectrum[0])
	return spectra

def cosine_distance(spectrum1, spectrum2):
	return (spectrum1.dot(spectrum2))/(math.sqrt(spectrum1.dot(spectrum1))*math.sqrt(spectrum2.dot(spectrum2)))

def rmsd(eigfreqs1, eigfreqs2):
	return math.sqrt(sum((eigfreqs1-eigfreqs2)**2)/len(eigfreqs1))

def save_spectra_as_figure(spectra, output_dir, molecule, outformat):
	"""Write the spectrum of all normal modes of a molecule as a PNG, PDF or SVG"""
	outformat_dir = output_dir + "/" + outformat.upper()
	if not Path(outformat_dir).is_dir():
		os.system("mkdir " + outformat_dir)
	output = outformat_dir + "/" + molecule + '.' + outformat
	
	spectra = normalize_spectra(spectra)

	colors     = itertools.cycle(('k', 'r', 'b'))
	linestyles = itertools.cycle(('-', '-', '-', ':', ':',':'))
	plt.figure(figsize=(10.8, 4.8))
	for spectrum in spectra:
		if not spectrum[3] == "Experimental data":
			cos_score       = cosine_distance(spectra[0][1], spectrum[1])
			statistics_file = output_dir + '/CSV/SINGLE/' + spectrum[3] + '_statistics.csv'
			with open(statistics_file, 'a') as csvfile:
			 	writer = csv.writer(csvfile, delimiter=',')
			 	writer.writerow([molecule, cos_score])
		if spectrum[3] == "OEP":
			spectrum[3] = "B3LYP/aug-cc-pVTZ"
		plt.plot(spectrum[0], spectrum[1], color=next(colors), linestyle=next(linestyles), label=spectrum[3])
	plt.legend(loc='upper right')
	plt.xlabel('Frequency, $cm^{-1}$')
	plt.ylabel('IR intensity')
	plt.yticks([])
	plt.savefig(output, format=outformat, bbox_inches='tight')
	plt.close()
	check_or_die(output, False)
	print('\n' + outformat.upper() + ' file saved at:', output)

def save_spectrum(exp_dir, qm_dir, qms, ff_dir, ffs, molecule, output_dir, gamma, png, pdf, svg):
	spectra = []
	
	exp_spectrum, start, stop, npoints = read_exp_data(exp_dir, molecule)
	spectra.append(exp_spectrum)
	
	for qm in qms:
		spectra.append(generate_spectrum(qm_dir, qm, molecule, None, start, stop, npoints, gamma))
	
	eigfreq_count = len(spectra[1][2])
	for ff in ffs:
		spectra.append(generate_spectrum(ff_dir, ff, molecule, eigfreq_count, start, stop, npoints, gamma))
	
	possible_formats         = ["png", "pdf", "svg"]
	desired_formats          = [ png,   pdf,   svg ]
	for i, desired_format in enumerate(desired_formats):
		if desired_format:
			selected_format = possible_formats[i]
			if selected_format in [ "png", "pdf", "svg" ]:
				save_spectra_as_figure(spectra, output_dir, molecule, selected_format)

def generate_spectrum_from_log(path, origin, start, stop, npoints, gamma):
	method_factors = {"G4": 0.965, "OEP": 0.968}
	scaling_factor = method_factors.get(origin, 1.0)
	print(scaling_factor)
	log = None
	for found_file in glob.glob('%s/*%s.log*' % (path, origin.lower())):
    		log = found_file
	if log:
		eigenfrequencies = []
		intensities      = []
		print('reading log file at:', log)
		if log.endswith(".gz"):
			lines = gzip.open(log, 'rt').readlines()
		else:
			lines = open(log, 'r').readlines()
		for line in lines:
			if "Frequencies" in line:
				values  = re.findall(r"[-+]?\d*\.\d+|\d+", line)
				for value in values:
					eigenfrequencies.append(float(value) * scaling_factor)
			elif "IR Inten" in line:
				values  = re.findall(r"[-+]?\d*\.\d+|\d+", line)
				for value in values:
					intensities.append(float(value))
	else:
        	sys.exit("The QM log file does not exist!!!")
	if eigenfrequencies and intensities:
		frequencies         = np.linspace(start, stop, npoints)
		intensity_all_modes = np.zeros(len(frequencies))
		for i, eigenfrequency in enumerate(eigenfrequencies):
			intensity_one_mode   = generate_cauchy_distribution(frequencies, eigenfrequency, gamma, intensities[i])
			intensity_all_modes += intensity_one_mode
		return [frequencies, intensity_all_modes, np.array(eigenfrequencies), origin]
	else:
		sys.exit("There are no frequencies and/or intensities in the QM log file!!!")

def read_exp_data(exp_dir, molecule):
	full_path = exp_dir + '/' + molecule + '.jdx'
	exp = None
	for found_file in glob.glob(full_path):
		exp = found_file
	if exp:
		intensities = []
		print('\n<EXP>\nreading experimental data file at:', exp)
		lines = open(exp, 'r').readlines()
		for line in lines:
			words  = line.split('=')
			if len(words) == 2:
				if "MAXX" in words[0]:
					stop = float(words[1])
				elif "MINX" in words[0]:
					start = float(words[1])
				elif "NPOINTS" in words[0]:
					npoints = int(words[1])
			else:
				words = line.split()
				if words[0].replace('.','',1).isdigit():
					for word in words[1:]:
						intensities.append(float(word))
		frequencies = np.linspace(start, stop, npoints)
		return [frequencies, np.array(intensities), None, "Experimental data"], start, stop, npoints
	else:
		sys.exit("The experimental data file does not exist!!!")
