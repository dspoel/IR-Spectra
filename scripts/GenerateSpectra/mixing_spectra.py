#!/usr/bin/env python3

import os, sys
from spectrum_functions import *
from pathlib import Path

exp_dir = "/home/alfred/JCAMP-DX/ABSORBANCE"
qm_dir = "/home/spoel/Liquids/MOLECULES"
qms = ["G4", "OEP"]
ff_dir = "/home/spoel/wd/THERMO"
ffs = ["CGenFF", "GAFF-BCC", "GAFF-ESP"]
gamma = 24

molecules = find_molecules(exp_dir, qm_dir, qms, ff_dir, ffs)

print('\nThe following number of molecules were found in all listed directories and will be processed:', len(molecules))

scores = []
statistics_file = '/home/alfred/IR-Spectra/results/CSV/SINGLE/Mixed_statistics.csv'
with open(statistics_file, 'w') as csvfile:
	writer = csv.writer(csvfile, delimiter=',')
	writer.writerow(['molecule', 'cos'])
	for molecule in molecules:
		print('\nNOW PROCESSING:', molecule)
	
		# Read experimental data
		exp_spectrum, start, stop, npoints = read_exp_data(exp_dir, molecule)
	

		# Read QM data
		intensities = []
		log = None
		query = '%s/G4/%s/*.log*' % (qm_dir, molecule)
		for found_file in glob.glob(query):
			log = found_file
		if log:
			print('reading log file at:', log)
			if log.endswith(".gz"):
				lines = gzip.open(log, 'rt').readlines()
			else:
				lines = open(log, 'r').readlines()
			for line in lines:
				if "IR Inten" in line:
					values  = re.findall(r"[-+]?\d*\.\d+|\d+", line)
					for value in values:
						intensities.append(float(value))
		eigfreq_count = len(intensities)


		# Read FF eigenfrequencies
		eigenfrequencies = extract_eigenfrequencies("%s/CGenFF/%s" % (ff_dir, molecule))
		eigenfrequencies = eigenfrequencies[-eigfreq_count:]

		# Generate mixed spectrum
		frequencies    = np.linspace(start, stop, npoints)
		mixed_intensities = np.zeros(len(frequencies))
		for i, eigenfrequency in enumerate(eigenfrequencies):
			mixed_intensities += generate_cauchy_distribution(frequencies, eigenfrequency, gamma, intensities[i])
		mixed_spectrum = [frequencies, mixed_intensities, None, "Mixed"]        
	
		# Normalize
		spectra = [exp_spectrum, mixed_spectrum]
		spectra = normalize_spectra(spectra)

		# Score the spectrum!
		score = cosine_distance(spectra[0][1],spectra[1][1])
		scores.append(score)

		output = "/home/alfred/IR-Spectra/results/MIXED/%s.png" % (molecule)

		colors     = itertools.cycle(('k', 'r'))
		linestyles = itertools.cycle(('-', '-'))
		plt.figure(figsize=(10.8, 4.8))
		for spectrum in spectra:
			plt.plot(spectrum[0], spectrum[1], color=next(colors), linestyle=next(linestyles), label=spectrum[3])
		plt.legend(loc='upper right')
		plt.xlabel('Frequency, $cm^{-1}$')
		plt.ylabel('IR intensity')
		plt.yticks([])
		plt.savefig(output, format="PNG", bbox_inches='tight')
		plt.close()
		check_or_die(output, False)
		writer.writerow([molecule, score])		

print("Average score for all molecules:", sum(scores)/len(scores))
