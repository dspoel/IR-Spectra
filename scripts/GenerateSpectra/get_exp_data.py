#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse, os, sys, urllib.request, sqlite3, csv, math
from dbutils import *
from mol_csv_api import *
from organic import is_organic
from enum import Enum
from spectrum_functions import *
from pathlib import Path

if __name__ == "__main__":

	DB                   = DbUtils('/home/spoel/Liquids/DATABASE/molecules.sqlite3.dat', False)
	molecules            = []
	CAS_registry_numbers = []
	for row in DB.cursor.execute("SELECT filename, cas FROM molecules"):
		molecules.append(str(row[0]).split('.')[0])
		CAS_registry_numbers.append(str(row[1]))
	print(CAS_registry_numbers)
	print(len(CAS_registry_numbers))	
	

	outpath = "/home/alfred/JCAMP-DX/"
	for i, molecule in enumerate(molecules):
		CAS_registry_number = CAS_registry_numbers[i]
		numbers             = CAS_registry_number.split(';')
		for number in numbers:
			url      = 'https://webbook.nist.gov/cgi/cbook.cgi?ID=%s&Units=SI&cIR=on' % (number.replace(" ",""))
			response = urllib.request.urlopen(url)
			data     = response.read()
			text     = data.decode('utf-8')
			lines    = text.splitlines()
			for line in lines:
				if "Download <a" in line:
					words = line.split('"')
					url   = "https://webbook.nist.gov" + words[1].replace("&amp;","&")
					urllib.request.urlretrieve(url, outpath + molecule + ".jdx")
