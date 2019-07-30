import csv, os

path_head     = "/home/alfred/IR-Spectra/results/"
path_tail_old = "_statistics.csv"
path_tail_new = "_statistics_filtered.csv"

# find nan rows
nan_list = []
for ff in ["CGenFF", "GAFF-BCC", "GAFF-ESP"]:
	path_old = path_head + ff + path_tail_old
	csv_file = open(path_old,'r')
	for row in csv_file.readlines():
		values   = row.rstrip('\n').split(",")
		molecule = values[0]
		for value in values:
			if value == "nan" and molecule not in nan_list:
				nan_list.append(molecule)

# create new csv without nan
for ff in ["CGenFF", "GAFF-BCC", "GAFF-ESP"]:
	path_old = path_head + ff + path_tail_old
	path_new = path_head + ff + path_tail_new 
	with open(path_old, 'r') as inp, open(path_new, 'w') as out:
		writer = csv.writer(out)
		for row in csv.reader(inp):
			if row[0] not in nan_list:
				writer.writerow(row)
	os.remove(path_old)
	os.system("mv " + path_new + " " + path_old)
			

