#!/usr/bin/env python3

import os, csv, itertools
import numpy as np
from numpy.random import choice 
import scipy.special
import scipy.stats
import math
from decimal import Decimal
from statistics import mean, stdev
import matplotlib.pyplot as plt
plt.switch_backend('agg')
plt.rcParams.update({'font.size': 11})

folders = [ "results_full", "results_medium", "results_short", "results_inverse" ]
myffs   = [ "G4", "OEP", "CGenFF", "GAFF-BCC", "GAFF-ESP", "OPLS" ]
methods = [ "pearson", "spearman" ]
directions = [ "", "_inv" ]
plotformat = 'pdf'

def read_correlations(folder, ff, method, direction):
    data_matrix = []
    with open("%s/CSV/%s_%s%s.csv" % (folder, ff, method, direction), "r") as infile:
        for line in infile.readlines():
            words = line.split("|")
            if (words[1][1] == ".") or (words[1][0] == "-"):
                data_matrix.append([float(word) for word in (words[1:-1])])
    if len(data_matrix) != len(data_matrix[-1]):
        print('Problem reading data')
        exit(1)
    return data_matrix

def count_ranks(matrix):
    ranks = []
    for i in range(len(matrix)):
       ranks.append(len(matrix) - scipy.stats.rankdata(matrix[i], method = 'max')[i]) 
    histo_data = []
    for i in range(len(matrix)):
        histo_data.append(list(ranks).count(i))
    return histo_data

def save_data(histos, myffs, outname):
   outfile = outname + ".csv"
   with open(outfile, 'w') as csvfile:
       writer = csv.writer(csvfile, delimiter = '|')
       firstline = ['Rank']
       for ff in myffs:
           firstline.append(ff)
       writer.writerow(firstline)
       for row in range(len(histos[myffs[0]])):
           line = []
           line.append(row + 1)
           for ff in myffs:
               line.append(histos[ff][row])
           writer.writerow(line)

def run_histogrammer(folders, myffs, methods, directions):
    histos = {}
    for folder in folders:
        for direction in directions:
            for method in methods:
                for ff in myffs:
                    full_matrix = read_correlations(folder, ff, method, direction)
                    histos[ff] = count_ranks(full_matrix)
                outname = folder + "/" + method + direction + "_rank_histo"
                save_data(histos, myffs, outname)
     
run_histogrammer(folders, myffs, methods, directions)

