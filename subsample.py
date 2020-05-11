#!/usr/bin/env python3

import os, csv, itertools
import numpy as np
from numpy.random import choice 
import scipy.special
import math
from decimal import Decimal
from statistics import mean, stdev
import matplotlib.pyplot as plt
plt.switch_backend('agg')
plt.rcParams.update({'font.size': 11})

folders = [ "results_full", "results_medium"]#, "results_short", "results_inverse" ]
myffs   = [ "G4", "OEP", "CGenFF", "GAFF-BCC", "GAFF-ESP", "OPLS" ]
methods = [ "pearson", "spearman" ]
directions = [ "", "_inv" ]
plotformat = 'pdf'

def probability_all_random(identified, database_size):
    p_e = Decimal(1.0 / database_size)
    prob = Decimal(0.0)
    prob_part1 = Decimal(scipy.special.binom(database_size, identified)) * p_e**identified * (1 - p_e)**(database_size - identified) 
    prob += prob_part1 
    return prob

def power_all(identified, database_size):
    powerfrac = Decimal(0.0)
    if database_size > 1:
        for ll in range(database_size - identified + 1):
            powerfrac += probability_all_random((database_size-ll), database_size)
        if powerfrac < 1E-300:
            powerval = 300.0
        else:
            powerval = -math.log10(powerfrac)
    else:
        powerval = 0.0
    return powerval

def create_subsamples(subsize, fullsize, number):
    if subsize > fullsize:
        print('Impossible size of subsample')
        exit(1)
    if number > scipy.special.binom(fullsize, subsize):
        print('Impossible number of subsamples')
        exit(1)
    i = 0
    newsub = []
    allsubs = []
    while i < number: 
        newsub = np.sort(choice(fullsize, subsize, replace=False))
        containcheck = False
        for k in range(len(allsubs)):
            if all(newsub == allsubs[k]):
                containcheck = True
        if not containcheck:
            allsubs.append(newsub)
            i += 1
    return allsubs

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

def create_submatrix(data_matrix, index_list):
    submatrix = []
    for i in index_list:
        submatrix.append(list(np.array(data_matrix[i])[index_list]))
    return submatrix

def count_matches(submatrix):
    matches = 0
    for i in range(len(submatrix)):
        if submatrix[i][i] == max(submatrix[i]):
            matches += 1
    return matches

def plot_results(xs, y1s, y2s, powery1s, powery2s, myffs, outname, form):
    lines      = itertools.cycle(('-', '--', '-.', ':', (0, (4, 1, 1, 1, 1, 1)), (0, (4, 4)) ))
    colors     = itertools.cycle(('k', '#1b9e77', '#d95f02', '#7570b3', '#e7298a', '#66a61e', '#e6ab02'))
    plt.figure(figsize=(10.8, 6.75))
    plt.rc('font', size=21)
    for ff in myffs:
        if ff == "OEP":
            figlabel = "B3LYP/aug-cc-pVTZ"
        elif ff == "G4":
            figlabel = "B3LYP/6-31G(2df,p)"
        else:
            figlabel = ff        
        x = np.array(xs[ff])
        y1 = np.array(y1s[ff])
        y2 = np.array(y2s[ff])
        lc = next(colors)
        plt.plot(x, y1, color=lc, linestyle=next(lines), label=figlabel)
        plt.fill_between(x, y1-y2, y1+y2, color=lc, alpha=0.25)#, label=figlabel)
    plt.xlabel('Subsample Size' , size = 24)
    plt.ylabel('Matches', size = 24)
    plt.legend(loc='upper left', fontsize=18) 
    outfile = outname + "_abs.pdf"
    plt.savefig(outfile, format=form, bbox_inches='tight')
    plt.close()
    colors     = itertools.cycle(('k', '#1b9e77', '#d95f02', '#7570b3', '#e7298a', '#66a61e', '#e6ab02'))
    lines      = itertools.cycle(('-', '--', '-.', ':', (0, (4, 1, 1, 1, 1, 1)), (0, (4, 4)) ))
    plt.figure(figsize=(10.8, 6.75))
    plt.rc('font', size=21)
    for ff in myffs:
        if ff == "OEP":
            figlabel = "B3LYP/aug-cc-pVTZ"
        elif ff == "G4":
            figlabel = "B3LYP/6-31G(2df,p)"
        else:
            figlabel = ff        
        x = np.array(xs[ff])
        y1 = np.array(y1s[ff])
        y2 = np.array(y2s[ff])
        lc = next(colors)
        plt.plot(x, y1/x, color=lc, linestyle=next(lines), label=figlabel)
        plt.fill_between(x, (y1-y2)/x, (y1+y2)/x, color=lc, alpha=0.25)#, label=figlabel)
    plt.xlabel('Subsample Size' , size = 24)
    plt.ylabel('Match Frequency', size = 24)
    plt.legend(loc='upper right', fontsize='small') 
    outfile = outname + "_rel.pdf"
    plt.savefig(outfile, format=form, bbox_inches='tight')
    plt.close()

    colors     = itertools.cycle(('k', '#1b9e77', '#d95f02', '#7570b3', '#e7298a', '#66a61e', '#e6ab02'))
    lines      = itertools.cycle(('-', '--', '-.', ':', (0, (4, 1, 1, 1, 1, 1)), (0, (4, 4)) ))
    plt.figure(figsize=(10.8, 6.75))
    plt.rc('font', size=21)
    for ff in myffs:
        if ff == "OEP":
            figlabel = "B3LYP/aug-cc-pVTZ"
        elif ff == "G4":
            figlabel = "B3LYP/6-31G(2df,p)"
        else:
            figlabel = ff        
        x = np.array(xs[ff])
        powery1 = np.array(powery1s[ff])
        powery2 = np.array(powery2s[ff])
        lc = next(colors)
        plt.plot(x, powery1, color=lc, linestyle=next(lines), label=figlabel)
        plt.fill_between(x, (powery1-powery2), (powery1+powery2), color=lc, alpha=0.25)#, label=figlabel)
    plt.xlabel('Subsample Size' , size = 24)
    plt.ylabel('Predictive Power', size = 24)
    plt.legend(loc='upper right', fontsize='small') 
    outfile = outname + "_pow.pdf"
    plt.savefig(outfile, format=form, bbox_inches='tight')
    plt.close()

def save_data(xs, y1s, y2s, powery1s, powery2s, myffs, outname):
   outfile = outname + "_data.csv"
   with open(outfile, 'w') as csvfile:
       writer = csv.writer(csvfile, delimiter = '|')
       firstline = ['Subsample Size']
       for ff in myffs:
           firstline.append(ff + 'm')
           firstline.append(ff + 'stdev')
           firstline.append(ff + 'power')
           firstline.append(ff + 'powerstdev')
       writer.writerow(firstline)
       for row in range(len(xs[myffs[0]])):
           line = []
           line.append(xs[myffs[0]][row])
           for ff in myffs:
               line.append(y1s[ff][row])
               line.append(y2s[ff][row])
               line.append(powery1s[ff][row])
               line.append(powery2s[ff][row])
           writer.writerow(line)

def run_submatch(folders, myffs, methods, directions, samplesize, plotformat):
    sizelists = {}
    meanlists = {}
    stdevlists= {}
    powermeans= {}
    powerstdevs= {}
    for folder in folders:
        for direction in directions:
            for method in methods:
                for ff in myffs:
                    full_matrix = read_correlations(folder, ff, method, direction)
                    sizelist = []
                    meanlist = []
                    stdevlist = []
                    powermean= []
                    powerstdev= []
                    for subsize in range(len(full_matrix)-1):
                        matchlist = []
                        powerlist = []
                        subsample_list = create_subsamples(subsize+1, len(full_matrix), samplesize)
                        for subsample in subsample_list:
                            submatrix = create_submatrix(full_matrix, subsample)
                            matches = count_matches(submatrix)
                            matchlist.append(matches)
                            powerlist.append(power_all(matches, subsize+1))
                        sizelist.append(subsize+1)
                        meanlist.append(mean(matchlist))
                        stdevlist.append(stdev(matchlist))
                        powermean.append(mean(powerlist))
                        powerstdev.append(stdev(powerlist))
                    sizelists[ff] = sizelist
                    meanlists[ff] = meanlist
                    stdevlists[ff] = stdevlist
                    powermeans[ff] = powermean
                    powerstdevs[ff] = powerstdev
                outname = folder + "/" + method + direction + "_subsampling"
                plot_results(sizelists, meanlists, stdevlists, powermeans, powerstdevs, myffs, outname, plotformat) 
                save_data(sizelists, meanlists, stdevlists, powermeans, powerstdevs,myffs, outname)
     
run_submatch(folders, myffs, methods, directions, 150, plotformat)

