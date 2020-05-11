#!/usr/bin/env python3

import os
import numpy as np
import scipy.special
import math
from decimal import Decimal

myffs   = [ "G4", "OEP", "CGenFF", "GAFF-BCC", "GAFF-ESP", "OPLS" ]
methods = [ "pearson", "spearman" ]
directions = [ "", "_inv" ]
latex   = {}
latex[methods[0]] = "Pearson $r$"
latex[methods[1]] = "Spearman $\\rho$"
metkey = {}
metkey[methods[0]] = "Pearson"
metkey[methods[1]] = "Spearman"

def probability_from_random(identified, class_size, database_size):
    #print(identified, class_size, database_size)
    p_e = Decimal(class_size)/Decimal(database_size)
    cl_fac = np.math.factorial(class_size)
    prob = Decimal(0.0)
    for i in range(database_size - class_size + 1):
        inclass = identified + i
        prob_part1 = Decimal(scipy.special.binom(database_size, inclass)) * p_e**inclass * (1 - p_e)**(database_size - inclass) 
        prob_part2 = Decimal(scipy.special.binom(class_size, identified)) * Decimal(np.math.factorial(inclass)) / Decimal(np.math.factorial(inclass - identified)) 
        prob_part3 = Decimal(np.math.factorial(database_size - inclass)) / Decimal(np.math.factorial(database_size + identified - inclass - class_size)) 
        prob_part4 = Decimal(np.math.factorial(database_size)) / Decimal(np.math.factorial(database_size - class_size))
        prob += prob_part1 * prob_part2 * prob_part3 / prob_part4
    return prob

def probability_all_random(identified, database_size):
    p_e = Decimal(1.0 / database_size)
    prob = Decimal(0.0)
    prob_part1 = Decimal(scipy.special.binom(database_size, identified)) * p_e**identified * (1 - p_e)**(database_size - identified) 
    prob += prob_part1 
    return prob


def read_class_data():
    alldata = {}
    for direction in directions:
        dirdata = {}
        for ff in myffs:
            ffdata = {}
            ffdata["all_molecules"] = {}
            for method in methods:
                with open("%s_%s%s_grouped.csv" % (ff, method, direction), "r") as infile:
                    for line in infile.readlines():
                        words = line.split("|")
                        if words[0] == "Class":
                            continue
                        if not words[0] in ffdata:
                            ffdata[words[0]] = {}
                        N = int(words[1])
                        ffdata[words[0]][method] = (N, int(words[2]), float(words[3]), int(words[4]), float(words[5]))
                with open("matching.txt", "r") as infile:
                    for line in infile.readlines():
                        if len(line) > 1:
                            words = line.split(" ")
                            if direction == "":
                                if words[0] == ff:
                                    if words[3] == metkey[method]:
                                        N = int(words[-8][:-1])
                                        ffdata["all_molecules"][method] = (N, int(words[-10]))
                            else:
                                if words[2] == ff:
                                    if words[3] == metkey[method]:
                                        N = int(words[-8][:-1])
                                        ffdata["all_molecules"][method] = (N, int(words[-10]))
            dirdata[ff] = ffdata
        alldata[direction] = dirdata
    return alldata # format: alldata[direction][ff][method]


def write_table(alldata, isClass):
    collectfile = open("group_correlations.tex", "w")
    printN = isClass
    database_size = alldata[""][myffs[0]]["all_molecules"][methods[0]][0]
    #db_fac = np.math.factorial(database_size)
    for ff in myffs:
        for direction in directions:
            fflabel = {}
            if ff == "OEP":
                fflabel[ff] = "B3LYP/aug-cc-pVTZ"
            elif ff == "G4":
                fflabel[ff] = "B3LYP/6-31G(2df,p)"
            else:
                fflabel[ff] = ff
            basename = ff + direction
            with open(basename+".tex", "w") as outf:
                collectfile.write("\\input{"+basename+".tex}\n")
                outf.write("\\begin{longtable}")
                if printN:
                    character = "c"
                else:
                    character = ""
                fields = (len(methods) * "|ccccc") 
                outf.write("{l%s%s}\n" % (character, fields))
                if direction == "_inv":
                    dirtext = "experimental spectra against database of theoretical spectra"
                else:
                    dirtext = "thoeretical spectra against database of experimental spectra"
                outf.write("\\caption{Correlation of %s for %s sorted according to compound classification.}\n" % (dirtext, fflabel[ff]))
                outf.write("\\label{%s}\\\\\n" % basename)
                outf.write("\\hline\n")
                for i in range(2):
                    outf.write("\\multirow{3}{*}{Class}")
                    if printN:
                        outf.write(" & \\multirow{3}{*}{N}")
                    for method in methods:
                        if method == methods[-1]:
                            outf.write("& \\multicolumn{5}{c}{%s}" % latex[method])
                        else:
                            outf.write("& \\multicolumn{5}{c|}{%s}" % latex[method])
                    outf.write("\\\\\n")
                    if printN:
                        outf.write(" & ")
                    for method in methods:
                        if method == methods[-1]:
                            outf.write("& \\multicolumn{2}{c}{Molecules} & \multicolumn{3}{c}{Classes}")
                        else:
                            outf.write("& \\multicolumn{2}{c}{Molecules} & \multicolumn{3}{c|}{Classes}")
                    outf.write("\\\\\n")
                    if printN:
                        outf.write(" & ")
                    for method in methods:
                        outf.write("& Matched & Fraction " * 2 + "& Power")
                    outf.write("\\\\\n")
                    outf.write("\\hline\n")
                    if i == 0:
                        outf.write("\\endfirsthead\n")
                    else:
                        outf.write("\\endhead\n")
                    outf.write("\\hline\n")
                    outf.write("\\endfoot\n")

                for kkk in range(2):
                    for mol in sorted(alldata[""][myffs[0]].keys()):
                        mmm = mol.replace("_", " ")
                        if (kkk == 0 and mmm != "all molecules"):
                            outf.write("%s" % (mmm ) )
                            if printN:
                                 N = alldata[""][myffs[0]][mol][methods[0]][0]
                                 outf.write("& %d " % N) 
                            for method in methods:
                                for k in range(2):
                                    valuenum  = alldata[direction][ff][mol][method][2*k+1]
                                    formatnum = "%d"
                                    valuefrac  = alldata[direction][ff][mol][method][2*k+2]
                                    formatfrac = "%.3f"
                                    formatstr = "& " + formatnum +" & " + formatfrac
                                    if k==1:
                                        powerfrac = Decimal(0.0)
                                        #powerfrac = valuefrac * alldata[""][myffs[0]]["all_molecules"][methods[0]][0] / alldata[""][myffs[0]][mol][methods[0]][0]  
                                        #mixedfrac = valuefrac * np.sqrt(alldata[""][myffs[0]]["all_molecules"][methods[0]][0] / alldata[""][myffs[0]][mol][methods[0]][0])
                                        #formatstr += " & %.2f & %.2f " 
                                        for ll in range(N - valuenum + 1):
                                            powerfrac += probability_from_random((N-ll), N, alldata[""][myffs[0]]["all_molecules"][methods[0]][0])
                                        if powerfrac < 1E-300:
                                            powerval = 300.0
                                        else:
                                            powerval = -math.log10(powerfrac)
                                        #print(powerval)
                                        formatstr += "& %.1f " 
                                        #outf.write(formatstr % (valuenum, valuefrac, powerfrac, mixedfrac))
                                        outf.write(formatstr % (valuenum, valuefrac, powerval))
                                    else:
                                        outf.write(formatstr % (valuenum, valuefrac))
                            outf.write("\\\\\n")
                        if (kkk == 1 and mmm == "all molecules"):
                            outf.write("%s" % mmm)
                            if printN:
                                outf.write("& %d " % alldata[""][myffs[0]][mol][methods[0]][0])
                            N = alldata[""][myffs[0]]["all_molecules"][methods[0]][0]
                            for method in methods:
                                valuenum = alldata[direction][ff][mol][method][1]
                                valuefrac = float(valuenum) / float(N)
                                powerfrac = Decimal(0.0)
                                for ll in range(N - valuenum + 1):
                                    powerfrac += probability_all_random((N-ll), N)
                                print((N-ll),N, powerfrac)
                                if powerfrac < 1E-300:
                                    powerval = 300.0
                                else:
                                    powerval = -math.log10(powerfrac)
                                outf.write("& %d & %.3f & --- & --- & %.1f " % (valuenum, valuefrac, powerval))
                            outf.write("\\\\\n")
                outf.write("\\hline\n")
                outf.write("\\end{longtable}\n")

alldata = read_class_data()
write_table(alldata, True)
