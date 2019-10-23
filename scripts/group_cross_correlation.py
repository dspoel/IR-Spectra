#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import sqlite3, csv, sys, argparse, math
from dbutils import *
from mol_csv_api import *
from organic import is_organic
from enum import Enum

def parseArguments():
    parser = argparse.ArgumentParser()
    parser.add_argument("-f", "--inputfile", dest="infile", help="Input file (database) for reading")
    parser.add_argument("-d", "--results_directory", dest="resdir", help="Directory containing correlation csv files")
    parser.add_argument("--verbose", help="write debug statements")
    args = parser.parse_args()
    return args

if __name__ == '__main__':
    args  = parseArguments()

    DB = DbUtils(args.infile, False)
    if (args.verbose):
        DB.set_debug()

    myclasses = []
    for row in DB.cursor.execute("select class from classification where website=1"):
        myclasses.append(str(row[0]))

    molsinclass = {}
    correctinclass = {}
    classinclass = {}

    for ff in [ "G4", "OEP", "GAFF-BCC", "GAFF-ESP", "CGenFF", "OPLS" ]:
        print("Processing %s" % ff)
        for prop in [ "pearson", "spearman" ]:
            for direction in [ "", "_inv" ]:
                #setting some parameters
                inname = args.resdir + "/CSV/" + ff + "_" + prop + direction + ".csv"
                outname = args.resdir + "/CSV/" + ff + "_" + prop + direction + "_grouped.csv" 
                outfile = open(outname, "w")
                outfile.write("Class|Molecules in Class|Correct Molecules|Fraction Molecules|Correct Class|Class Fraction\n")
                #resetting counters
                for cls in myclasses:
                    molsinclass[cls] = 0
                    correctinclass[cls] = 0
                    classinclass[cls] = 0

                #comparison
                lines = open(inname, "r").readlines()
                for line in lines:
                    words = line.split("|")
                    if (words[0][0:3] == "Exp") or (words[0][0:3] == "The"):
                        continue
                    else:
                        molecule1 = words[0]
                        molecule2 = words[-1][:-1]
                        mol1classes = []
                        mol2classes = []
                        for row in DB.cursor.execute(("""select cls.class
                                                         from   molecules      as mol, 
                                                                classification as cls,
                                                                link_mol_class as lmc
                                                         where  cls.classid  = lmc.classid  and 
                                                                lmc.molid    = mol.molid    and 
                                                                mol.filename = '%s.sdf'     and 
                                                                cls.website  = 1""" % molecule1)):
                            mol1classes.append(row[0])
                        for row in DB.cursor.execute(("""select cls.class
                                                         from   molecules      as mol, 
                                                                classification as cls,
                                                                link_mol_class as lmc
                                                         where  cls.classid  = lmc.classid  and 
                                                                lmc.molid    = mol.molid    and 
                                                                mol.filename = '%s.sdf'     and 
                                                                cls.website  = 1""" % molecule2)):
                            mol2classes.append(row[0])
                        for cls in mol1classes:
                            molsinclass[cls] += 1
                            if molecule1 == molecule2:
                                correctinclass[cls] += 1
                            if cls in mol2classes:
                                classinclass[cls] += 1

                #output
                for cls in myclasses: 
                    if molsinclass[cls] > 0:    
                        correctfrac = float(correctinclass[cls]) / float(molsinclass[cls])
                        classfrac = float(classinclass[cls]) / float(molsinclass[cls])
                        outfile.write("%s|%d|%d|%f|%d|%f\n" % (cls, molsinclass[cls], correctinclass[cls], correctfrac, classinclass[cls], classfrac))
                outfile.close()

