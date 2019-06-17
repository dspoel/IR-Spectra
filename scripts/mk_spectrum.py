#!/usr/bin/env python3

import os, glob

dirs = [ "G4", "CGenFF", "GAFF-ESP" ]

def get_mols():
    mols = {}
    for dir in dirs:
        os.chdir(dir)
        for m in glob.glob("*"):
            spectrum = m + "/spectrum.xvg"
            if os.path.isfile(spectrum):
                if m in mols:
                    mols[m] += 1
                else:
                    mols[m] = 1
        os.chdir("..")
    return mols

def make_plot(mol, target):
    if (os.path.isfile(target)):
        return
    xtmp = "x.xvg"
    os.system("cp template.xvg %s" % xtmp)
    i = 0
    with open(xtmp, "a") as xxx:
        for d in dirs:
            xxx.write("@target g%d.s0\n" % i)
            xxx.write("@type xy\n")
            with open(d + "/" + mol + "/spectrum.xvg", "r") as inf:
                for l in inf:
                    xxx.write(l)
            xxx.write("&\n")
            i += 1
    if i == len(dirs):
        os.system("xmgrace -hardcopy -hdevice PNG -printfile %s %s" % ( target, xtmp) )

base_dir =  "../../VCHEM/SPECTRUM/"
os.makedirs(base_dir, exist_ok=True)
mols = get_mols()
for m in mols:
    if mols[m] == len(dirs):
        target = base_dir + m + ".png"
        make_plot(m, target)
