#!/usr/bin/env python3

import os

# Update this
gmx  = "/Users/spoel/software-release-2019/bin/gmx_mpi_d"

def check_or_die(filenm, die):
    if not os.path.isfile(filenm):
        print("ERROR: generating " + filenm)
        if die:
            exit(1)

def run_one_nm(die, sigma, scale_factor, output):
    conf = "conf.gro"
    for minimize in [ "cg" ]:
        tpr  = minimize + ".tpr"
        os.system((gmx + " grompp -f ../../MDP/%s.mdp -o %s -v -maxwarn 1 -c conf.gro" % ( minimize, tpr)))
        check_or_die(tpr, die)
        conf = ("after_%s.g96" % minimize )
        os.system((gmx + " mdrun -s %s -c %s -v" % ( tpr, conf)))
        check_or_die(conf, die)

    tpr = "nm.tpr"
    os.system((gmx + " grompp -c %s -f ../../MDP/nm -o %s -v" % ( conf, tpr)))
    check_or_die(tpr, die)

    mtz = "nm.mtx"
    os.system(gmx + " mdrun -s %s -v -mtx %s" % ( tpr, mtz ))
    check_or_die(mtz, die )
    os.system(("%s nmeig -last 1000 -s %s -f %s -sigma %d > %s" % (gmx, tpr, mtz, sigma, output)))

run_one_nm(True, 1, 1, "output")
