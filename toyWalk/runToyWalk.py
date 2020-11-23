#!/usr/bin/env python3

"""
Main code for submitting a walk
"""

import sys, time
sys.path.append('../smodels')
sys.path.append('../')
from walker.walkingWorker import main

t0 = time.time()
seed = 135
main(0, 1, "default", cheatcode = 0, dbpath = './toy-database',
    rundir = './results', maxsteps = 200, nevents = 100000,
    seed = seed, catchem = False, record_history = True )

print('Done in %1.1f min' %((time.time()-t0)/60.))
