#!/usr/bin/env python3

"""
Main code for submitting a walk
"""

import sys, time
sys.path.append('smodels')
sys.path.append('protomodels')
from walker.walkingWorker import main

t0 = time.time()
seed = 135
main(0, 1, "default", cheatcode = 0, dbpath = './toyWalk/toy-database',
    rundir = './toyWalk', maxsteps = 200, nevents = 100000,
    seed = seed, catchem = False )

print('Done in %1.1f min' %((time.time()-t0)/60.))