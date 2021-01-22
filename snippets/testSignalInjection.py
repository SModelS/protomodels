#!/usr/bin/env python3

""" just test that the fake signal gets injected in the right part of the UL maps """

import math
from smodels.experiment.databaseObj import Database
from smodels.tools.physicsUnits import GeV

# ./ptools/fetchFromClip.py -R rundir.frozen1 --database
db=Database("default.pcl")
print ( "db", db.databaseVersion )
er = db.getExpResults( [ "CMS-SUS-19-006" ] )[0]
ds= er.datasets[0]
print ( ds.txnameList )
txn = ds.txnameList[6]
print ( txn )

def distance ( mass ):
    # return math.sqrt ( (mass[0] - 735.)**2 + (mass[1]-162.6)**2  )
    return math.sqrt ( (mass[0] - 1166.)**2 + (mass[1]-162.6)**2  )

masses = []
for mLSP in range ( 100, 240, 50 ):
    for msquark in range ( 850, 1300, 50 ):
        masses.append (  [msquark,mLSP] )
        # masses.append (  [[msquark*GeV,mLSP*GeV],[msquark*GeV,mLSP*GeV]] )
for mass in masses:
    mvec = [ [ mass[0]*GeV, mass[1]*GeV ], [ mass[0]*GeV, mass[1]*GeV ] ]
    oUL=txn.getULFor(mvec,expected=False)
    eUL=txn.getULFor(mvec,expected=True)
    congr=False
    if abs(eUL - oUL ) / ( eUL + oUL ) < 1e-5:
        congr=True
    scongr = "not injected"
    if congr == False:
        scongr = "injected"
    print ( f"mass {mass}: {distance(mass)}: {scongr}" )
