#!/usr/bin/env python3

"""
Collect signficances and p values of UL results 
"""
    
from smodels.experiment.databaseObj import Database
from smodels.tools.physicsUnits import GeV
import random, pickle, numpy
import scipy.stats

def collect():
    db = Database ( "../../smodels-database" )
    ers = db.getExpResults ( dataTypes = [ "upperLimit" ], onlyWithExpected=True )
    allSs = []
    for er in ers:
        txnlist = er.datasets[0].txnameList
        for txn in txnlist:
            for i in range(1000):
                m1 = random.uniform ( 100, 2500 )
                m3 = random.uniform ( 0, m1 )
                m = [[ m1*GeV, m3*GeV], [ m1*GeV, m3*GeV] ]
                ul,eul=None,None
                try:
                    ul = txn.getULFor(m, False )
                    eul = txn.getULFor(m, True )
                except Exception:
                    m2 = random.uniform ( 0, 2500 )
                    m = [[ m1*GeV, m2*GeV, m3*GeV], [ m1*GeV, m2*GeV, m3*GeV] ]
                    try:
                        ul = txn.getULFor(m, False )
                        eul = txn.getULFor(m, True )
                    except Exception:
                        pass
                if type(ul) == type(None) or type(eul) == type(None):
                    continue
                sigma = eul / 1.96
                S = float ( ( ul - eul ) / sigma )
                allSs.append ( S )
                # print ( "->", er.globalInfo.id, txn, S )
    print ("all", allSs )
    print ("all", max(allSs) )
    f=open("ulSs.pcl","wb")
    pickle.dump(allSs,f)
    f.close()

def read():
    f=open("ulSs.pcl","rb")
    allSs=pickle.load(f)
    f.close()
    return allSs

def computeP ( allSs ):
    ret = []
    for s in allSs:
        ret.append ( scipy.stats.norm.cdf ( s ) )
    return ret

def plotS ( allSs ):
    from matplotlib import pyplot as plt
    plt.hist ( allSs, bins=numpy.arange(-2,3,.1) )
    plt.xlabel ( "significance Z" )
    plt.savefig ( "ulSs.png") 

def plotP ( ps ):
    from matplotlib import pyplot as plt
    plt.clf()
    plt.hist ( ps, bins=numpy.arange(0,1.01,.05) )
    plt.xlabel ( "p-values" )
    plt.savefig ( "ulPs.png") 

if __name__ == "__main__":
    # collect()
    allSs = read()
    ps = computeP ( allSs )
    plotS ( allSs )
    plotP ( allSs )
