#!/usr/bin/env python3

"""
Collect signficances and p values of UL results 
"""
    
import smodels.experiment.txnameObj ## gridpoints!!
smodels.experiment.txnameObj.TxNameData._keep_values = True
from smodels.experiment.databaseObj import Database
from smodels.base.physicsUnits import GeV, fb, pb
import random, pickle, sys, time
import numpy as np
import scipy.stats

def collect():
    db = Database ( "./database.pcl" ) # , force_load = "txt" )
    ers = db.getExpResults ( dataTypes = [ "upperLimit" ], onlyWithExpected=True )
    allSs = []
    for er in ers:
        txnlist = er.datasets[0].txnameList
        for txn in txnlist:
            ct=0
            origdata = eval(txn.txnameData.origdata)
            for point in origdata:
                m = point[0]
                rul = point[1]
                ul,eul=None,None
                try:
                    ul = txn.getULFor(m, False )
                    eul = txn.getULFor(m, True )
                except Exception:
                    pass
                if type(ul) == type(None) or type(eul) == type(None):
                    continue
                sigma = eul / 1.96
                S = float ( ( ul - eul ) / sigma )
                if (S < -1.8 or S > 3.5) and ct<3:
                # if S > 10. and ct<3:
                    print ( )
                    print ( "S=%.2f for ul=%s, eul=%s sigma=%s" % ( S, ul, eul, sigma ) )
                    print ( "  at ", er.globalInfo.id, txn.txName, m, "rul", rul )
                    ct += 1
                allSs.append ( S )
                # print ( "->", er.globalInfo.id, txn, S )
    print ("all", min(allSs), np.mean(allSs), max(allSs) )
    f=open("ulSs.pcl","wb")
    pickle.dump(allSs,f)
    f.close()
    sys.exit()

def read():
    f=open("ulSs.pcl","rb")
    allSs=pickle.load(f)
    f.close()
    return allSs

def computeP ( allSs ):
    return scipy.stats.norm.cdf ( allSs )
    """
    ret = []
    for s in allSs:
        ret.append ( scipy.stats.norm.cdf ( s ) )
    return ret
    """

def plotS ( allSs ):
    from matplotlib import pyplot as plt
    xmin, xmax = -1.5, 2.5
    bins = np.arange( xmin,xmax ,.1)
    clipped = np.clip( allSs, bins[0], bins[-1] )
    print ( "stats for Z:", np.mean(clipped),"+-",np.std(clipped) )
    r = plt.hist ( clipped, bins=bins )
    mr0 = max(r[0])
    plt.plot ( bins, [ mr0*scipy.stats.norm.pdf(x) for x in bins ] )
    plt.title ( "significances from upper limits, SModelS 2.0.0-beta" )
    plt.text(xmax+(xmax-xmin)*0.032,.1,time.asctime(),c="grey", rotation=90 )
    plt.plot ( [ 0., 0. ], [ 0, mr0 ], linestyle="-." )
    plt.xlabel ( "significance Z" )
    plt.savefig ( "ulSs.png") 

def plotP ( ps ):
    from matplotlib import pyplot as plt
    plt.clf()
    plt.hist ( ps, bins=np.arange(0,1.01,.05) )
    plt.title ( "$p$-values from upper limits, SModelS 2.0.0-beta" )
    plt.text(1.08,.1,time.asctime(),c="grey", rotation=90 )
    plt.xlabel ( "$p$-values" )
    plt.savefig ( "ulPs.png") 

if __name__ == "__main__":
    # collect()
    allSs = read()
    print ( f"{len(allSs)} points total" )
    nAS = np.array ( allSs )
    print ( f"{len(nAS[nAS>0])} points > 0" )
    print ( f"{len(nAS[nAS>2])} points > 2" )
    print ( f"{len(nAS[nAS>2.5])} points > 2.5" )
    print ( f"{len(nAS[nAS>3])} points > 3" )
    ps = computeP ( allSs )
    plotS ( allSs )
    plotP ( allSs )
