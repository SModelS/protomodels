#!/usr/bin/env python3

"""
Collect signficances and p values of UL results 
"""
    
import smodels.experiment.txnameObj ## gridpoints!!
smodels.experiment.txnameObj.TxNameData._keep_values = True
from smodels.experiment.databaseObj import Database
from smodels.tools.physicsUnits import GeV, fb, pb
import random, pickle, numpy, sys
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
    print ("all", min(allSs), numpy.mean(allSs), max(allSs) )
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
    ret = []
    for s in allSs:
        ret.append ( scipy.stats.norm.cdf ( s ) )
    return ret

def plotS ( allSs ):
    from matplotlib import pyplot as plt
    plt.hist ( allSs, bins=numpy.arange(-2,3,.1) )
    plt.title ( "significances from upper limits, SModelS 2.0.0-beta" )
    plt.xlabel ( "significance Z" )
    plt.savefig ( "ulSs.png") 

def plotP ( ps ):
    from matplotlib import pyplot as plt
    plt.clf()
    plt.hist ( ps, bins=numpy.arange(0,1.01,.05) )
    plt.title ( "$p$-values from upper limits, SModelS 2.0.0-beta" )
    plt.xlabel ( "$p$-values" )
    plt.savefig ( "ulPs.png") 

if __name__ == "__main__":
    # collect()
    allSs = read()
    print ( f"{len(allSs)} points total" )
    nAS = numpy.array ( allSs )
    print ( f"{len(nAS[nAS>0])} points > 0" )
    print ( f"{len(nAS[nAS>2])} points > 2" )
    print ( f"{len(nAS[nAS>2.5])} points > 2.5" )
    print ( f"{len(nAS[nAS>3])} points > 3" )
    ps = computeP ( allSs )
    plotS ( allSs )
    plotP ( allSs )
