#!/usr/bin/env python3

""" small script to create stats that show that starting with the initialiser is
better than starting with the SM / any first 5 steps
"""

def createStatsForInit():
    dbpath = "official.pcl"
    # from smodels.experimental.databaseObj import Database
    # db = Database ( dbpath )
    dbpath = "official"
    from walker.initialiser import Initialiser
    dictfile = "signal_database.dict"
    allowN1N1Prod = True
    initialiser = Initialiser ( walkerid = "stats",
        dictfile = dictfile, allowN1N1Prod = allowN1N1Prod,
        dbpath = dbpath )
    Kvalues = []
    for i in range(10):
        ret = initialiser.bestOfN(5)
        K = ret["K"]
        print ( "K {K}" )
        Kvalues.append( K )
    with open ( "init.stats", "wt" ) as f:
        f.write ( Kvalues + "\n" )
        f.close()

def plotStats():
    with open ( "init.stats", "rt" ) as f:
        Kvalues = eval ( f.read() )

if __name__ == "__main__":
    createStatsForInit()
