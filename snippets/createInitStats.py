#!/usr/bin/env python3

""" small script to create stats that show that starting with the initialiser is
better than starting with the SM / any first 5 steps
"""

def writeModels( models ):
    with open ( "init.stats", "wt" ) as f:
        f.write ( f"{models}\n" )
        f.close()

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
    models = []
    for i in range(1000):
        ret = initialiser.bestOfN(5)
        K, TL = ret["K"], ret["TL"]
        print ( f"K={K} TL={TL}" )
        models.append( ret )
        writeModels ( models )

def plotStats():
    with open ( "init.stats", "rt" ) as f:
        models = eval ( f.read() )

if __name__ == "__main__":
    createStatsForInit()
