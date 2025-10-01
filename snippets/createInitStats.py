#!/usr/bin/env python3

""" small script to create stats that show that starting with the initialiser is
better than starting with the SM / any first 5 steps
"""

from ptools.helpers import py_dumps
from base.locker import lock, unlock

def writeModel( model : dict, outfile : str = "init5.stats" ):
    """ append model to outfile """
    models = []
    if os.path.exists ( outfile ):
        with open ( outfile, "rt" ) as f:
            models = eval(f.read())
    models.append ( model )
    lock ( outfile )
    with open ( outfile, "wt" ) as f:
        ds = py_dumps ( models, indent=4 )
        f.write ( f"{ds}\n" )
        f.close()
    unlock ( outfile )

def createStatsForInit( dbpath : str = "official.pcl", 
        outfile : str = "init5.stats", bestOfN : int = 5 ):
    # dbpath = "official.pcl"
    # from smodels.experimental.databaseObj import Database
    # db = Database ( dbpath )
    dbpath = "official"
    from walker.initialiser import Initialiser
    dictfile = "signal_database.dict"
    from base.runEnviron import RunEnviron
    environ = RunEnviron()
    initialiser = Initialiser ( walkerid = "stats",
        dictfile = dictfile, environ = environ )
    for i in range(1000):
        ret = initialiser.bestOfN(bestOfN)
        K, TL = ret["K"], ret["TL"]
        print ( f"K={K} TL={TL}" )
        writeModel ( ret, outfile )

if __name__ == "__main__":
    import argparse
    argparser = argparse.ArgumentParser(
                        description='small script to gather stats for initialiser comparison' )
    argparser.add_argument ( '-n', '--bestOfN',
            help='best of how many [5]',
            type=int, default=5 )
    argparser.add_argument ( '-o', '--outfile',
            help='output file ["init@@N@@.stats"]',
            type=str, default="init@@N@@.stats" )
    args = argparser.parse_args()
    outfile = args.outfile.replace( "@@N@@", str(args.bestOfN) )
    import os
    dirname = os.path.dirname ( os.path.abspath ( __file__ ) )
    os.chdir ( dirname ) # got to protomodels/snippets
    createStatsForInit( outfile = outfile, bestOfN = args.bestOfN )
