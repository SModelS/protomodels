#!/usr/bin/env python3

""" small script to create stats that show that starting with the initialiser is
better than starting with the SM / any first 5 steps
"""

from ptools.helpers import py_dumps
from base.locker import lock, unlock
import time

def writeModel( model : dict, outfile : str = "init5.stats" ):
    """ append model to outfile """
    models = []
    if os.path.exists ( outfile ):
        with open ( outfile, "rt" ) as f:
            models = eval(f.read())
    models.append ( model )
    lock ( outfile )
    with open ( outfile, "wt" ) as f:
        print ( f"[createInitStats] writing model to {outfile}" )
        ds = py_dumps ( models, indent=4 )
        f.write ( f"{ds}\n" )
        f.close()
    unlock ( outfile )

def createStatsForInit( args : dict ):
    """ create the stats for the initialiser

    args:
        - dbpath (os.PathLike): path to database
        - outfile (os.PathLike): write to this file
        - bestOfN (int): best of how many attempts
    """
    dbpath = args["dbpath"]
    outfile = args["outfile"]
    bestOfN = args["bestOfN"]
    outfile = outfile.replace( "@@N@@", str(bestOfN) )
    dictfile = args["dictfile"]
    # dbpath = "official.pcl"
    #from smodels.experimental.databaseObj import Database
    # db = Database ( dbpath )
    # dbpath = "official"
    from walker.initialiser import Initialiser
    # dictfile = "signal_database.dict"
    from base.runEnviron import RunEnviron
    environ = RunEnviron()
    initialiser = Initialiser ( walkerid = "stats",
        dictfile = dictfile, environ = environ )
    for i in range(1000):
        t0 = time.time()
        ret = initialiser.bestOfN(bestOfN)
        dt = time.time() - t0
        K, TL = ret["K"], ret["TL"]
        ret[ "dt" ] = dt
        print ( f"K={K} TL={TL}" )
        writeModel ( ret, outfile )

def createRunDict():
    """ create a simple, default run.dict file """
    from base.runEnviron import RunEnviron
    RunEnviron.create()

if __name__ == "__main__":
    import argparse
    argparser = argparse.ArgumentParser(
                        description='small script to gather stats for initialiser comparison' )
    argparser.add_argument ( '-n', '--bestOfN',
            help='best of how many [5]',
            type=int, default=5 )
    argparser.add_argument ( '-c', '--create_rundict',
            help='create a default run.dict, then exit', action="store_true" )
    argparser.add_argument ( '-d', '--dbpath',
            help='databasse path ["official"]',
            type=str, default="official" )
    argparser.add_argument ( '-D', '--dictfile',
            help='path to database dict file ["signal_database.dict"]',
            type=str, default="signal_database.dict" )
    argparser.add_argument ( '-o', '--outfile',
            help='output file ["init@@N@@.stats"]',
            type=str, default="init@@N@@.stats" )
    args = argparser.parse_args()
    if args.create_rundict:
        createRunDict()
        import sys; sys.exit()
    import os
    dirname = os.path.dirname ( os.path.abspath ( __file__ ) )
    os.chdir ( dirname ) # got to protomodels/snippets
    createStatsForInit( vars(args) )
