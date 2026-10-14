#!/usr/bin/env python3

import os, sys

""" module to contain helper methods for the multiverse
"""


def checkSignalModel ( signal_model : str ):
    if signal_model == "":
        print ( f"[logZs] no signal model defined." )
        sys.exit()
    if not os.path.exists ( signal_model ):
        print ( f"[logZs] {signal_model} does not exist." )
        sys.exit()

def obtainTestStatForPred ( pred : "TheoryPrediction" ) -> dict:
    """ for the given theory prediction, get test stats
    :returns: dictionary with all test statistics
    """
    from ptools.helpers import zFromNLLs
    nll_sm = pred.nll ( mu= 0. )
    nll_bsm = pred.nll ( mu= 1. )
    Z =  zFromNLLs ( nll_sm, nll_bsm )
    ret = { "nll_sm": nll_sm, "nll_bsm": nll_bsm, "Z": Z }
    return ret

def obtainTestStatsForAllPreds ( preds : list ) -> dict:
    """ obtain test stats for all theory predictions
    :returns: dictionary with dataset id as key
    """
    ret = {}
    for pred in preds:
        dId = pred.dataId()
        ts = obtainTestStatForPred ( pred )
        ret[dId]=ts
    return ret

def writeDict ( D : dict, outfile : str ):
    from smodels_utils.helper.various import py_dumps
    with open ( outfile, "wt" ) as f:
        d = py_dumps (  D )
        f.write ( f"{d}\n" )
    f.close()

def obtainDictForDM ( dm : float, ma : "Manipulator", pr : "predictor",
       origmasses : dict, environ : "RunEnviron" ) -> dict:
    """ obtain our comprehensive dictionary for a specific mass change
    """
    environ.database.subs[0].setActiveExpResults( analysisIDs=[ "CMS-EXO-20-004" ] )
    ma.M.masses[1000022] = origmasses[1000022]+dm
    ma.M.masses[1000023] = origmasses[1000023]-dm
    predictions = pr.obtainPredictions ( ma.M )
    D = obtainTestStatsForAllPreds ( predictions )
    return D

def logZs ( signal_model : str = "rundir2_stable.dict",
        outfile : str = "Zs.dict",
        interactive : bool = False ) -> dict:
    """ create the truth.dict file
    :param signal_model: the input signal_model.dict file

    :returns: the truth as a dictionary
    """
    checkSignalModel ( signal_model )
    from builder.manipulator import Manipulator
    # from builder.protomodel import ProtoModel
    from tester.predictor import Predictor
    from pbase.runEnviron import RunEnviron
    environ = RunEnviron ( "run_Zs.dict" )

    walkerid = "logz"
    ma = Manipulator( signal_model, walkerid = walkerid, environ = environ )
    pr = Predictor( walkerid, environ )
    print ( f"[logZs] starting" )
    import copy
    origmasses = copy.deepcopy( ma.M.masses )
    ret = {}
    for dm in range(-10,30):
        D = obtainDictForDM ( dm, ma, pr, origmasses, environ )
        ret[dm]=D
    writeDict ( ret, outfile )

    if interactive:
        import sys, IPython; IPython.embed( colors = "neutral" ); sys.exit()

if __name__ == "__main__":
    import argparse
    argparser = argparse.ArgumentParser(
        description='log the Z values' )
    argparser.add_argument ( '-i', '--interactive', action="store_true",
        help="interactive" )
    args = argparser.parse_args()
    logZs( **vars(args) )
