#!/usr/bin/env python3

"""
.. module:: hiscoreCLI
   :synposis: a command line interface to hiscore lists. to browser, interact,
    experiment.

.. moduleauthor:: Wolfgang Waltenberger <wolfgang.waltenberger@gmail.com>

"""

import os, sys
pmpath =os.path.abspath ( f"{os.path.dirname (  __file__ )}/../../" )
sys.path.insert(0,pmpath)

from typing import Union
from smodels_utils.helper.terminalcolors import LIGHTGREEN, RESET, CYAN

def pprint ( *args ):
    """ logging """
    print ( f"{CYAN}[hiscoreCLI]{RESET} {' '.join(map(str,args))}" )

def cli( infile : str = "hiscores_global.dict",
         dbpath : str = "official", do_srcombine : bool= True,
         walkerid : Union[str,int] = 0 ):
    """ fire up the interactive shell, preconfigured!

    :param infile: read hiscore from infile, can contain a single hiscore as a
    dictionary, or a list of hiscores.
    :param dbpath: path to database
    :param do_srcombine: if true, do sr-combinations (when would you not?)
    :param walkerid: log as walker #walkerid

    Example usage (at the interactive shell)

    .. code-block:: python3

    >>> # change a mass, see how the prediction changes
    >>> ma.M.masses[1000024] = 70.
    >>> pr.predict ( ma.M, keep_predictions = True )
    """
    import sys
    sys.path.insert(0,"../")
    sys.path.insert(0,"../../")
    import csetup
    csetup.setup()
    from smodels_utils.helper.terminalcolors import RED, GREEN, YELLOW
    pprint ( "starting interactive session." )
    import copy, numpy, scipy, scipy.stats
    pprint ( f"        python: {RED}copy, numpy, scipy, scipy.stats, math{RESET}" )
    from smodels.base.physicsUnits import pb, fb, GeV, TeV
    pprint ( f"     Constants: {RED}pb, fb, GeV, TeV{RESET}" )
    from ptools.hiscoreTools import fetchHiscoresObj
    from builder import manipulator
    from walker import hiscores
    from tester import combiner, predictor
    from ptools import helpers
    pprint ( f"       Modules: {RED}manipulator, hiscores, combiner, predictor, helpers{RESET}" )
    from walker.hiscores import Hiscores
    from base.runEnviron import RunEnviron
    from builder.protomodel import ProtoModel
    from builder.manipulator import Manipulator
    from tester.combiner import Combiner
    from tester.predictor import Predictor
    from tester.critic import Critic
    from ptools.sparticleNames import SParticleNames
    from smodels.experiment.databaseObj import Database
    environ = RunEnviron()
    pprint ( f"       Classes: {RED}ProtoModel, Combiner, Predictor, Hiscores, Database,{RESET}" )
    pprint ( f"                {RED}SParticleNames{RESET}" )
    hi = fetchHiscoresObj ( infile, None, environ = environ, walkerid = walkerid )
    pprint ( f"{RED}hi = fetchHiscoresObj ('{infile}', ... ) # Hiscore {RESET}" )
    namer = SParticleNames()
    from importlib import reload
    pprint ( f"{RED}namer = SParticleNames(){RESET}" )
    protomodel = hi.hiscores[0]
    protomodel.walkerid = walkerid
    pprint ( f"{RED}protomodel = hi.hiscores[0]{RESET}" )
    ma = Manipulator ( protomodel, environ )
    pprint ( f"{RED}ma = Manipulator ( protomodel ){RESET}" )
    ma.M.createNewSLHAFileName()
    pprint ( f"{RED}co = Combiner ( protomodel ){RESET}" )
    co = Combiner( walkerid ) # instantiate for convenience
    pprint ( f"{RED}pr = Predictor ( ){RESET}" )
    pr = Predictor( walkerid, environ = environ ) # instantiate for convenience
    pprint ( f"{RED}cr = Critic ( ){RESET}" )
    cr = Critic ( walkerid, environ = environ )
    cr = Critic ( walkerid, environ = environ )
    pprint ( f"{YELLOW}pr.predict(ma,keep_predictions=True,force_computation_K=False,{RESET}" )
    pprint ( f"{YELLOW}           keep_slhafile=True{RESET}" )
    pr.predict( ma, keep_predictions=True, force_computation_K=False,
                keep_slhafile=True )

    # print ( f"[hiscoreCLI] Instantiations: {RED}ma, co, hi, pr{RESET}" )

    if args.execute not in [ "", None ]:
        if os.path.exists ( args.execute ):
            with open ( args.execute, "rt" ) as f:
                print ( f"[hiscoreCLI] {GREEN}executing {args.execute}{RESET}" )
                globals().update(locals())
                exec ( f.read(), globals() )

    if not args.nointeractive:
        import IPython
        IPython.embed( using=False )
    ma.M.delCurrentSLHA()

if __name__ == "__main__":
    import argparse, os, sys
    argparser = argparse.ArgumentParser(
            description='interactive session with hiscore list loaded' )
    argparser.add_argument ( '-f', '--infile',
            help='Hiscore file. [hiscores_global.dict]',
            type=str, default="hiscores_global.dict" )
    argparser.add_argument ( '-i', '--walkerid',
            help="walker id ['cli']",
            type=str, default='cli' )
    argparser.add_argument ( '-d', '--dbpath',
            help='Database path. If auto, get from run.dict. [auto]',
            type=str, default=None )
    argparser.add_argument ( '-n', '--nointeractive',
            help='Dont start interactive shell',
            action = "store_true" )
    argparser.add_argument ( '-D', '--dont_srcombine',
            help='Do NOT combine results',
            action = "store_true" )
    argparser.add_argument ( '-x', '--execute',
            help='execute python script EXECUTE before going interactive [None]',
            type=str, default=None )
    args = argparser.parse_args()
    do_srcombine = True
    if args.dont_srcombine:
        do_srcombine = False
    #else:
    #    if args.do_srcombine == False:
    #        print ( f"[hiscoreCLI] really? no srcombine? will anyhow set to true" )
    #        args.do_srcombine = True
    if not os.path.exists ( args.infile ):
        pprint ( f"error: input file {args.infile} does not exist." )
        sys.exit()
    if os.path.exists ( "./run.dict" ):
        pprint ( f"found run.dict file. will use its values." )
        with open ( "./run.dict", "rt" ) as f:
            txt = f.read()
            f.close()
            d = eval(txt)
            if "do_srcombine" in d:
                do_srcombine = d["do_srcombine"]
            if args.dbpath in [ "auto", None ] and "dbpath" in d:
                pprint ( f"setting dbpath to {d['dbpath']}" )
                args.dbpath = d['dbpath']
    cli ( args.infile, args.dbpath, do_srcombine, args.walkerid )
