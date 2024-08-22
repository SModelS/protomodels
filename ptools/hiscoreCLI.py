#!/usr/bin/env python3

""" 
.. module:: hiscoreCLI
   :synposis: a command line interface to hiscore lists. to browser, interact,
    experiment.

.. moduleauthor:: Wolfgang Waltenberger <wolfgang.waltenberger@gmail.com>

"""

def cli( infile : str = "hiscores.dict", 
         dbpath : str = "official", do_srcombine : bool= True ):
    """ fire up the interactive shell, preconfigured!

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
    from colorama import Fore as ansi
    print ( "[hiscoreCLI] starting interactive session." )
    import copy, numpy, scipy, scipy.stats
    print ( f"[hiscoreCLI]         python: {ansi.RED}copy, numpy, scipy, scipy.stats, math{ansi.RESET}" )
    from smodels.base.physicsUnits import pb, fb, GeV, TeV
    print ( f"[hiscoreCLI]      Constants: {ansi.RED}pb, fb, GeV, TeV{ansi.RESET}" )
    from ptools.hiscoreTools import fetchHiscoresObj
    from builder import manipulator
    from walker import hiscores
    from tester import combiner, predictor
    from ptools import helpers
    print ( f"[hiscoreCLI]        Modules: {ansi.RED}manipulator, hiscores, combiner, predictor, helpers{ansi.RESET}" )
    from walker.hiscores import Hiscores
    from builder.protomodel import ProtoModel
    from builder.manipulator import Manipulator
    from tester.combiner import Combiner
    from tester.predictor import Predictor
    from tester.critic import Critic
    from ptools.sparticleNames import SParticleNames
    from smodels.experiment.databaseObj import Database
    print ( f"[hiscoreCLI]        Classes: {ansi.RED}ProtoModel, Combiner, Predictor, Hiscores, Database," )
    print ( f"                             SParticleNames{ansi.RESET}" )
    hi = fetchHiscoresObj ( infile, None, dbpath )
    print ( f"[hiscoreCLI] {ansi.RED}hi = fetchHiscoresObj ('{infile}', ... ) # Hiscore {ansi.RESET}" )
    namer = SParticleNames()
    from importlib import reload
    print ( f"[hiscoreCLI] {ansi.RED}namer = SParticleNames(){ansi.RESET}" )
    protomodel = hi.hiscores[0]
    print ( f"[hiscoreCLI] {ansi.RED}protomodel = hi.hiscores[0]{ansi.RESET}" )
    ma = Manipulator ( protomodel )
    print ( f"[hiscoreCLI] {ansi.RED}ma = Manipulator ( protomodel ){ansi.RESET}" )
    ma.M.createNewSLHAFileName()
    print ( f"[hiscoreCLI] {ansi.RED}co = Combiner ( protomodel ){ansi.RESET}" )
    co = Combiner() # instantiate for convenience
    print ( f"[hiscoreCLI] {ansi.RED}pr = Predictor ( ){ansi.RESET}" )
    pr = Predictor( 0, do_srcombine=do_srcombine, dbpath=dbpath ) # instantiate for convenience
    print ( f"[hiscoreCLI] {ansi.RED}cr = Critic ( ){ansi.RESET}" )
    cr = Critic ( 0, do_srcombine=do_srcombine, dbpath=dbpath )

    # print ( f"[hiscoreCLI] Instantiations: {ansi.RED}ma, co, hi, pr{ansi.RESET}" )

    if args.execute not in [ "", None ]:
        if os.path.exists ( args.execute ):
            with open ( args.execute, "rt" ) as f:
                exec ( f.read() )

    if not args.nointeractive:
        import IPython
        IPython.embed( using=False )
    ma.M.delCurrentSLHA()

if __name__ == "__main__":
    import argparse, os, sys
    argparser = argparse.ArgumentParser(
            description='interactive session with hiscore list loaded' )
    argparser.add_argument ( '-f', '--infile',
            help='Hiscore file. [hiscores.dict]',
            type=str, default="hiscores.dict" )
    argparser.add_argument ( '-d', '--dbpath',
            help='Database path. [official]',
            type=str, default="official" )
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
        print ( f"[hiscoreCLI] error: input file {args.infile} does not exist." )
        sys.exit()
    if os.path.exists ( "./run.dict" ):
        print ( f"[hiscoreCLI] found run.dict file. will use its values." )
        with open ( "./run.dict", "rt" ) as f:
            txt = f.read()
            f.close()
            d = eval(txt)
            if "do_srcombine" in d:
                do_srcombine = d["do_srcombine"]
    cli ( args.infile, args.dbpath, do_srcombine )
