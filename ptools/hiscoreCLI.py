#!/usr/bin/env python3

""" 
.. module:: hiscoreCLI
   :synposis: a command line interface to hiscore lists. to browser, interact,
    experiment.

.. moduleauthor:: Wolfgang Waltenberger <wolfgang.waltenberger@gmail.com>

"""

def cli( infile : str = "hiscores.dict", 
         dbpath : str = "official", do_combine : bool= True ):
    """ fire up the interactive shell, preconfigured!

    Example usage (at the interactive shell)

    .. code-block:: python3

    >>> # change a mass, see how the prediction changes
    >>> ma.M.masses[1000024] = 70.
    >>> pr.predict ( ma.M, keep_predictions = True )
    """
    from colorama import Fore as ansi
    print ( "[hiscoreCLI] starting interactive session." )
    import copy, numpy, scipy, scipy.stats
    print ( f"[hiscoreCLI]         python: {ansi.RED}copy, numpy, scipy, scipy.stats, math{ansi.RESET}" )
    from smodels.tools.physicsUnits import pb, fb, GeV, TeV
    print ( f"[hiscoreCLI]      Constants: {ansi.RED}pb, fb, GeV, TeV{ansi.RESET}" )
    from ptools.hiscoreTools import createHiscoreObj
    from builder import manipulator
    from walker import hiscore
    from tester import combiner, predictor
    from ptools import helpers
    print ( f"[hiscoreCLI]        Modules: {ansi.RED}manipulator, hiscore, combiner, predictor, helpers{ansi.RESET}" )
    from walker.hiscore import Hiscore
    from builder.protomodel import ProtoModel
    from builder.manipulator import Manipulator
    from tester.combiner import Combiner
    from tester.predictor import Predictor
    from smodels.experiment.databaseObj import Database
    print ( f"[hiscoreCLI]        Classes: {ansi.RED}ProtoModel, Combiner, Predictor, Hiscore, Database{ansi.RESET}" )
    hi = createHiscoreObj ( infile, None, dbpath )
    print ( f"[hiscoreCLI] {ansi.RED}hi = createHiscoreObj ('{infile}', ... ) # Hiscore {ansi.RESET}" )
    #print ( f"[hiscoreCLI]      Variables: {ansi.RED}protomodel{ansi.RESET}" )
    protomodel = hi.hiscores[0]
    print ( f"[hiscoreCLI] {ansi.RED}protomodel = hi.hiscores[0]{ansi.RESET}" )
    ma = Manipulator ( protomodel )
    print ( f"[hiscoreCLI] {ansi.RED}ma = Manipulator ( protomodel ){ansi.RESET}" )
    ma.M.createNewSLHAFileName()
    print ( f"[hiscoreCLI] {ansi.RED}co = Combiner ( protomodel ){ansi.RESET}" )
    co = Combiner() # instantiate for convenience
    print ( f"[hiscoreCLI] {ansi.RED}pr = Predictor ( ){ansi.RESET}" )
    pr = Predictor( 0, do_combine=do_combine, dbpath=dbpath ) # instantiate for convenience
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
    argparser.add_argument ( '-D', '--do_combine',
            help='Do combine results',
            action = "store_true" )
    argparser.add_argument ( '-x', '--execute',
            help='execute python script EXECUTE before going interactive [None]',
            type=str, default=None )
    args = argparser.parse_args()
    if not os.path.exists ( args.infile ):
        print ( f"[hiscoreCLI] error: input file {args.infile} does not exist." )
        sys.exit()
    if os.path.exists ( "./run.dict" ):
        print ( f"[hiscoreCLI] found run.dict file. will use its values." )
        with open ( "./run.dict", "rt" ) as f:
            txt = f.read()
            f.close()
            d = eval(txt)
            if "do_combine" in d:
                args.do_combine = d["do_combine"]
    cli ( args.infile, args.dbpath, args.do_combine )
