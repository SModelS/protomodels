#!/usr/bin/env python3

""" 
.. module:: hiscoreCLI
   :synposis: a command line interface to hiscore lists. to browser, interact,
    experiment.

.. moduleauthor:: Wolfgang Waltenberger <wolfgang.waltenberger@gmail.com>

"""

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
    from colorama import Fore as ansi
    print ( "[hiscoreCLI] starting interactive session." )
    from smodels.tools.physicsUnits import pb, fb, GeV, TeV
    print ( f"[hiscoreCLI]      Constants: {ansi.RED}pb, fb, GeV, TeV{ansi.RESET}" )
    from walker.hiscore import Hiscore
    from ptools.hiscoreTools import createHiscoreObj
    hi = createHiscoreObj ( args.infile, None, args.dbpath )
    # hi = Hiscore ( 0, False, args.infile )
    protomodels = hi.hiscores
    print ( f"[hiscoreCLI]      Variables: {ansi.RED}protomodels{ansi.RESET}" )
    import builder
    protomodel = protomodels
    from builder.protomodel import ProtoModel
    from builder.manipulator import Manipulator
    from tester.combiner import Combiner
    from tester.predictor import Predictor
    print ( f"[hiscoreCLI]         python: {ansi.RED}copy, numpy, scipy, scipy.stats, math{ansi.RESET}" )
    print ( f"[hiscoreCLI]        Modules: {ansi.RED}manipulator, hiscore, combiner, predictor, helpers{ansi.RESET}" )
    print ( f"[hiscoreCLI]        Classes: {ansi.RED}ProtoModel, Combiner, Predictor, Hiscore, Database{ansi.RESET}" )
    # so we can also use Andre's pcl files
    if type(protomodels)==ProtoModel:
        protomodel = protomodels
    else:
        protomodel = protomodels[0]
    ma = Manipulator ( protomodel )
    ma.M.createNewSLHAFileName()
    co = Combiner() # instantiate for convenience
    pr = Predictor( 0, do_combine=args.do_combine, dbpath=args.dbpath ) # instantiate for convenience
    print ( f"[hiscoreCLI] Instantiations: {ansi.RED}ma, co, hi, pr{ansi.RESET}" )
    from tester import combiner
    from walker import hiscore
    from tester import predictor
    from tester.combiner import Combiner
    from tester.predictor import Predictor
    from builder.protomodel import ProtoModel
    from smodels.experiment.databaseObj import Database
    import copy, numpy, scipy, scipy.stats, math
    from ptools import helpers
    # import hiscore #Keep it for convenience

    if args.execute not in [ "", None ]:
        if os.path.exists ( args.execute ):
            with open ( args.execute, "rt" ) as f:
                exec ( f.read() )

    if not args.nointeractive:
        import IPython
        IPython.embed( using=False )
    ma.M.delCurrentSLHA()
