#!/usr/bin/env python3

""" A class that centralizes access to the hiscore list over multiple threads.
"""

import pickle, subprocess, sys, os
from colorama import Fore as ansi
from scipy import stats
sys.path.insert(0,"../")
from ptools.csetup import setup
setup()
from builder.manipulator import Manipulator
from walker.hiscore import Hiscore
from typing import Union, Dict
from argparse import Namespace

def count ( protomodels ):
    return len(protomodels)-protomodels.count(None)

def sortByZ ( protomodels ):
    protomodels.sort ( reverse=True, key = lambda x: x.Z )
    return protomodels[:5] ## only 5

def sortByK ( protomodels ):
    protomodels.sort ( reverse=True, key = lambda x: x.K )
    return protomodels[:5] ## only 5

def storeList ( protomodels, savefile ):
    """ store the best protomodels in another hiscore file """
    h = Hiscore ( 0, True, savefile, backup=True, hiscores = protomodels )
    h.hiscores = protomodels
    print ( "[hiscore] saving %d protomodels to %s" % \
            ( count(protomodels), savefile ) )
    if savefile.endswith ( ".pcl" ) or savefile.endswith( ".hi" ):
        h.writeListToPickle ( check=False )
        if "states" in savefile: ## do both for the states
            h.writeListToDictFile()
    else: ## assume a dict file
        h.writeListToDictFile()

def discuss ( protomodel, name ):
    print ( "Currently %7s K=%.3f, Z=%.3f [%d/%d particles, %d predictions] (walker #%d)" % \
            (name, protomodel.K, protomodel.Z, len(protomodel.unFrozenParticles()),len(protomodel.masses.keys()),len(protomodel.bestCombo), protomodel.walkerid ) )

def discussBest ( protomodel, detailed ):
    """ a detailed discussion of number 1 """
    p = 2. * ( 1. - stats.norm.cdf ( protomodel.Z ) ) ## two times because one-sided
    print ( "Current      best K=%.3f, Z=%.3f, p=%.2g [%d/%d particles, %d predictions] (walker #%d)" % \
            ( protomodel.K, protomodel.Z, p, len(protomodel.unFrozenParticles()),len(protomodel.masses.keys()),len(protomodel.bestCombo), protomodel.walkerid ) )
    if detailed:
        print ( "Solution was found in step #%d" % protomodel.step )
        for i in protomodel.bestCombo:
            print ( "  prediction in best combo: %s (%s)" % ( i.analysisId(), i.dataType() ) )

def printProtoModels ( protomodels, detailed, nmax=10 ):
    names = { 0: "highest", 1: "second", 2: "third" }
    for c,protomodel in enumerate(protomodels):
        if c >= nmax:
            break
        if protomodel == None:
            break
        sc = "%dth" % (c+1)
        if c in names.keys():
            sc = names[c]
        if c==0:
            discussBest ( protomodel, detailed )
        else:
            discuss ( protomodel, sc )

def pprintEvs ( protomodel ):
    """ pretty print number of events """
    if protomodel.nevents > 1000:
        return "%dK evts" % ( protomodel.nevents/1000 )
    return str(protomodel.nevents)+ " evts"

def compileList( nmax ):
    """ compile the list from individual H*hi
    """
    import glob
    files = glob.glob ( "H*.hi" )
    allprotomodels=[]
    import progressbar
    pb = progressbar.ProgressBar(widgets=["file #",progressbar.Counter(),
            "/%d " % len(files), progressbar.Percentage(),
            progressbar.Bar( marker=progressbar.RotatingMarker() ),
            progressbar.AdaptiveETA()])
    pb.maxval = len(files)
    pb.start()
    for ctr,fname in enumerate(files):
        pb.update(ctr)
        try:
            with open( fname,"rb+") as f:
                protomodels = pickle.load ( f )
                try:
                    pickle.load(f)
                except EOFError:
                    pass
                ## add protomodels, but without the Nones
                f.close()
                allprotomodels += list ( filter ( None.__ne__, protomodels ) )
                allprotomodels = sortByK ( allprotomodels )
        except ( AttributeError, Exception, IOError, OSError, FileNotFoundError, EOFError, UnicodeDecodeError, pickle.UnpicklingError ) as e:
            cmd = "rm -f %s" % fname
            print ( "[hiscore] could not open %s (%s). %s." % ( fname, e, cmd ) )
            subprocess.getoutput ( cmd )
    pb.finish()
    if nmax > 0:
        while len(allprotomodels)<nmax:
            allprotomodels.append ( None )
    return allprotomodels

def updateHiscoreHi ( args : Namespace ) -> Dict:
    """ the function that updates the hiscore.hi file
    :param args: detailed, outfile, infile, print, fetch, nmax,
                 check, interactive, nevents.
                 see "if __main__" part below.
    :returns: { "Z": highest significance,
                "step": step, "model": model, "K": bayesian_K  }
    """
    ret =  { "Z": 0., "step": 0, "model": None, "K": -100. }

    if args.detailed:
        args.print = True
    if args.outfile.lower() in [ "none", "", "false" ]:
        args.outfile = None
    infile = args.infile
    if type(infile) is str and infile.lower() in [ "none", "" ]:
        infile = None
    trundir = None
    if hasattr ( args, "rundir" ):
        trundir = args.rundir
    rundir = setup( trundir )
    if infile == "default":
        infile = f"{rundir}/hiscore.hi"
    if args.outfile == infile:
        print ( "[hiscore] outputfile is same as input file. will assume that you do not want me to write out at all." )
        args.outfile = None

    if infile is None:
        print ( f"[hiscore] compiling a hiscore list with {args.nmax} protomodels" )
        protomodels = compileList( args.nmax ) ## compile list from H<n>.hi files
    else:
        if not os.path.exists ( infile ):
            if os.path.exists ( "hiscores.dict" ):
                print ( f"[hiscoreTools] {infile} does not exist, but hiscores.dict does! rebuild {infile}." )
                hi = Hiscore.fromDictionaryFile ( "hiscores.dict", dbpath=args.dbpath )
                hi.writeListToPickle ( infile )
            else:
                print ( f"[hiscoreTools] neither {infile} nor hiscores.dict exists. abort!" )
                sys.exit()

        with open(infile,"rb") as f:
            try:
                protomodels = pickle.load ( f )
                try:
                    pickle.load ( f )
                except EOFError:
                    pass
                f.close()
            except (BlockingIOError,OSError) as e:
                print ( "file handling error on %s: %s" % ( infile, e ) )
                ## make sure we dont block!
                raise e

    if protomodels[0] == None:
        print ( "[hiscore] error, we have an empty hiscore list" )
        return ret

    sin = infile
    if sin == None:
        sin = "H*.hi"
    pevs = pprintEvs ( protomodels[0] )
    print ( "[hiscore] hiscore from %s[%d] is at K=%.3f, Z=%.3f (%s)" % \
            ( sin, protomodels[0].walkerid, protomodels[0].K, protomodels[0].Z, pevs ) )

    if type(args.outfile)==str and (".pcl" in args.outfile or ".hi" in args.outfile ):
        if not hasattr ( protomodels[0], "analysisContributions" ):
            print ( "[hiscore] why does the winner not have analysis contributions?" )
            ma = Manipulator ( protomodels[0] )
            hi = Hiscore( 0, True, f"{rundir}/hiscore.hi" )
            hi.computeAnalysisContributions(ma)
            protomodels[0]=ma.M
            hi.writeListToPickle()
        if not hasattr ( protomodels[0], "particleContributions" ):
            print ( "[hiscore] why does the winner not have particle contributions?" )
            ma = Manipulator ( protomodels[0] )
            from tester.predictor import Predictor
            predictor = None
            dbpath = args.dbpath
            if not "/" in dbpath:
                dbpath =  f"{rundir}/{args.dbpath}"
            if hasattr ( args, "dbpath" ):
                dbpath = args.dbpath
            if os.path.exists ( dbpath ):
                predictor = Predictor ( 0, dbpath = dbpath, expected = False, 
                        select= "all", do_combine = args.do_combine )
            hi = Hiscore( 0, True, f"{rundir}/hiscore.hi", predictor = predictor )
            hi.computeParticleContributions(ma)
            protomodels[0]=ma.M
            hi.writeListToPickle()

    if args.outfile is not None:
        storeList ( protomodels, args.outfile )

    if hasattr ( args, "check" ) and args.check:
        protomodel = protomodels[0]
        protomodel.predict()
        print ( "[hiscore] args.check, implement" )

    if args.print:
        printProtoModels ( protomodels, args.detailed, min ( 10, args.nmax ) )

    if len(protomodels)>0 and protomodels[0] != None:
        ret["Z"]=protomodels[0].Z
        ret["K"]=protomodels[0].K
        ret["step"]=protomodels[0].step
        ret["model"]=protomodels[0]
        return ret
    return ret

if __name__ == "__main__":
    import argparse
    argparser = argparse.ArgumentParser(
            description='interactive session with hiscore list loaded' )
    argparser.add_argument ( '-f', '--infile',
            help='Hiscore file. [hiscore.hi]',
            type=str, default="hiscore.hi" )
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
        print ( f"[hiscoreTools] error: input file {args.infile} does not exist." )
        sys.exit()
    if os.path.exists ( "./run.dict" ):
        print ( f"[hiscoreTools] found run.dict file. will use its values." )
        with open ( "./run.dict", "rt" ) as f:
            txt = f.read()
            f.close()
            d = eval(txt)
            if "do_combine" in d:
                args.do_combine = d["do_combine"]
    hi = Hiscore ( 0, False, args.infile )
    protomodels = hi.hiscores
    import builder
    protomodel = protomodels
    from builder.protomodel import ProtoModel
    # so we can also use Andre's pcl files
    if type(protomodels)==ProtoModel:
        protomodel = protomodels
    else:
        protomodel = protomodels[0]
    ma = Manipulator ( protomodel )
    ma.M.createNewSLHAFileName()
    from smodels.tools.physicsUnits import pb, fb, TeV
    print ( "[hiscoreTools] starting interactive session." )
    print ( f"[hiscoreTools]      Constants: {ansi.RED}pb, fb, TeV{ansi.RESET}" )
    print ( f"[hiscoreTools]      Variables: {ansi.RED}protomodels{ansi.RESET}" )
    print ( f"[hiscoreTools]         python: {ansi.RED}copy, numpy, scipy, scipy.stats, math{ansi.RESET}" )
    print ( f"[hiscoreTools]        Modules: {ansi.RED}manipulator, hiscore, combiner, predictor, helpers{ansi.RESET}" )
    print ( f"[hiscoreTools]        Classes: {ansi.RED}ProtoModel, Combiner, Predictor, Hiscore, Database{ansi.RESET}" )
    print ( f"[hiscoreTools] Instantiations: {ansi.RED}ma, co, hi, pr{ansi.RESET}" )
    from tester import combiner
    from walker import hiscore
    from tester import predictor
    from tester.combiner import Combiner
    from tester.predictor import Predictor
    from builder.protomodel import ProtoModel
    from smodels.base.physicsUnits import pb, fb, GeV, TeV
    from smodels.experiment.databaseObj import Database
    import copy, numpy, scipy, scipy.stats, math
    co = Combiner() # instantiate for convenience
    pr = Predictor( 0, do_combine=args.do_combine, dbpath=args.dbpath ) # instantiate for convenience
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
