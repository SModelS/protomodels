#!/usr/bin/env python3

""" A class that centralizes access to the hiscore list over multiple threads.
"""

__all__ = [ "updateHiscoreHi", "hiscoreHiNeedsUpdate", "createHiscoreObj" ]

import pickle, subprocess, sys, os
from colorama import Fore as ansi
from scipy import stats
sys.path.insert(0,"../../")
from protomodels.csetup import setup
setup()
from builder.manipulator import Manipulator
from builder.protomodel import ProtoModel
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

def hiscoreHiNeedsUpdate ( dictfile : str = "hiscores.dict", 
                           picklefile : str = "hiscore.hi",
                           entrynr : Union[None,int] = 0 ) -> bool:
    """ is hiscore.hi behind hiscores.dict, so it needs an update?
    :param entrynr: check for for this entry, 0 is first. 
    If None, check all

    :returns: true if update is needed
    """
    if not os.path.exists ( dictfile ):
        # we dont even have a dict file, nothing to update
        return False
    if not os.path.exists ( picklefile ):
        # simple case, we dont have a pickle file, lets create one!
        return True
    with open ( dictfile, "rt" ) as f:
        txt = f.read()
        dictcontent = eval(txt)
        if type(entrynr) != type(None) and entrynr >= len(dictcontent):
            print ( f"[hiscoreTools] entry #{entrynr} not existent" )
            return False
        f.close()
    from walker.hiscore import Hiscore
    hi = Hiscore ( 0, False, picklefile )
    def compare ( dentry, pentry ):
        ## compare one dictentry with one pickleentry
        newK = dentry["K"]
        oldK = pentry.K
        if ( newK - oldK ) > 1e-4:
            print ( f"[hiscoreTools] top K value changed {newK:.3f}..{oldK:.3f}" )
            return True
        return False

    if type(entrynr) == int:
        if entrynr >= len(hi.hiscores):
            return True
        dictentry = dictcontent[entrynr]
        pickleentry = hi.hiscores[entrynr]
        return compare ( dictentry, pickleentry )
    # entrynr is None, check all!
    if len(dictcontent)!=len(hi.hiscores):
        return True
    for i in range(len(dictentry)):
        ci = compare ( dictcontent[i], hi.hiscores[i] )
        if ci:
            return True
    return False

def createHiscoreObj ( dictfile : str = "hiscores.dict", 
                       picklefile : Union[None,str] = None,
                       dbpath : str = "official" ) -> Hiscore:
    """ create Hiscore object from hiscore.hi file. 
    update hiscore.hi file before, if needed.

    :param dictfile: dictionary to update hiscore.hi from, if needed.
    :param picklefile: the cached hiscore pickle file. if None,
    then dictfile but replace hiscores.dict with hiscore.hi

    :returns: hiscore object
    """
    if picklefile is None:
        picklefile = dictfile.replace("hiscores.dict","hiscore.hi" )
    if not hiscoreHiNeedsUpdate ( dictfile, picklefile, 0 ):
        print ( f"[hiscoreTools] can reuse cache: {picklefile}" )
        return Hiscore ( 0, False, picklefile )
    print ( f"[hiscoreTools] updating cache: {picklefile}" )
    hi = Hiscore.fromDictionaryFile ( dictfile, dbpath=dbpath )
    hi.writeListToPickle ( picklefile )
    return hi
        
def updateHiscoreHi ( args : Namespace ) -> Dict:
    """ the function that updates the hiscore.hi file
    :param args: detailed, outfile, infile, print, fetch, nmax,
                 check, interactive, nevents.
                 see "if __main__" part below.
    :returns: { "Z": highest significance,
                "step": step, "model": model, "K": bayesian_K  }
    """
    print ( "FIXME obsolete!!!" )
    # FIXME replace!
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
        infile = "hiscore.hi"
    dictfile = infile.replace( "hiscore.hi","hiscores.dict" )
    hi = createHiscoreObj ( dictfile )
    protomodels = hi.hiscores

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
