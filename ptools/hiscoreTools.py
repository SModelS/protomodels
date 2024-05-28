#!/usr/bin/env python3

""" A module that contains all the functions for meddling with hiscores
"""

__all__ = [ "hiscoreHiNeedsUpdate", "fetchHiscoresObj" ]

import pickle, subprocess, sys, os, time
import numpy as np
from colorama import Fore as ansi
from scipy import stats
sys.path.insert(0,"../../")
from protomodels.csetup import setup
setup()
from builder.manipulator import Manipulator
from builder.protomodel import ProtoModel
from walker.hiscores import Hiscores
from typing import Union, Dict, List, Set
from argparse import Namespace
from os import PathLike

def count ( protomodels : List[ProtoModel] ) -> int:
    return len(protomodels)-protomodels.count(None)

def sortByTL ( protomodels : List[ProtoModel], n : int = 5 ) -> List:
    protomodels.sort ( reverse=True, key = lambda x: x.TL )
    return protomodels[:n] ## only n

def sortByK ( protomodels : List[ProtoModel], n : int = 5 ) -> List:
    protomodels.sort ( reverse=True, key = lambda x: x.K )
    return protomodels[:n] ## only n

def discuss ( protomodel, name ):
    print ( f"Currently {name:7s} K={protomodel.K:.3f}, TL={protomodel.TL:.3f} [{len(protomodel.unFrozenParticles())}/{len(protomodel.masses.keys())} particles, {len(protomodel.bestCombo)} predictions] (walker #{protomodel.walkerid})" )

def discussBest ( protomodel, detailed ):
    """ a detailed discussion of number 1 """
    p = 2. * ( 1. - stats.norm.cdf ( protomodel.TL ) ) ## two times because one-sided
    print ( "Current      best K=%.3f, TL=%.3f, p=%.2g [%d/%d particles, %d predictions] (walker #%d)" % \
            ( protomodel.K, protomodel.TL, p, len(protomodel.unFrozenParticles()),len(protomodel.masses.keys()),len(protomodel.bestCombo), protomodel.walkerid ) )
    if detailed:
        print ( f"Solution was found in step #{protomodel.step}" )
        for i in protomodel.bestCombo:
            print ( f"  prediction in best combo: {i.analysisId()} ({i.dataType})" )

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
        return f"{protomodel.nevents/1000}K evts"
    return str(protomodel.nevents)+ " evts"

def obtainHiscore ( number : int,
        hiscorefile : PathLike = "hiscores.dict" ) -> ProtoModel:
    """ obtain hiscore number <number> from >hiscorefile>

    :returns: protomodel object
    """
    hi = fetchHiscoresObj ( hiscorefile )
    TL = hi.hiscores[number].TL
    K = hi.hiscores[number].K
    print ( f"[hiscoreTools] obtaining #{number}: K={K:.3f}" )
    ret = hi.hiscores[ number ]
    return ret

def hiscoreHiNeedsUpdate ( dictfile : str = "hiscores.dict",
                           picklefile : str = "hiscores.cache",
                           entrynr : Union[None,int] = 0 ) -> bool:
    """ is hiscores.cache behind hiscores.dict, so it needs an update?
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
        if type(dictcontent)==dict: # make it work with single models also
            dictcontent = [ dictcontent ]
        if type(entrynr) != type(None) and entrynr >= len(dictcontent):
            print ( f"[hiscoreTools] entry #{entrynr} not existent" )
            return False
        f.close()
    from walker.hiscores import Hiscores
    hi = Hiscores ( 0, False, picklefile )
    def compare ( dentry, pentry ) -> bool:
        ## compare one dictentry with one pickleentry,
        ## true, if things are different
        newV = dentry["K"] + dentry["TL"] + sum(dentry["masses"].values()) + \
               sum(dentry["ssmultipliers"].values())
        if pentry == None:
            return True

        oldV = pentry.K + pentry.TL + sum(pentry.masses.values()) + \
               sum(pentry.ssmultipliers.values())
        if 2. * abs( newV - oldV ) / ( newV + oldV ) > 1e-4:
            # print ( f"[hiscoreTools] top V value changed {newV:.3f}..{oldV:.3f}" )
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

def fetchHiscoresObj ( dictfile : str = "hiscores.dict",
                       picklefile : Union[None,str] = None,
                       dbpath : str = "official" ) -> Hiscores:
    """ create Hiscores object from hiscores.cache file.
    update hiscores.cache file before, if needed.

    :param dictfile: dictionary to update hiscores.cache from, if needed.
    :param picklefile: the cached hiscore pickle file. if None,
    then dictfile but replace hiscores.dict with hiscores.cache

    :returns: hiscore object
    """
    if picklefile is None:
        picklefile = dictfile.replace(".dict",".cache" )
    from ptools import helpers
    shortname = helpers.simplifyUnixPath ( picklefile )
    if not hiscoreHiNeedsUpdate ( dictfile, picklefile, 0 ):
        print ( f"[hiscoreTools] can reuse cache: {shortname}" )
        return Hiscores ( 0, False, picklefile )
    print ( f"[hiscoreTools] updating cache: {shortname}" )
    hi = Hiscores.fromDictionaryFile ( dictfile, dbpath=dbpath )
    hi.writeListToPickle ( picklefile )
    return hi

def mergeTwoModels ( model1 : str, model2: str ) -> Union[None,Dict]:
    """ merge two models, add all particles from both models.
    If a particle appears in both models, its mass will be the average
    of the two, etc.

    :returns: merged model
    """
    if not os.path.exists ( model1 ):
        print ( f"[hiscoreTools] {model1} does not exist" )
        return None
    if not os.path.exists ( model2 ):
        print ( f"[hiscoreTools] {model2} does not exist" )
        return None
    f=open ( model1, "rt" )
    txt=f.read()
    f.close()
    dict1 = eval(txt)
    f=open ( model2, "rt" )
    txt=f.read()
    f.close()
    dict2 = eval(txt)
    import copy
    ret = copy.deepcopy(dict1)
    for pid,m in dict2["masses"].items():
        if not pid in ret["masses"]:
            # ok, we just copy <pid> from dict2
            ret["masses"][pid]=m
            ret["decays"][pid]=dict2["decays"][pid]
            for pidpair,ssms in dict2["ssmultipliers"].items():
                if pid in pidpair or -pid in pidpair:
                    ret["ssmultipliers"][pidpair]=ssms
        else:
            ret["masses"][pid]=.5*m + .5*ret["masses"][pid]
    for drop in [ "K", "TL", "xsecs[fb]", "walkerid", "step" ]:
        if drop in ret:
            ret.pop( drop )
    ret["timestamp"] = time.asctime()
    return ret

def mergeNModels ( models : List[Dict] ) -> Union[None,Dict]:
    """ merge two models, add all particles from both models.
    If a particle appears in both models, its mass will be the average
    of the two, etc.

    :returns: merged model
    """
    if len(models)== 0:
        print ( f"[hiscoreTools] provided empty list of models" )
        return None
    if len(models)==1: # trivial merge
        models[0]["timestamp"] = time.asctime()
        return models[0]
    ret = { "masses": {}, "decays": {}, "ssms": {} }

    def collectPids ( models : List[Dict] ) -> Set:
        """ collect all pids in all models """
        pids = set()
        for m in models:
            for k in m["masses"].keys():
                pids.add ( k )
        return pids
    def computeAverageMassesForPid ( pid : int, models : List[Dict] ) -> Dict:
        masses=[]
        for model in models:
            if pid in model["masses"]:
                masses.append ( model["masses"][pid] )
        return np.mean ( masses )

    def computeAverageDecaysForPid ( pid : int, models : List[Dict] ) -> Dict:
        decays = {}
        Stot = 0.
        for model in models:
            if pid in model["decays"]:
                mdecays = model["decays"][pid]
                for daughterpids, br in mdecays.items():
                    if not daughterpids in decays:
                        decays[daughterpids] = 0.
                    decays[daughterpids] += br
                    Stot += br
        for daughterpids, br in decays.items():
            decays[daughterpids]/=Stot
        return decays

    def computeAverageSSMs ( models : List[Dict] ) -> Dict:
        ssms = {}
        for model in models:
            mssms = model["ssmultipliers"]
            for k,v in mssms.items():
                ssms[k]=v # just bluntly take them over.
        return ssms

    # first, we collect all pids
    pids = collectPids ( models )
    # next we compute average masses
    masses, decays, ssms = {}, {}, {}
    for pid in pids:
        mass = computeAverageMassesForPid ( pid, models )
        masses[pid]=mass
        piddecays = computeAverageDecaysForPid ( pid, models )
        decays[pid] = piddecays
    ssms = computeAverageSSMs ( models )
    ret["masses"]=masses
    ret["decays"]=decays
    ret["ssmultipliers"]=ssms
    ret["timestamp"] = time.asctime()
    return ret

