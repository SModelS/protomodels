#!/usr/bin/env python3

""" A module that contains all the functions for meddling with hiscores
"""

__all__ = [ "hiscoreHiNeedsUpdate", "fetchHiscoresObj" ]

import pickle, subprocess, sys, os, time
import numpy as np
from colorama import Fore as ansi
from scipy import stats
sys.path.insert(0,"../")
sys.path.insert(0,"../../")
from protomodels.csetup import setup
setup()
from builder.manipulator import Manipulator
from builder.protomodel import ProtoModel
from walker.hiscores import Hiscores
from typing import Union, Dict, List, Set, Tuple
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
    print ( f"Current      best K={protomodel.K:.3f}, TL={protomodel.TL:.3f}, p={p:.2g} [{len(protomodel.unFrozenParticles())}/{len(protomodel.masses.keys())} particles, {len(protomodel.bestCombo)} predictions] (walker #{int(protomodel.walkerid)})" )
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
        sc = f"{int(c + 1)}th"
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
    return f"{protomodel.nevents!s} evts"

def obtainHiscore ( number : int,
        hiscorefile : PathLike = "hiscores_global.dict",
        walkerid : Union[str,int] = 0,
        dbpath : PathLike = "official" ) -> ProtoModel:
    """ obtain hiscore number <number> from <hiscorefile>

    :param walkerid: log everything as walker #walkerid

    :returns: protomodel object
    """
    hi = fetchHiscoresObj ( hiscorefile, walkerid = walkerid,
           dbpath = dbpath )
    TL = hi.hiscores[number].TL
    K = hi.hiscores[number].K
    sK = "K=None" if K==None else f"K={K:.3f}"
    print ( f"[hiscoreTools] obtaining #{number}: {sK}" )
    ret = hi.hiscores[ number ]
    return ret

def hiscoreHiNeedsUpdate ( dictfile : str = "hiscores_global.dict",
                           picklefile : str = "hiscores_global.cache",
                           entrynr : Union[None,int] = 0,
                           walkerid : Union[str,int] = 0 ) -> bool:
    """ is hiscores_global.cache behind hiscores_global.dict, so it needs an update?
    :param entrynr: check for for this entry, 0 is first.
    If None, check all
    :param walkerid: log everything as walker #walkerid

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
        txt = txt.replace ( "null","float('nan')" )
        dictcontent = eval(txt)
        if type(dictcontent)==dict: # make it work with single models also
            dictcontent = [ dictcontent ]
        if type(entrynr) != type(None) and entrynr >= len(dictcontent):
            print ( f"[hiscoreTools] entry #{entrynr} not existent" )
            return False
        f.close()
    from walker.hiscores import Hiscores
    hi = Hiscores ( walkerid, False, picklefile )

    def compare ( dentry, pentry ) -> Tuple[bool,str]:
        """ compare dentry with pentry

        :returns: Tuple[bool,str] bool is true, if different!
        str gives explanation
        """
        ## compare one dictentry with one pickleentry,
        ## true, if things are different
        if pentry == None:
            return True
        if pentry.K == None: ## picklefile is not working
            # so, update!
            return True, "pentry.K is None"
        if not "K" in dentry or not "TL" in dentry:
            return True, "no K or TL in dentry"
        newV = sum(dentry["masses"].values()) + \
               sum(dentry["ssmultipliers"].values())

        oldV = sum(pentry.masses.values()) + \
               sum(pentry.ssmultipliers.values())

        if "K" in dentry and dentry["K"] is not None:
            newV += dentry["K"]
            oldV += pentry.K
        if "TL" in dentry and dentry["TL"] is not None:
            newV += dentry["TL"]
            oldV += pentry.TL
        delta = 2. * abs( newV - oldV ) / ( newV + oldV ) 
        if delta > 1e-4:
            # print ( f"[hiscoreTools] top V value changed {newV:.3f}..{oldV:.3f}" )
            return True, "Vs are different by {delta}"
        return False, "entries are the same"

    if type(entrynr) == int:
        if entrynr >= len(hi.hiscores):
            return True
        dictentry = dictcontent[entrynr]
        pickleentry = hi.hiscores[entrynr]
        ret, reason = compare ( dictentry, pickleentry )
        return ret
    # entrynr is None, check all!
    if len(dictcontent)!=len(hi.hiscores):
        return True
    for i in range(len(dictentry)):
        ci = compare ( dictcontent[i], hi.hiscores[i] )
        if ci:
            return True
    return False

def fetchHiscoresObj ( dictfile : str = "hiscores_global.dict",
                       picklefile : Union[None,str] = None,
                       dbpath : str = "official",
                       walkerid : Union[str,int] = 0 ) -> Hiscores:
    """ create Hiscores object from hiscores_global.cache file.
    update hiscores_global.cache file before, if needed.

    :param dictfile: dictionary to update hiscores_global.cache from, if needed.
    :param picklefile: the cached hiscore pickle file. if None,
    then dictfile but replace hiscores_global.dict with hiscores_global.cache
    :param walkerid: log everything as walker #walkerid

    :returns: hiscore object
    """
    if picklefile is None:
        picklefile = dictfile.replace(".dict",".cache" )
    from ptools import helpers
    shortname = helpers.simplifyUnixPath ( picklefile )
    if not hiscoreHiNeedsUpdate ( dictfile, picklefile, walkerid=walkerid ):
        print ( f"[hiscoreTools] can reuse cache: {shortname}" )
        return Hiscores ( 0, False, picklefile )
    print ( f"[hiscoreTools] updating cache: {shortname} ... " )
    hi = Hiscores.fromDictionaryFile ( dictfile, dbpath=dbpath, walkerid = walkerid )
    hi.writeListToPickle ( picklefile )
    print ( f"[hiscoreTools] cache {shortname} updated!" )
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

