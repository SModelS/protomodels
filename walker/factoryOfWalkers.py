#!/usr/bin/env python3

"""
.. module:: factoryOfWalkers
   :synopsis: facility that creates armies of randomWalkers

.. moduleauthor:: Wolfgang Waltenberger <wolfgang.waltenberger@gmail.com>  

"""

__all__ = [ "createWalkers" ]

import os, sys
from os import PathLike
from typing import Union, Dict, List
try:
    from torch import multiprocessing
except:
    import multiprocessing

def _run ( walker, catch_exceptions, seed ):
    if seed is not None:
        from ptools import helpers
        helpers.seedRandomNumbers( seed + walker.walkerid )
        print ( f"[factoryOfWalkers] setting random seed to {seed}" )
    if not catch_exceptions:
        walker.walk()
        return
    try:
        walker.walk(catch_exceptions)
    except Exception as e:
        import time
        with open("exceptions.log","a") as f:
            f.write ( f"time {time.asctime()}\n" )
            f.write ( f"walker {walker.walkerid} threw: {e}\n" )
            import traceback
            f.write ( f"traceback: {str(traceback.format_exc())}\n" )
            if hasattr ( walker.manipulator.M, "currentSLHA" ):
                f.write ( f"slha file was {walker.manipulator.M.currentSLHA}\n" )
        from colorama import Fore as ansi
        print ( f"{ansi.RED}walker {walker.walkerid} threw: {e}{ansi.RESET}\n" )

def startWalkers ( walkers : List, catch_exceptions : bool = False,
                   seed : Union[None,int] = None ) -> int:
    """ start the walkers

    :param catch_exceptions: If True will catch the exceptions and exit.
    :param seed: random seed number (optional)
    :returns: number of started walkers
    """
    processes=[]
    print("[factoryOfWalkers] startWalkers")
    if len(walkers) == 1: 
        _run ( walkers[0], catch_exceptions, seed )
        return 1
    for walker in walkers:
        p = multiprocessing.Process ( target=_run, args=( walker, catch_exceptions, seed ) )
        p.start()
        processes.append(p)
    for p in processes:
        p.join()
    return len(processes)

def writeMetaInfo ( rundir : str, meta : Dict ):
    """ write meta info of factory run to run.dict. complain if data is 
    different from previous info.
    """
    dictfile = f"{rundir}/run.dict"
    if os.path.exists ( dictfile ):
        oldmeta = {}
        with open ( dictfile, "rt" ) as f:
            txt = f.read()
            f.close()
            if len(txt)>0:
                oldmeta = eval(txt)
            for k,v in oldmeta.items():
                if not k in meta:
                    print ( f"[factoryOfWalkers] run's meta info changed: {k} was {v} now not in meta" )  
                    continue
                if meta[k] != v:
                    print ( f"[factoryOfWalkers] run's meta info changed: {k} was {v} not {meta[k]}" )
    else:
        with open ( dictfile, "wt" ) as f:
            from ptools.helpers import py_dumps
            ds = py_dumps ( meta, indent = 4 )
            #ds = ds.replace ( "false", "False" )
            #ds = ds.replace ( "true", "True" )
            f.write ( ds + "\n" )
            # f.write ( f"{meta!s}\n" )
            f.close()

def createWalkers ( rvars: dict ):
    """ a worker node to set up to run walkers

    rvars ( dict ):
      - nmin: the walker id of the first walker
      - nmax: the walker id + 1 of the last walker
      - continueFrom: start with protomodels given in the pickle file or hiscore dictionary file
      - cheatcode: in case this is not 0 or "no_cheat", we wish to start from a cheat model
      - rundir: overrride default rundir, if None use default
      - maxsteps: maximum number of steps to be taken
      - seed: random seed number (optional)
      - test_param_space: if true, run with constant K and TL (=1.0)
      - run_mcmc: if true, run mcmc walk without changing dimensions
      - catch_exceptions: If True will catch the exceptions and exit.
      - select: select only subset of results (all for all, em for efficiency 
        maps only, ul for upper limits only, alternatively select for txnames via
        e.g. "txnames:T1,T2", short names are recognized, e.g.
        "txnames:electroweakinos_offshell,T1"
      - cap_ssm: set the maximum value for all signal strength multipliers (default=100)
      - do_srcombine: if true, then also perform combinations, either via
                       simplified likelihoods or via pyhf
      - record_history: if True, then use history recorders
      - update_hiscores: if True, then finish your run and
                            after that run hiscore updater
      - stopTeleportationAfter: integer, stop teleportation after this step has 
        been reached. -1 or None means, dont run teleportation at all.
      - forbiddenparticles: an optional list of particles we wont touch in this
        run
      - templateSLHA: the template file that is used
      - allowN1N1Prod: allow N1 N1 production mode
      - susy_mode: susy mode, dont touch ssms
      - use_initialiser: if string, then interpret it as path to database
    """
    globals().update ( rvars ) # doesnt work for all
    dbpath = rvars["dbpath"]
    use_initialiser = rvars["use_initialiser"]
    #jmax = rvars["jmax"]
    continueFrom = rvars["continueFrom"]
    cheatcode = rvars["cheatcode"]
    do_srcombine = rvars["do_srcombine"]
    meta = { "dbpath": dbpath, "select": select, "do_srcombine": do_srcombine,
             "forbidden": forbiddenparticles, "templateSLHA": templateSLHA,
             "allowN1N1Prod": allowN1N1Prod, "susy_mode": susy_mode,
             "use_initialiser": use_initialiser }
    from builder.manipulator import Manipulator
    from ptools.moreHelpers import namesForSetsOfPids
    Manipulator.forbiddenparticles = namesForSetsOfPids ( forbiddenparticles )
    writeMetaInfo ( rundir, meta )

    if rundir != None and "<rundir>" in dbpath:
        dbpath=dbpath.replace("<rundir>", f"{rundir}/" )
    pfile, states = None, None
    if continueFrom == "default":
        continueFrom = f"{rundir}/states.dict" 
        if not os.path.exists ( continueFrom ):
            continueFrom = "default"
    if continueFrom.lower() not in [ "none", "" ]:
        if not os.path.exists ( continueFrom ):
            print ( f"[factoryOfWalkers] error: supplied a save states file ,,{continueFrom}'', but it doesnt exist" )
        else:
            import pickle
            try:
                if continueFrom.endswith ( ".dict" ):
                    with open( continueFrom, "rt" ) as f:
                        states = eval ( f.read() )
                else:
                    with open ( continueFrom, "rb" ) as f:
                        states = pickle.load ( f )
                pfile = continueFrom
            except Exception as e:
                print ( f"error when trying to load file {continueFrom}: {e}" )
                pfile = None
    walkers = []
    #Set random seed
    from walker.randomWalker import RandomWalker
    for i in range(nmin,nmax):
        if pfile is None:
            import time
            import socket
            hostname = socket.gethostname().replace(".cbe.vbc.ac.at","")
            atime = time.strftime('%H:%M:%S')
            label = f"[factoryOfWalkers:{hostname};{atime}]"
            print ( f"{label} starting {i} @ {rundir} with cheatcode {cheatcode}" )
            rvars["walkerid"]= i
            w = RandomWalker( rvars )
            walkers.append ( w )
        elif pfile.endswith(".hi") or pfile.endswith(".pcl"):
            nstates = len(states )
            ctr = i % nstates
            print ( f"[factoryOfWalkers] fromModel {i}: loading {ctr}/{nstates}" )
            rvars["walkerid"]=i
            w = RandomWalker.fromProtoModel ( states[ctr], rvars )
            walkers.append ( w )
        else:
            nstates = len(states )
            ctr = i % nstates
            print ( f"[factoryOfWalkers] fromDict {i}: loading {ctr}/{nstates}" )
            rvars["walkerid"]=i
            w = RandomWalker.fromDictionary ( states[ctr], rvars )
            walkers.append ( w )
    #start running walkers
    startWalkers ( walkers, catch_exceptions=catch_exceptions, seed=seed )
    if update_hiscores:
        import time
        from ptools import updateHiscores
        ctAttempts = 0 ## count how often we tried
        succeeded = False
        while ctAttempts < 7:
            steps = updateHiscores.countSteps( writeSubmitFile = False )
            if not type(steps)==tuple:
                print ( f"[factoryOfWalkers] been asked to update hiscores, but dont understand steps {steps}" )
                sys.exit(-1)
            print ( f"[factoryOfWalkers] been asked to update hiscores: {steps[0]} == {nmax*maxsteps}" )
            ctAttempts += 1
            if steps[0] == nmax*maxsteps: ## are we last?
                updateHiscores.loop ( rundir = rundir, maxruns=1,
                                      doPlots=False, uploadTo="latest" )
                succeeded = True
                break
            else:
                time.sleep ( (ctAttempts**2+1)*180 )
        if succeeded:
            print ( "[factoryOfWalkers] ran updater successfully." )
        else:
            print ( f"[factoryOfWalkers] tried more {ctAttempts} times. stop trying." )

if __name__ == "__main__":
    import sys
    sys.path.insert(0,"../")
    sys.path.insert(0,f"{os.environ['HOME']}/git/smodels/")
    sys.path.insert(0,"../../")
    from walker.randomWalker import RandomWalker
    s = "txnames:TChiWZ,TChiWZoff,TChiWW,TChiWWoff,TChiWH,TChiH,TChiZZ,TSlepSlep"
    s = "all"
    dbpath = "./default.pcl"
    dbpath = "official"
    w = RandomWalker( walkerid=0, nsteps = 200, 
                      dbpath=dbpath, cheatcode="no_cheat", select=s,
                      rundir="./", seed = None )
    w.walk()

