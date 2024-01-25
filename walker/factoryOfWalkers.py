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
            f.write ( "time %s\n" % time.asctime() )
            f.write ( "walker %d threw: %s\n" % ( walker.walkerid, e ) )
            if hasattr ( walker.model, "currentSLHA" ):
                f.write ("slha file was %s\n" % walker.model.currentSLHA )
        import colorama
        print ( "%swalker %d threw: %s%s\n" % ( colorama.Fore.RED, walker.walkerid, e, colorama.Fore.RESET ) )

def startWalkers ( walkers : List, catch_exceptions : bool = False,
                   seed : Union[None,int] = None ) -> int:
    """ start the walkers

    :param catch_exceptions: If True will catch the exceptions and exit.
    :param seed: random seed number (optional)
    :returns: number of started walkers
    """
    processes=[]
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
            f.write ( str(meta)+"\n" )
            f.close()

def createWalkers( nmin : int , nmax : int, continueFrom : PathLike,
          dbpath : PathLike = "official", cheatcode : int = 0, 
          rundir : Union[None,str] = None, maxsteps : int = 10000,
          nevents : int = 100000, seed : Union[None,int] = None, 
          catch_exceptions : bool = True, select : str = "all",
          do_combine : bool = False, record_history : bool = False, 
          update_hiscores : bool = False, stopTeleportationAfter : int = -1,
          forbiddenparticles : List[int] = [] ):
    """ a worker node to set up to run walkers

    :param nmin: the walker id of the first walker
    :param nmax: the walker id + 1 of the last walker
    :param continueFrom: start with protomodels given in the pickle file or hiscore dictionary file
    :param cheatcode: in case we wish to start from a cheat model
    :param rundir: overrride default rundir, if None use default
    :param maxsteps: maximum number of steps to be taken
    :param nevents: number of MC events when computing cross-sections
    :param seed: random seed number (optional)
    :param catch_exceptions: If True will catch the exceptions and exit.
    :param select: select only subset of results (all for all, em for efficiency 
    maps only, ul for upper limits only, alternatively select for txnames via
    e.g. "txnames:T1,T2", short names are recognized, e.g.
    "txnames:electroweakinos_offshell,T1"

    :param do_combine: if true, then also perform combinations, either via
                       simplified likelihoods or via pyhf
    :param record_history: if True, then use history recorders
    :param update_hiscores: if True, then finish your run and
                            after that run hiscore updater
    :param stopTeleportationAfter: integer, stop teleportation after this step has 
    been reached. -1 or None means, dont run teleportation at all.
    :param forbiddenparticles: an optional list of particles we wont touch in this
    run
    """
    meta = { "dbpath": dbpath, "select": select, "do_combine": do_combine }
    from builder.manipulator import Manipulator
    Manipulator.forbiddenparticles = forbiddenparticles
    writeMetaInfo ( rundir, meta )

    if rundir != None and "<rundir>" in dbpath:
        dbpath=dbpath.replace("<rundir>","%s/" % rundir )
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
            print ( f"[factoryOfWalkers:{hostname};{time.strftime('%H:%M:%S')}] starting {i} @ {rundir} with cheatcode {cheatcode}" )
            w = RandomWalker( walkerid=i, nsteps = maxsteps,
                              dbpath=dbpath, cheatcode=cheatcode, select=select,
                              rundir=rundir, nevents=nevents, do_combine = do_combine,
                              record_history=record_history, seed=seed,
                              stopTeleportationAfter = stopTeleportationAfter )
            walkers.append ( w )
        elif pfile.endswith(".hi") or pfile.endswith(".pcl"):
            nstates = len(states )
            ctr = i % nstates
            print ( "[factoryOfWalkers] fromModel %d: loading %d/%d" % ( i, ctr, nstates ) )
            w = RandomWalker.fromProtoModel ( states[ctr], strategy = "aggressive",
                    walkerid = i, nsteps = maxsteps,
                    expected = False, select = select, dbpath = dbpath,
                    rundir = rundir, do_combine = do_combine, seed = seed,
                    stopTeleportationAfter = stopTeleportationAfter )
            walkers.append ( w )
        else:
            nstates = len(states )
            ctr = i % nstates
            print ( "[factoryOfWalkers] fromDict %d: loading %d/%d" % ( i, ctr, nstates ) )
            w = RandomWalker.fromDictionary ( states[ctr], nsteps = maxsteps,
                    strategy = "aggressive", walkerid = i, dbpath = dbpath, 
                    expected = False, select = select, rundir = rundir, 
                    nevents = nevents, do_combine = do_combine, 
                    seed = seed, stopTeleportationAfter = stopTeleportationAfter )
            walkers.append ( w )
    startWalkers ( walkers, catch_exceptions=catch_exceptions, seed=seed )
    if update_hiscores:
        import time
        from ptools import updateHiscores
        ctAttempts = 0 ## count how often we tried
        succeeded = False
        while ctAttempts < 7:
            steps = updateHiscores.countSteps( writeSubmitFile = False )
            if not type(steps)==tuple:
                print ( "[factoryOfWalkers] been asked to update hiscores, but dont understand steps %s" % steps )
                sys.exit(-1)
            print ( "[factoryOfWalkers] been asked to update hiscores: %d == %d" % \
                    ( steps[0], nmax*maxsteps ) )
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
    from walker.randomWalker import RandomWalker
    s = "txnames:TChiWZ,TChiWZoff,TChiWW,TChiWWoff,TChiWH,TChiH,TChiZZ,TSlepSlep"
    s = "all"
    dbpath = "./default.pcl"
    dbpath = "official"
    # dbpath = "~/git/smodels-database"
    w = RandomWalker( walkerid=0, nsteps = 200, 
                      dbpath=dbpath, cheatcode=0, select=s,
                      rundir="./", nevents=1000, seed = None )
    w.walk()

