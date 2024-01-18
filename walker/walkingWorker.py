#!/usr/bin/env python3

import os, sys
from os import PathLike
try:
    from torch import multiprocessing
except:
    import multiprocessing

def _run ( walker, catchem, seed ):
    if seed is not None:
        from ptools import helpers
        helpers.seedRandomNumbers( seed + walker.walkerid )
        print ( f"[walkingWorker] setting random seed to {seed}" )
    if not catchem:
        walker.walk()
        return
    try:
        walker.walk(catchem)
    except Exception as e:
        import time
        with open("exceptions.log","a") as f:
            f.write ( "time %s\n" % time.asctime() )
            f.write ( "walker %d threw: %s\n" % ( walker.walkerid, e ) )
            if hasattr ( walker.model, "currentSLHA" ):
                f.write ("slha file was %s\n" % walker.model.currentSLHA )
        import colorama
        print ( "%swalker %d threw: %s%s\n" % ( colorama.Fore.RED, walker.walkerid, e, colorama.Fore.RESET ) )

def startWalkers ( walkers, catchem=False, seed = None ):

    processes=[]
    for walker in walkers:
        p = multiprocessing.Process ( target=_run, args=( walker, catchem, seed ) )
        p.start()
        processes.append(p)
    for p in processes:
        p.join()


def main( nmin, nmax, continueFrom : PathLike,
          dbpath : PathLike = "<rundir>/database.pcl",
          cheatcode = 0, dump_training = False, rundir=None, maxsteps = 10000,
          nevents = 100000, seed = None,  catchem=True, select="all",
          do_combine = False, record_history = False, update_hiscores = False,
          stopTeleportationAfter : int = -1 ):
    """ a worker node to set up to run walkers
    :param nmin: the walker id of the first walker
    :param nmax: the walker id + 1 of the last walker
    :param continueFrom: start with protomodels given in the pickle file or hiscore dictionary file
    :param cheatcode: in case we wish to start from a cheat model
    :param dump_training: dump training data for the NN
    :param rundir: overrride default rundir, if None use default
    :param maxsteps: maximum number of steps to be taken
    :param nevents: number of MC events when computing cross-sections
    :param seed: random seed number (optional)
    :param catchem: If True will catch the exceptions and exit.
    :param select: select only subset of results (all for all, em for efficiency maps only,
            ul for upper limits only, alternatively select for txnames via
            e.g. "txnames:T1,T2"
    :param do_combine: if true, then also perform combinations, either via
                       simplified likelihoods or via pyhf
    :param record_history: if True, then use history recorders
    :param update_hiscores: if True, then finish your run and
                            after that run hiscore updater
    :param stopTeleportationAfter: integer, stop teleportation after this step has been
            reached. -1 or None means, dont run teleportation at all.
    """

    if rundir != None and "<rundir>" in dbpath:
        dbpath=dbpath.replace("<rundir>","%s/" % rundir )
    pfile, states = None, None
    if continueFrom == "default":
        continueFrom = f"{rundir}/states.dict" 
        if not os.path.exists ( continueFrom ):
            continueFrom = "default"
    if continueFrom.lower() not in [ "none", "" ]:
        if not os.path.exists ( continueFrom ):
            print ( f"[walkingWorker] error: supplied a save states file ,,{continueFrom}'', but it doesnt exist" )
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
            print ( "[walkingWorker] starting %d @ %s with cheatcode %d" % ( i, rundir, cheatcode ) )
            w = RandomWalker( walkerid=i, nsteps = maxsteps,
                              dump_training = dump_training,
                              dbpath=dbpath, cheatcode=cheatcode, select=select,
                              rundir=rundir, nevents=nevents, do_combine = do_combine,
                              record_history=record_history, seed=seed,
                              stopTeleportationAfter = stopTeleportationAfter )
            walkers.append ( w )
        elif pfile.endswith(".hi") or pfile.endswith(".pcl"):
            nstates = len(states )
            ctr = i % nstates
            print ( "[walkingWorker] fromModel %d: loading %d/%d" % ( i, ctr, nstates ) )
            w = RandomWalker.fromProtoModel ( states[ctr], strategy = "aggressive",
                    walkerid = i, nsteps = maxsteps, dump_training=dump_training,
                    expected = False, select = select, dbpath = dbpath,
                    rundir = rundir, do_combine = do_combine, seed = seed,
                    stopTeleportationAfter = stopTeleportationAfter )
            walkers.append ( w )
        else:
            nstates = len(states )
            ctr = i % nstates
            print ( "[walkingWorker] fromDict %d: loading %d/%d" % ( i, ctr, nstates ) )
            w = RandomWalker.fromDictionary ( states[ctr], nsteps = maxsteps,
                    strategy = "aggressive", walkerid = i,
                    dump_training=dump_training, dbpath = dbpath, expected = False,
                    select = select, rundir = rundir, nevents = nevents,
                    do_combine = do_combine, 
                    seed = seed, stopTeleportationAfter = stopTeleportationAfter )
            walkers.append ( w )
    startWalkers ( walkers, catchem=catchem, seed=seed )
    if update_hiscores:
        import time
        from ptools import updateHiscores
        ctAttempts = 0 ## count how often we tried
        succeeded = False
        while ctAttempts < 7:
            steps = updateHiscores.countSteps( writeSubmitFile = False )
            if not type(steps)==tuple:
                print ( "[walkingWorker] been asked to update hiscores, but dont understand steps %s" % steps )
                sys.exit(-1)
            print ( "[walkingWorker] been asked to update hiscores: %d == %d" % \
                    ( steps[0], nmax*maxsteps ) )
            ctAttempts += 1
            if steps[0] == nmax*maxsteps: ## are we last?
                updateHiscores.main ( rundir = rundir, maxruns=1,
                                      doPlots=False, uploadTo="latest" )
                succeeded = True
                break
            else:
                time.sleep ( (ctAttempts**2+1)*180 )
        if succeeded:
            print ( "[walkingWorker] ran updater successfully." )
        else:
            print ( f"[walkingWorker] tried more {ctAttempts} times. stop trying." )

if __name__ == "__main__":
    import sys
    sys.path.insert(0,"../")
    sys.path.insert(0,"../../")
    from walker.randomWalker import RandomWalker
    s = "txnames:TChiWZ,TChiWZoff,TChiWW,TChiWWoff,TChiWH,TChiH,TChiZZ,TSlepSlep"
    s = "all"
    dbpath = "./default.pcl"
    dbpath = "official"
    # dbpath = "~/git/smodels-database"
    w = RandomWalker( walkerid=0, nsteps = 200, dump_training = False,
                      dbpath=dbpath, cheatcode=0, select=s,
                      rundir="./", nevents=1000, seed = None )
    w.walk()

