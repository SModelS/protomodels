#!/usr/bin/env python3

import os
try:
    from torch import multiprocessing
except:
    import multiprocessing

def _run ( walker, catchem, seed=None ):
    from ptools import helpers
    #Set random seed
    if seed is not None:
        helpers.seedRandomNumbers(seed)
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

def startWalkers ( walkers, seed=None,  catchem=False):

    processes=[]
    for walker in walkers:
        p = multiprocessing.Process ( target=_run, args=( walker, catchem, seed ) )
        p.start()
        processes.append(p)
    for p in processes:
        p.join()


def main( nmin, nmax, cont,
          dbpath = "<rundir>/database.pcl",
          cheatcode = 0, dump_training = False, rundir=None, maxsteps = 10000,
          nevents = 100000, seed = None,  catchem=True, select="all",
          do_combine = False, record_history = False ):
    """ a worker node to set up to run walkers
    :param nmin: the walker id of the first walker
    :param nmax: the walker id + 1 of the last walker
    :param cont: start with protomodels given in the pickle file 'cont'
    :param cheatcode: in case we wish to start from a cheat model
    :param dump_training: dump training data for the NN
    :param rundir: overrride default rundir, if None use default
    :param maxsteps: maximum number of steps to be taken
    :param nevents: number of MC events when computing cross-sections
    :param seed: random seed number (optional)
    :param catchem: If True will catch the exceptions and exit.
    :param select: select only subset of results (all for all, em for efficiency maps only, ul for upper limits only, alternatively select for txnames via e.g. "txnames:T1,T2"
    :param do_combine: if true, then also perform combinations, either via
                       simplified likelihoods or via pyhf
    :param record_history: if True, then use history recorders
    """

    if rundir != None and "<rundir>" in dbpath:
        dbpath=dbpath.replace("<rundir>","%s/" % rundir )
    pfile, states = None, None
    if cont == "default":
        cont = "%s/states.dict" % rundir
        if not os.path.exists ( cont ):
            cont = "default"
    if cont.lower() not in [ "none", "" ]:
        if not os.path.exists ( cont ):
            print ( "[walkingWorker] error: supplied a save states file ,,%s'', but it doesnt exist" % cont )
        else:
            import pickle
            try:
                if cont.endswith ( ".dict" ):
                    with open( cont, "rt" ) as f:
                        states = eval ( f.read() )
                else:
                    with open ( cont, "rb" ) as f:
                        states = pickle.load ( f )
                pfile = cont
            except Exception as e:
                print ( "error when trying to load pickle file %s: %s" % ( cont, e ) )
                pfile = None
    # print ( "[walkingWorker] called main with cont='%s', pfile='%s'." % ( cont, pfile ) )

    # print ( "[walkingWorker] I am already inside the python script! Hostname is", socket.gethostname()  )
    walkers = []
    from walker.randomWalker import RandomWalker
    for i in range(nmin,nmax):
        if pfile is None:
            print ( "[walkingWorker] starting %d @ %s with cheatcode %d" % ( i, rundir, cheatcode ) )
            w = RandomWalker( walkerid=i, nsteps = maxsteps, 
                              dump_training = dump_training,
                              dbpath=dbpath, cheatcode=cheatcode, select=select, 
                              rundir=rundir, nevents=nevents, do_combine = do_combine,
                              record_history=record_history )
            walkers.append ( w )
        elif pfile.endswith(".hi") or pfdile.endswith(".pcl"):
            nstates = len(states )
            ctr = i % nstates
            print ( "[walkingWorker] fromModel %d: loading %d/%d" % ( i, ctr, nstates ) )
            w = RandomWalker.fromProtoModel ( states[ctr], strategy = "aggressive",
                    walkerid = i, nsteps = maxsteps, dump_training=dump_training, 
                    expected = False, select = select, dbpath = dbpath, 
                    rundir = rundir, do_combine = do_combine, record_history = record_history )
            walkers.append ( w )
        else:
            nstates = len(states )
            ctr = i % nstates
            print ( "[walkingWorker] fromDict %d: loading %d/%d" % ( i, ctr, nstates ) )
            w = RandomWalker.fromDictionary ( states[ctr], nsteps = maxsteps, 
                    strategy = "aggressive", walkerid = i, 
                    dump_training=dump_training, dbpath = dbpath, expected = False, 
                    select = select, rundir = rundir, nevents = nevents,
                    do_combine = do_combine, record_history = record_history )
            walkers.append ( w )
    startWalkers ( walkers, seed=seed, catchem=catchem )

if __name__ == "__main__":
    import sys
    sys.path.insert(0,"../")
    from walker.randomWalker import RandomWalker
    s = "txnames:TChiWZ,TChiWZoff,TChiWW,TChiWWoff,TChiWH,TChiH,TChiZZ,TSlepSlep"
    s = "all"
    w = RandomWalker( walkerid=0, nsteps = 10, dump_training = False,
                      dbpath="./default.pcl", cheatcode=0, select=s, 
                      rundir="./", nevents=1000 )
    w.walk()
    
