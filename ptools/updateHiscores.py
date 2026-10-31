#!/usr/bin/env python3

""" simple script that perpetually updates the hiscores_global.cache file,
and the plots """

__all__ = [ "loop", "didSRCombine" ]

import time, types, sys, os, subprocess
from os import PathLike
from typing import Union, Dict

from pbase.runEnviron import RunEnviron

def setup( rundir = None ):
    if "CODEDIR" in os.environ:
        codedir = os.environ['CODEDIR']
        sys.path.insert(0,f"{codedir}smodels/" )
        sys.path.insert(0,f"{codedir}smodels-utils/" )
        sys.path.insert(0,f"{codedir}/protomodels/" )
    if rundir != None:
        rundir = rundir.replace ( "~", os.environ["HOME"] )
        os.chdir ( rundir )
        return rundir
    rundir = os.environ['RUNDIR']
    # rundir = "/mnt/hephy/pheno/ww/rundir/"
    # rundir = "./"
    if os.path.exists ( "./rundir.conf" ):
        with open ( "./rundir.conf" ) as f:
            rundir = f.read().strip()
    rundir = rundir.replace ( "~", os.environ["HOME"] )
    os.chdir ( rundir )
    return rundir

def didSRCombine ( rundir : str ) -> Union[None,bool]:
    """ find out whether or not they did sr-combine in this run """
    ret = None
    fname = f"{rundir}/run.dict"
    if os.path.exists ( fname ):
        with open ( fname, "rt" ) as f:
            txt = f.read()
            d = eval(txt)
            if "do_srcombine" in d:
                ret = d["do_srcombine"]
    return ret

def countSteps( printout = True, writeSubmitFile = False, doSubmit = False ):
    """ count the number of steps taken accoring to walker logs
    :param printout: print out statistics
    :param writeSubmitFile: write a submit file for the non-finished jobs
    :param doSubmit: if True, then do submit jobs without further asking
    """
    import glob
    files = glob.glob("walker*log")
    steps = {}
    g = None
    if writeSubmitFile:
        g = open ( "submit.sh", "wt" )
        g.write ( "#!/bin/sh\n\n" )
    files.sort()
    for f in files:
        nr = int ( f.replace("walker","").replace(".log","") )
        h = open ( f, "rt" )
        lines = h.readlines()
        h.close()
        slurmid = 0
        if "slurm jobid" in lines[0]:
            p = lines[0].rfind ( "jobid " )
            slurmid = int ( lines[0][p+6:] )
        for cl, line in enumerate ( lines[::-1] ):
            if "Step" in line:
                laststep = line[line.find("Step")+5:]
                for c in [ "/", ":", " has", " " ]:
                    if c in laststep:
                        laststep = laststep[:laststep.find(c)]
                laststep = int ( laststep.strip() )
                slurmfile = ""
                if slurmid > 0:
                    slurmfile = f"{os.environ['HOME']}/outputs/walk-{slurmid}.out"
                #print ( nr, laststep )
                steps[nr]=laststep
                if writeSubmitFile and laststep < 1000:
                    rundir = os.getcwd()
                    if rundir.endswith("/"):
                        rundir=rundir[:-1]
                    p = rundir.rfind("/")
                    rundir = rundir[p+1:]
                    g.write ( f"rm -rf {os.environ['HOME']}/{rundir}/H{nr}.hi\n" )
                    g.write ( f"./slurm_walk.py -R {rundir} -n {int(nr)} -N {int(nr + 1)} -M 1000\n" )
                    if slurmfile != "":
                        g.write ( f"rm -rf {slurmfile}\n" )
                break
    keys = list ( steps.keys() )
    keys.sort()
    tots = 0
    finished = []
    for k in keys:
        tots += steps[k]
        if printout and steps[k] < 1000:
            print ( f"walker {int(k)}: {int(steps[k])}" )
        if steps[k] == 1000:
            finished.append ( k )
    if printout and len(finished)>0:
        print ( f"{len(finished)} walkers finished: {','.join(list(map(str, finished)))}" )
    if printout:
        print ( f"we have {len(keys)} entries, total of {tots} steps." )
    if writeSubmitFile:
        keys = list ( steps.keys() )
        keys.sort()
        for k in range(0,50):
            if not k in keys:
                rundir = os.getcwd()
                if rundir.endswith("/"):
                    rundir=rundir[:-1]
                p = rundir.rfind("/")
                rundir = rundir[p+1:]
                g.write ( f"./slurm_walk.py -R {rundir} -n {int(k)} -N {int(k + 1)} -M 1000\n" )
        g.close()
        os.chmod ( "submit.sh", 0o755 )
        cmd = f"cp submit.sh {os.environ['HOME']}"
        subprocess.getoutput ( cmd )
        if doSubmit:
            cmd = f"cd {os.environ['HOME']}; ./submit.sh; cd -"
            a = subprocess.getoutput ( cmd )
            print ( a )
    return tots,steps

def updateHiscores( dictfile : os.PathLike = "{rundir}/hiscores_global.dict",
                cachefile : os.PathLike = "{rundir}/hiscores_global.cache",
                environ : Union[RunEnviron,None ] = None,
                walkerid : Union[str,int] = 0,
                hiscore_nr : int = 0 ) -> Dict:
    """ update the hiscores FIXME +
    :param hiscore_nr: which hiscore, 0 is the firsst
    """
    assert environ != None, "set RunEnviron"
    assert type(environ) != str, "set RunEnviron, not str"
    dictfile = dictfile.replace("{rundir}",environ.rundir)
    cachefile = cachefile.replace("{rundir}",environ.rundir)
    from ptools import hiscoreTools
    hi = hiscoreTools.fetchHiscoresObj ( dictfile, cachefile, environ = environ,
            walkerid = walkerid )

    from builder.manipulator import Manipulator
    D = Manipulator ( hi.hiscores[hiscore_nr], environ=environ ).writeDictFile ( None )
    D["model"]=hi.hiscores[hiscore_nr]
    return D

def plot( TL : float, K : float, environ : RunEnviron, upload : str ="230",
          verbose : bool = False,
          dictfile : str = "{rundir}/hiscores_global.dict",
          walkerid : Union[str,int] = 0,
          git_commit : bool = False ):
    """ create all hiscore plots

    :param upload: the "label" of the upload. determines the directory name at
    https://smodels.github.io/protomodels/
    Typically: latest, official, 230, ...
    :param dbpath: path to database, look for default.pcl in rundir by default
    :param verbose: be verbose, if true
    """
    dictfile = dictfile.replace("{rundir}",environ.rundir)
    from plotting import plotHiscore
    from argparse import Namespace
    args = Namespace()
    args.upload = upload
    args.number = 0
    args.verbose = verbose
    args.detailed = False
    args.destinations = False
    args.hiscorefile = dictfile
    args.environ = environ
    args.verbosity = "info"
    args.horizontal = False
    args.html = True
    args.ruler = False # True
    args.decays = False # True
    args.masses_decays = True # True
    args.predictions = True
    args.tex = False
    args.keep = False
    args.commit = git_commit
    args.walkerid = walkerid
    plotHiscore.runPlotting ( args )

# default options, shared by updateOnce and loop. Any of these can be
# overridden via keyword arguments, see updateOnce
hiscoresOptions = {
    "createPlots" : True,
    "uploadTo"    : "temp",
    "verbose"     : False,
    "environ"     : None,
    "dictfile"    : "{rundir}/hiscores_global.dict",
    "cachefile"   : "{rundir}/hiscores_global.cache",
    "walkerid"    : 0,
    "git_commit"  : True,
    "hiscore_nr"  : 0,
}

def updateOnce( rundir : str, Kold : float, **kwargs ) -> float:
    """ perform a single update of hiscores_global.cache, and create
    the plots

    :param rundir: the run directory
    :param Kold: the value of K from the previous iteration
    :param kwargs: optional settings, any key of `hiscoresOptions` can be
    overridden. In particular: createPlots (suppress plotting if False),
    uploadTo (upload plots to "~/git/smodels.github.io/<uploadTo>"),
    verbose (verbosity), environ (the run environ), walkerid (log
    everything as walker #walkerid), hiscore_nr (which hiscore, 0 is the
    first), dictfile, cachefile, git_commit
    :returns: the updated value of Kold

    As a side effect, a dictionary holding the timestamp, the
    dictfile, and a human readable time is written to "uphi.log".
    """
    opts = dict ( hiscoresOptions )
    opts.update ( kwargs )
    createPlots = opts["createPlots"]
    TLfile = f"{rundir}/TLold.conf"
    Kfile = f"{rundir}/Kold.conf"
    D = updateHiscores( opts["dictfile"], opts["cachefile"], opts["environ"],
            walkerid=opts["walkerid"], hiscore_nr = opts["hiscore_nr"] )
    TL, step, K = float("nan"),0,float("nan")
    model = D["model"]
    if "TL" in D:
        TL = D["TL"]
    if "step" in D:
        step = D["step"]
    if "K" in D:
        K = D["K"]
    if K is not None and K > Kold + 1e-10: #  + .001:
        from builder.manipulator import Manipulator
        m = Manipulator ( model, opts["environ"] )
        T=str(int(time.time()))
        m.writeDictFile ( f"pmodel-{T}.dict", comment="history keeper" )
        with open ( f"{rundir}history.txt", "at" ) as f:
            f.write ( f"{time.asctime()}, step={step}, TL={TL:.4f}, K={K:.4f}, t={T}\n" )
            f.close()
        with open ( TLfile, "wt" ) as f:
            f.write ( f"{str(TL)}\n" )
            f.close()
        with open ( Kfile, "wt" ) as f:
            f.write ( f"{str(K)}\n" )
            f.close()
        Kold = K
    if createPlots:
        plot ( TL, K, opts["environ"], opts["uploadTo"], opts["verbose"],
               dictfile = opts["dictfile"], walkerid = opts["walkerid"],
               git_commit = opts["git_commit"] )
    else:
        print ( "[updateHiscores] was not asked to create plots" )
    if os.path.exists ( Kfile ): ## so we can meddle from outside
        with open ( Kfile, "rt" ) as f:
            try:
                Kold = float ( f.read().strip() )
            except ValueError as e:
                Kold = -90.
    now = time.time()
    info = { "timestamp" : int(now), "dictfile" : opts["dictfile"],
             "time" : time.asctime ( time.localtime ( now ) ) }
    from ptools.helpers import py_dumps
    with open ( "uphi.log", "wt" ) as f:
        f.write ( py_dumps ( info ) + "\n" )
    return Kold

def loop( maxruns : Union[None,int] = 3, **kwargs ):
    """ loop (maxruns times) that updates hiscores_global.cache

    :param maxruns: maximally iterate that many times, if None then loop endlessly
    :param kwargs: optional settings, handed on to `updateOnce`. Any key of
    `hiscoresOptions` can be overridden, in particular environ (the run
    environ, set to None for default (run.dict)) and createPlots (if
    False, suppress plotting). A rundir keyword overrides the rundir of
    the environ.
    """
    opts = dict ( hiscoresOptions )
    opts.update ( kwargs )
    if "doPlots" in opts: ## backwards compatibility, doPlots -> createPlots
        opts["createPlots"] = opts.pop ( "doPlots" )
    environ = opts["environ"]
    if environ == None:
        environ = RunEnviron()
    if type(environ)==str:
        environ = RunEnviron ( environ )
    opts["environ"] = environ
    rundir = opts.pop ( "rundir", None )
    if rundir == None:
        rundir = environ.rundir
    rundir = setup ( rundir )
    i = 0
    Kold = -90.
    Kfile = f"{rundir}/Kold.conf"
    if os.path.exists ( Kfile ):
        with open ( Kfile, "rt" ) as f:
            try:
                Kold = float ( f.read().strip() )
            except ValueError as e:
                pass
    while True:
        i+=1
        if maxruns != None and i > maxruns:
            break
        if i>1:
            time.sleep(60.)
        Kold = updateOnce( rundir, Kold, **opts )

if __name__ == "__main__":
    import argparse
    argparser = argparse.ArgumentParser('run the updateHiscores loop')
    argparser.add_argument ( '--maxruns',
            help='maxruns', type=int, default=3 )
    argparser.add_argument ( '--createPlots',
            help='create the plots', action="store_true" )
    argparser.add_argument ( '--uploadTo',
            help='where to upload', type=str, default=None )
    argparser.add_argument ( '--environ',
            help='path to environment file', type=str,
            default=None )
    argparser.add_argument ( '--dictfile',
            help='dictfile', type=str, default="{rundir}/hiscores_global.dict" )
    args = argparser.parse_args()
    vargs=vars(args)
    loop( **vargs )
