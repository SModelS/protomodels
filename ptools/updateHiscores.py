#!/usr/bin/env python3

""" simple script that perpetually updates the hiscores_global.cache file,
and the plots """

__all__ = [ "loop", "didSRCombine" ]

import time, types, sys, os, subprocess
from os import PathLike
from typing import Union, Dict

from base.runEnviron import RunEnviron

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
                    g.write ( f"./slurm.py -R {rundir} -n {int(nr)} -N {int(nr + 1)} -M 1000\n" )
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
                g.write ( f"./slurm.py -R {rundir} -n {int(k)} -N {int(k + 1)} -M 1000\n" )
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

def loop( maxruns : Union[None,int] = 3, createPlots : bool=True,
          uploadTo : str = "temp", environ : Union[RunEnviron, None] = None,
          verbose : bool = False, dictfile = "{rundir}/hiscores_global.dict",
          cachefile = "{rundir}/hiscores_global.cache",
          walkerid : Union[str,int] = 0, git_commit : bool = True,
          hiscore_nr : int = 0 ):
    """ loop (maxruns times) that updates hiscores_global.cache

    :param maxruns: maximally iterate that many times, if None then loop endlessly
    :param createPlots: if False, suppress plotting
    :param uploadTo: upload plots to directory "~/git/smodels.github.io/<uploadTo>"
    :param verbose: verbosity
    :param environ: the run environ, set to None for default (run.dict)
    :param walkerid: log everything as walker #walkerid
    :param hiscore_nr: which hiscore, 0 is the firsst
    """
    if environ == None:
        environ = RunEnviron()
    rundir = setup( environ.rundir )
    i = 0
    TL, TLold, step, K, Kold = 0., 0., 0, -90., -90.
    TLfile = f"{rundir}/TLold.conf"
    if os.path.exists ( TLfile ):
        with open ( TLfile, "rt" ) as f:
            TLold = float ( f.read().strip() )
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
        D = updateHiscores( dictfile, cachefile, environ,
                walkerid=walkerid, hiscore_nr = hiscore_nr )
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
            m = Manipulator ( model, environ )
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
            TLold = TL
            Kold = K
        if createPlots:
            plot ( TL, K, environ, uploadTo, verbose,
                   dictfile = dictfile, walkerid = walkerid,
                   git_commit = git_commit )
        else:
            print ( "[updateHiscores] was not asked to create plots" )
        if os.path.exists ( Kfile ): ## so we can meddle from outside
            with open ( Kfile, "rt" ) as f:
                try:
                    Kold = float ( f.read().strip() )
                except ValueError as e:
                    Kold = -90.

if __name__ == "__main__":
    loop()
