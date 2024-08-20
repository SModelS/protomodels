#!/usr/bin/env python3

""" simple script that perpetually updates the hiscores.cache file,
and the plots """

__all__ = [ "loop", "didSRCombine" ]

import time, types, sys, os, subprocess
from os import PathLike
from typing import Union, Dict

def setup( rundir = None ):
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
                    g.write ( "./slurm.py -R %s -n %d -N %d -M 1000\n" % \
                              ( rundir, nr, nr+1 ) )
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
            print ( "walker %d: %d" % ( k, steps[k] ) )
        if steps[k] == 1000:
            finished.append ( k )
    if printout and len(finished)>0:
        print ( "%d walkers finished: %s" % \
                (len(finished), ",".join( list(map(str,finished ) ) ) ) )
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
                g.write ( "./slurm.py -R %s -n %d -N %d -M 1000\n" % \
                          ( rundir, k, k+1 ) )
        g.close()
        os.chmod ( "submit.sh", 0o755 )
        cmd = "cp submit.sh %s" % os.environ["HOME"]
        subprocess.getoutput ( cmd )
        if doSubmit:
            cmd = "cd %s; ./submit.sh; cd -" % os.environ["HOME"]
            a = subprocess.getoutput ( cmd )
            print ( a )
    return tots,steps

def updateHiscores( rundir : Union[None,PathLike] = None,
                    dbpath : Union[None,PathLike] = None,
                    do_srcombine : bool = False ) -> Dict:
    dictfile = f"{rundir}/hiscores.dict"
    cachefile = f"{rundir}/hiscores.cache"
    from ptools import hiscoreTools
    hi = hiscoreTools.fetchHiscoresObj ( dictfile, cachefile, dbpath )
    from builder.manipulator import Manipulator
    D = Manipulator ( hi.hiscores[0] ).writeDictFile ( None )
    D["model"]=hi.hiscores[0]
    return D

def plot( TL : float, K : float, rundir : os.PathLike, upload : str ="230",
          dbpath : str = "official", verbose : bool = False ):
    """ create all hiscore plots

    :param upload: the "label" of the upload. determines the directory name at
    https://smodels.github.io/protomodels/
    Typically: latest, official, 230, ...
    :param dbpath: path to database, look for default.pcl in rundir by default
    :param verbose: be verbose, if true
    """
    from plotting import plotHiscore
    from argparse import Namespace
    args = Namespace()
    args.upload = upload
    args.number = 0
    args.verbose = verbose
    args.detailed = False
    args.destinations = False
    args.hiscorefile = f"{rundir}/hiscores.dict"
    args.dbpath = dbpath.replace("@rundir@",rundir )
    # args.dbpath = f"{rundir}/default.pcl"
    args.rundir = rundir
    args.verbosity = "info"
    args.horizontal = False
    args.html = True
    args.ruler = True
    args.decays = True
    args.predictions = True
    args.tex = False
    args.keep = False
    args.commit = False
    if K > 4.0:
        args.commit = True
    plotHiscore.runPlotting ( args )

def loop( rundir : Union[None,os.PathLike] = None,
          maxruns : Union[None,int] = 3, createPlots : bool=True,
          uploadTo : str = "temp", dbpath : str = "official",
          verbose : bool = False, do_srcombine : bool = False ):
    """ loop (maxruns times) that updates hiscores.cache

    :param maxruns: maximally iterate that many times, if None then loop endlessly
    :param createPlots: if False, suppress plotting
    :param uploadTo: upload plots to directory "~/git/smodels.github.io/<uploadTo>"
    :param dbpath: path to database, @rundir@ will replaced with actual rundir
    :param verbose: verbosity
    :param do_srcombine: if we need to reconstruct .hi file, reconstruct the proper
    way.
    """
    rundir = setup( rundir )
    i = 0
    TL, TLold, step, K, Kold = 0., 0., 0, -90., -90.
    TLfile = f"{rundir}/TLold.conf"
    if os.path.exists ( TLfile ):
        with open ( TLfile, "rt" ) as f:
            TLold = float ( f.read().strip() )
    Kfile = f"{rundir}/Kold.conf"
    if os.path.exists ( Kfile ):
        with open ( Kfile, "rt" ) as f:
            Kold = float ( f.read().strip() )
    while True:
        i+=1
        if maxruns != None and i > maxruns:
            break
        D = updateHiscores( rundir, dbpath, do_srcombine )
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
            m = Manipulator ( model )
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
            plot ( TL, K, rundir, uploadTo, dbpath, verbose )
        else:
            print ( "[updateHiscores] was not asked to create plots" )
        time.sleep(60.)
        if os.path.exists ( Kfile ): ## so we can meddle from outside
            with open ( Kfile, "rt" ) as f:
                Kold = float ( f.read().strip() )

if __name__ == "__main__":
    loop()
