#!/usr/bin/env python3

""" simple script that perpetually updates hiscore list
    from H<n>.pcl """

import time, types, sys, os, subprocess

def setup( rundir = None ):
    # codedir = "/mnt/hephy/pheno/ww/git/"
    codedir = "/scratch-cbe/users/wolfgan.waltenberger/git/"
    sys.path.insert(0,"%ssmodels/" % codedir )
    sys.path.insert(0,"%ssmodels-utils/" % codedir )
    sys.path.insert(0,"%s/protomodels/" % codedir )
    sys.path.insert(0,"%ssmodels-utils/prototools/" % codedir )
    if rundir != None:
        rundir = rundir.replace ( "~", os.environ["HOME"] )
        os.chdir ( rundir )
        return rundir
    rundir = "/scratch-cbe/users/wolfgan.waltenberger/rundir/"
    # rundir = "/mnt/hephy/pheno/ww/rundir/"
    # rundir = "./"
    if os.path.exists ( "./rundir.conf" ):
        with open ( "./rundir.conf" ) as f:
            rundir = f.read().strip()
    rundir = rundir.replace ( "~", os.environ["HOME"] )
    os.chdir ( rundir )
    return rundir

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
                    slurmfile = f"/scratch-cbe/users/wolfgan.waltenberger/outputs/walk-{slurmid}.out"
                #print ( nr, laststep )
                steps[nr]=laststep
                if writeSubmitFile and laststep < 1000:
                    rundir = os.getcwd()
                    if rundir.endswith("/"):
                        rundir=rundir[:-1]
                    p = rundir.rfind("/")
                    rundir = rundir[p+1:]
                    g.write ( f"rm -rf /scratch-cbe/users/wolfgan.waltenberger/{rundir}/H{nr}.hi\n" )
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

def updateHiscores( rundir=None ):
    args = types.SimpleNamespace()
    args.print = True
    args.interactive = False
    args.detailed = False
    args.fetch = False
    args.check = False
    args.nmax = 1
    args.dbpath = "%s/default.pcl" % rundir
    args.outfile = "hiscore.hi"
    if rundir != None:
        args.outfile = "%s/hiscore.hi" % rundir
    args.infile = None
    args.rundir = rundir
    # args.maxloss = .01
    # args.nevents = 50000
    from ptools import hiscoreTools
    import socket
    hostname = socket.gethostname().replace(".cbe.vbc.ac.at","")
    print ( "[updateHiscores] now update %s on %s:%s" % \
            ( args.outfile, hostname, rundir ) )
    D = hiscoreTools.main ( args )
    return D

def updateStates( rundir=None):
    args = types.SimpleNamespace()
    args.print = True
    args.rundir = rundir
    args.detailed = False
    args.interactive = False
    args.fetch = False
    args.check = False
    args.dbpath = "%s/default.pcl" % rundir
    args.nmax = 20
    args.outfile = "states.dict"
    if rundir != None:
        args.outfile = "%s/states.dict" % rundir
        args.rundir = rundir
    args.infile = None
    # args.maxloss = .003
    # args.nevents = 50000
    import hiscoreTools
    print ( )
    print ( "[updateHiscores] now update %s" % args.outfile )
    hiscoreTools.main ( args )
    print ( "[updateHiscores] done updating %s" % args.outfile )
    print ( )

def plot( Z, K, rundir, upload="latest" ):
    from plotting import plotHiscore
    from argparse import Namespace
    args = Namespace()
    args.upload = upload
    args.number = 0
    args.detailed = False
    args.destinations = False
    args.picklefile = "%s/hiscore.hi" % rundir
    args.dbpath = "%s/default.pcl" % rundir
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
    if K > 5.0:
        args.commit = True
    plotHiscore.runPlotting ( args )

def main( rundir = None, maxruns=3, doPlots=True, uploadTo="latest" ):
    """ eternal loop that updates hiscore.hi and states.dict
    :param maxruns: maximally iterate that many times
    :param doPlots: if False, suppress plotting
    :param uploadTo: upload plots to directory "~/git/smodels.github.io/<uploadTo>"
    """
    rundir = setup( rundir )
    i = 0
    Z, Zold, step, K, Kold = 0., 0., 0, -90., -90.
    Zfile = "%s/Zold.conf" % rundir
    if os.path.exists ( Zfile ):
        with open ( Zfile, "rt" ) as f:
            Zold = float ( f.read().strip() )
    Kfile = "%s/Kold.conf" % rundir
    if os.path.exists ( Kfile ):
        with open ( Kfile, "rt" ) as f:
            Kold = float ( f.read().strip() )
    while True:
        i+=1
        if maxruns != None and i > maxruns:
            break
        D = updateHiscores( rundir )
        Z,step,model,K = D["Z"],D["step"],D["model"],D["K"]
        if K > Kold + .001:
            from builder.manipulator import Manipulator
            m = Manipulator ( model )
            T=str(int(time.time()))
            m.writeDictFile ( "pmodel-%s.py" % T, comment="history keeper" )
            with open ( "%shistory.txt" % rundir, "at" ) as f:
                f.write ( "%s, step=%d, Z=%.4f, K=%.4f, t=%s\n" % ( time.asctime(),step,Z,K,T) )
                f.close()
            with open ( Zfile, "wt" ) as f:
                f.write ( "%s\n" % str(Z) )
                f.close()
            with open ( Kfile, "wt" ) as f:
                f.write ( "%s\n" % str(K) )
                f.close()
            Zold = Z
            Kold = K
        if doPlots:
            plot ( Z, K, rundir, uploadTo )
        updateStates( rundir )
        time.sleep(60.)
        if os.path.exists ( Kfile ): ## so we can meddle from outside
            with open ( Kfile, "rt" ) as f:
                Kold = float ( f.read().strip() )

if __name__ == "__main__":
    main()
