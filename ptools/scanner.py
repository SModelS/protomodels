#!/usr/bin/env python3

""" draw Z as a function of a model parameter """

import numpy, sys, os, copy, time, subprocess, glob
from csetup import setup
setup()
from smodels.tools.wrapperBase import WrapperBase
WrapperBase.defaulttempdir="./" ## keep the temps in our folder
from builder.manipulator import Manipulator
from smodels.tools.runtime import nCPUs
from tester.predictor import Predictor
from ptools.sparticleNames import SParticleNames

def getHiscore( force_copy = False, rundir = None ):
    """ get the hiscore from the picklefile
    :param force_copy: if True, force a cp command on the pickle file
    """
    from walker.hiscore import Hiscore
    # spids = str(pids).replace("[","").replace("]","").replace(" ","").replace(",","").replace("0","")
    picklefile =rundir + "hiscore2.hi" # % spids
    backupfile = rundir+"hiscore.hi"
    # picklefile =rundir + "hiscore.hi" # % spids
    ## do this always
    h2Outdated = False
    if os.path.exists ( picklefile ) and os.path.exists ( backupfile ):
        if os.stat ( picklefile ).st_mtime < os.stat ( backupfile ).st_mtime:
            h2Outdated = True
    if force_copy or (not os.path.exists ( picklefile )) or h2Outdated:
        cmd = "cp %s %s" % ( backupfile, picklefile )
        import subprocess
        o = subprocess.getoutput ( cmd )
        print ( "[scanner] %s: %s" % ( cmd, o ) )
    import socket
    hostname = socket.gethostname().replace(".cbe.vbc.ac.at","")
    print ( "[scanner] retrieving hiscore object %s on %s .... " % \
             ( picklefile, hostname ) )
    hi = Hiscore( walkerid=0, save_hiscores=False, picklefile = picklefile )
    Z=hi.hiscores[0].Z
    K=hi.hiscores[0].K
    print ( "[scanner] done retrieving hiscore object, highest at K=%.2f, Z=%.2f" % \
             (K, Z ) )
    return hi

def predProcess ( args ):
    """ one thread that computes predictions for masses given in mrange
    """
    i = args["i"]
    import time
    # time.sleep(5*i) ## would that help??
    print ( "[scanner:%d] starting thread" % ( i ) )
    model = args["model"]
    model.walkerid = 100000+10000*i + model.walkerid
    pid = args["pid"]
    predictor = args["predictor"]
    nevents = args["nevents"]
    dry_run = args["dry_run"]
    rundir = args["rundir"]
    mrange = args["mrange"]
    preserve_xsecs = args["preserve_xsecs"]
    ret = {}
    right_xsecs = copy.deepcopy ( model._stored_xsecs )
    for ctr,m in enumerate(mrange):
        model.createNewSLHAFileName ( prefix = "s%dp%d%.2f" % ( i, pid, m ) )
        model.masses[pid]=m
        if preserve_xsecs:
            ## fake a computation of xsecs, update the masses check
            model._stored_xsecs = copy.deepcopy ( right_xsecs )
            model._xsecMasses = dict([[pid,m] for pid,m in model.masses.items()])
        ts = time.strftime("%H:%M:%S" )
        print ( "[scanner:%d-%s] start with %d/%d, m=%.1f (%d events)" % \
                ( i, ts, ctr, len(mrange), m, nevents ) )
        #model.predict ( nevents = nevents, check_thresholds=False )
        if dry_run:
            print ( "[scanner:%d-%s] dry-run, not doing anything" % ( i, ts ) )
            model.K = 0.
            model.Z = 0.
            model.rvalues = [ 0, 0, 0 ]
        else:
            predictor.predict ( model ) # , nevents = nevents, check_thresholds=False )
        ret[m]=(model.Z,model.rvalues[0],model.K)
        # model.delCurrentSLHA()
    return ret

def printCombo ( combo, comment="" ):
    """ pretty print a theory pred combo """
    print ( "combination %s" % comment )
    for tp in combo:
        print( " `- %s" % tp.analysisId() )

def printXSecs ( xsecs, comment="" ):
    """ pretty print the list of xsecs """
    print ( "xsecs %s" % comment )
    for xsec in xsecs:
        print ( " `- %s, %s [%.30s]" % ( xsec.info, xsec.value, xsec.pid ) )

def printModel ( model, comment="" ):
    print ( "model %s" % comment )
    txt=""
    for pid,m in model.masses.items():
        if m>5e5:
            continue
        txt+="%d:%d, " % ( pid, m )
    print ( " `- m %s" % txt[:-2] )
    txt=""
    for pids,ssm in model.ssmultipliers.items():
        if abs(ssm-1)<1e-2:
            continue
        txt+="%s:%.2f, " % ( pids, ssm )
    print ( " `- ssm %s" % txt[:-2] )


def ssmProcess ( args ):
    """ one thread that computes predictions for ssms given in ssmrange
    """
    i = args["i"]
    import time
    # time.sleep(5*i) ## would that help??
    print ( "[scanner:%d] starting thread" % ( i ) )
    model = args["model"]
    pids = args["pids"]
    predictor = args["predictor"]
    nevents = args["nevents"]
    ssmrange = args["ssmrange"]
    ssm = args["ssm"]
    model.walkerid = 200000+10000*i + model.walkerid
    model.createNewSLHAFileName ( prefix = "ssm%dp%d%d%.2f" % ( i, pids[0], pids[1], ssm ) )
    if not pids in model.ssmultipliers:
        print ( "[scanner:%d] error cannot find pids %s" % (i, str(pids) ) )
        return
    ret = {}
    ts = time.strftime("%H:%M:%S" )
    model.delXSecs()
    # model.predict ( nevents = nevents, recycle_xsecs = True )
    predictor.predict ( model )
    print ( "[scanner:%d-%s] before we begin, Z is %.3f" % ( i, ts, model.Z ) )

    for ctr,ssm in enumerate(ssmrange):
        ssmold = model.ssmultipliers[pids]
        print ( "[scanner:%d] we change the ssm from %.3f to %.3f" % \
                ( i, ssmold, ssm ) )
        ma = Manipulator ( model )
        ma.changeSSM ( pids, ssm )
        model = ma.M
        ts = time.strftime("%H:%M:%S" )
        print ( "[scanner:%d-%s] start with %d/%d, ssm=%.2f (%d events)" % \
                ( i, ts, ctr, len(ssmrange), ssm, nevents ) )
        # model.predict ( nevents = nevents, recycle_xsecs = True )
        predictor.predict ( model ) # #nevents = nevents, recycle_xsecs = True )
        print ( "[scanner:%d-%s]   `- Z=%.3f" % ( i, ts, model.Z ) )
        ret[ssm]=(model.Z,model.rvalues[0],model.K)
    return ret

def produce( hi, pid=1000022, nevents = 100000, dry_run=False,
             nproc=5, fac = 1.008, rundir = "", preserve_xsecs = False ):
    """ produce pickle files of mass scans for pid, with nevents
    :param hi: hiscore list object
    :param nproc: number of processes
    :param fac: factor with which to multiply interval
    :param preserve_xsecs: adjust the SSMs to that xsecs are preserved
    """
    if type(pid) in [ list, tuple, set ]:
        pid = set(map(abs,pid))
        for p in pid:
            produce ( hi, p, nevents, dry_run, nproc, fac, rundir = rundir,
                      preserve_xsecs = preserve_xsecs )
        return
    pid = abs(pid)
    model = hi.hiscores[0]
    if fac == None:
        fac = 1.008
        if model.masses[pid]<150.:
            fac = 1.007
        if model.masses[pid]<80.:
            fac = 1.006
    if preserve_xsecs and not hasattr ( model, "stored_xsecs" ):
        print ( "[scanner] preserve_xsec mode, so computing the xsecs now" )
        model.computeXSecs( keep_slha = True )
    if model == None:
        print ( "[scanner] cannot find a model in %s" % hi.pickleFile )
    apid = abs(pid)
    mass = model.masses[pid]
    if mass > 9e5:
        print ( "mass %d too high. Wont produce." % mass )
        return
    # model.createNewSLHAFileName ( prefix = "scan%s" % pid )
    Zs = {}
    fm = .6 ## lower bound (relative) on mass
    # mrange = numpy.arange ( mass * fm, mass / fm, .008*mass )
    mrangetot = [ mass ]
    m1,m2 = mass, mass
    dm = fac
    while m1 > fm * mass:
        m1 = mass/dm
        m2 = mass*dm
        dm = dm * fac
        mrangetot.append( m1 )
        if pid in [ 1000006, 2000006 ] and m2 > 1550.:
            continue ## there is nothing to see beyond
        mrangetot.append( m2 )
    mrangetot.sort()
    mranges = [ mrangetot[i::nproc] for i in range(nproc) ]
    print ( "[scanner] start scanning with m(%d)=%.1f with %d procs, %d mass points, %d events" % \
            ( pid, mass, nproc, len(mrangetot), nevents ) )
    expected = False
    select = "all"
    dbpath = rundir + "/default.pcl"
    predictor =  Predictor( 0, dbpath=dbpath,
                            expected=expected, select=select )
    import multiprocessing
    pool = multiprocessing.Pool ( processes = len(mranges) )
    args = [ { "model": model, "rundir": rundir, "pid": pid, "nevents": nevents, "predictor": predictor, "preserve_xsecs": preserve_xsecs, "dry_run": dry_run,
               "i": i, "mrange": x } for i,x in enumerate(mranges) ]
    Zs={}
    tmp = pool.map ( predProcess, args )
    # model.delCurrentSLHA()
    for r in tmp:
        Zs.update(r)
    if dry_run:
        return
    import pickle
    with open ( "scanM%s.pcl" % pid, "wb" ) as f:
        pickle.dump ( Zs, f )
        pickle.dump ( mass, f )
        pickle.dump ( nevents, f )
        pickle.dump ( time.asctime(), f )
        pickle.dump ( preserve_xsecs, f )
        f.close()

def produceSSMs( hi, pid1, pid2, nevents = 100000, dry_run=False,
             nproc=5, fac = 1.008, rundir= "" ):
    """ produce pickle files for ssm scan, for (pid1,pid2), with nevents
    :param hi: hiscore list object
    :param nproc: number of processes
    :param fac: factor with which to multiply interval
    """
    if fac == None:
        fac = 1.008
    print ( "produceSSMs", pid1, pid2 )
    model = hi.hiscores[0]
    pids = ( pid1, pid2 )
    if pid2 < pid1:
        pids = ( pid2, pid1 )
    if not pids in model.ssmultipliers:
        print ( "[scanner] could not find pids %s in multipliers" % ( str(pids) ) )
        print ( "only", model.ssmultipliers )
        return
    ssm = model.ssmultipliers[pids]
    # print ( "[scanner] starting with %s: %.2f" % ( pids, ssm ) )
    Zs = {}
    fm = .6 ## lower bound (relative) on mass
    # mrange = numpy.arange ( ssm * fm, ssm / fm, .01*ssm )
    ssmrangetot = [ ssm ]
    ssm1,ssm2 = ssm, ssm
    dssm = fac
    while ssm1 > fm * ssm:
        ssm1 = ssm/dssm
        ssm2 = ssm*dssm
        ssmrangetot.append( ssm1 )
        ssmrangetot.append( ssm2 )
        dssm = dssm * fac
    ssmrangetot.sort()
    if nproc > len(ssmrangetot):
        nproc = len(ssmrangetot)
    ssmranges = [ ssmrangetot[i::nproc] for i in range(nproc) ]
    print ( "[scanner] start scanning with ssm(%d,%d)=%.2f with %d procs, %d ssm points, %d events" % \
            ( pid1, pid2, ssm, nproc, len(ssmrangetot), nevents ) )
    import multiprocessing
    pool = multiprocessing.Pool ( processes = len(ssmranges) )
    expected = False
    select = "all"
    dbpath = rundir + "/default.pcl"
    predictor =  Predictor( 0, dbpath=dbpath,
                            expected=expected, select=select )
    args = [ { "model": model, "pids": pids, "nevents": nevents, "ssm": ssm,
               "predictor": predictor, "rundir": rundir, "dry_run": dry_run,
               "i": i, "ssmrange": x } for i,x in enumerate(ssmranges) ]
    Zs={}
    tmp = pool.map ( ssmProcess, args )
    for r in tmp:
        Zs.update(r)
    if dry_run:
        return
    import pickle
    filename = "ssm%d%d.pcl" % (pids[0],pids[1])
    if os.path.exists ( filename ):
        subprocess.getoutput ( "cp %s %sold" % ( filename, filename ) )
    with open ( filename, "wb" ) as f:
        pickle.dump ( Zs, f )
        pickle.dump ( ssm, f )
        pickle.dump ( nevents, f )
        pickle.dump ( time.asctime(), f )
        f.close()

def findPidPairs ( rundir ):
    """ search for ssm*pcl files, report the corresponding pid pairs.
    :returns: list of tuple of pids
    """
    ret = []
    files = glob.glob("ssm*pcl")
    files += glob.glob("%s/ssm*pcl" % rundir )
    for f in files:
        p = f.find("ssm")
        s = f[p+3:]
        s = s.replace(".pcl","")
        split = 7
        if s[0]=="-":
            split = 8
        pids = ( ( int(s[:split]), int(s[split:]) ) )
        ret.append ( pids )
    return ret

def findPids ( rundir ):
    """ search for scanM*pcl files, report the corresponding pids.
    :returns: set of pids
    """
    ret = set()
    files = glob.glob("scanM*pcl")
    files += glob.glob("%s/scanM*pcl" % rundir )
    for f in files:
        p = f.find("scanM")
        s = f[p+5:]
        s = s.replace(".pcl","")
        ret.add ( int(s) )
    return ret

def draw( pid= 1000022, interactive=False, pid2=0, copy=False,
          drawtimestamp = True, rundir = None, plotrmax=False,
          rthreshold = 1.3, upload = "latest" ):
    """ draw plots
    :param copy: copy final plots to ../../smodels.github.io/protomodels/latest
    :param drawtimestamp: if True, put a timestamp on it
    :param plotrmax: if True, plot also rmax curve
    :param upload: upload directory, default is "latest"
    """
    if pid2 == 0: ## means all
        pidpairs = findPidPairs( rundir )
        if len(pidpairs) == 0:
            print ( "[scanner] could not find ssm*pcl files. Maybe run ./fetchFromClip.py --ssms" )
            return
        for pids in pidpairs:
            try:
                draw ( pids[0], interactive, pids[1], copy, drawtimestamp, \
                       rundir, plotrmax, upload = upload )
            except Exception as e:
                print ( "[scanner] %s" % e )
        return

    def isSSMPlot():
        ## is this an ssm or a mass plot
        return pid2!=-1

    import matplotlib
    matplotlib.use("Agg")
    from matplotlib import pyplot as plt
    import pickle
    namer = SParticleNames ( susy = False )
    #if False:
    #    rundir = ""
    picklefile = "%sscanM%s.pcl" % (rundir, pid )
    if isSSMPlot():
        picklefile = "%sssm%s%d.pcl" % ( rundir, pid, pid2 )
    with open ( picklefile, "rb" ) as f:
        Zs = pickle.load( f )
        cmass = pickle.load ( f ) ## cmass is pids
        nevents = pickle.load ( f )
        timestamp = pickle.load ( f )
    x = list(Zs.keys())
    x.sort()
    y, yr, ydashed = [], [], []
    rs = []
    rsarea = []
    for i in x:
        y_ = Zs[i]
        y0=y_
        if type(y_)==tuple:
            idx = 2
            if len(y_)==2:
                idx = 0
            y0 = y_[idx]
            if y_[1] > rthreshold+.05 and plotrmax:
                rsarea.append ( y_[1] )
                y0 = -1.
            else:
                rsarea.append ( 0. )
            rs.append ( y_[1] )
            ydashed.append ( Zs[i][idx] )
        y2_ = y0
        if y2_ < 0.:
            y2_ = float("nan")
        y.append ( y0 )
        yr.append ( y2_ )
    pname = namer.texName ( pid, addDollars=False )
    if isSSMPlot():
        #pname = namer.texName ( pid, addDollars=False, addSign=True )+","+\
        #        namer.texName ( pid2, addDollars=False, addSign=True )
        pname = namer.texName ( pid2, addDollars=False, addSign=True )+\
                namer.texName ( pid, addDollars=False, addSign=True )
    fig,ax1 = plt.subplots()
    plt.plot ( x, ydashed, linewidth=.3, c="tab:blue", zorder=0 )
    plt.plot ( x, yr, linewidth=2., label="$K(%s)$" % ( pname ), 
               c="tab:blue", zorder=0 )
    #plt.plot ( x, yr, linewidth=2., label="$K(%s)$, %dk events" % ( pname, nevents/1000 ), 
    #           c="tab:blue", zorder=0 )
    ax1.tick_params ( axis="y", labelcolor="tab:blue", labelsize=12, labelleft=True )
    ax1.tick_params ( axis="x", labelsize=12 )
    ax1.set_ylabel ( "K", c="tab:blue", fontsize=15 )
    # ax1.set_xlabel ( "m [GeV]", fontsize=13 )
    ax1.set_xlabel ( "$m(%s)$ [GeV]" % pname, fontsize=16 )
    maxyr = numpy.nanmax(ydashed)
    # print ( "ydashed", ydashed )
    ax1.set_ylim ( bottom = 2., top=8.4 )
    #ax1.set_ylim ( bottom = 2., top=maxyr*1.03 )
    rsarea[0]=0.
    rsarea[-1]=0.
    if len(rs) == len(x) and plotrmax:
        ax2 = ax1.twinx()
        ax1.plot ([], [], label="$r_\mathrm{max}$", c="tab:red", zorder=1 )
        ax2.plot ( x, rs, label="$r_\mathrm{max}$", c="tab:red", zorder=2 )
        ax2.tick_params ( axis="y", labelcolor="tab:red", labelsize=12 )
        #ax2.set_ylim ( bottom=min(rs)*.7, top = 1.9 )
        ax2.set_ylim ( bottom=0.6, top = 1.9 )
        ax2.set_ylabel ( "$r_\mathrm{max}$", c="tab:red", fontsize=16 )
    if len(rsarea) == len(x) and plotrmax:
        # ax3 = ax1.twinx()
        ax2.fill ( x, rsarea, lw=0, edgecolor="white", alpha=.2, 
                   facecolor="tab:red", zorder=-1 )
    ymax = max(y)
    imax = y.index ( ymax )
    xmax = x[imax]
    param="%d GeV" % xmax
    if isSSMPlot():
        param="%.3f" % xmax
    # label = "maximum K\n K(%s)=%.2f" % (param, ymax )
    label = "maximum K"
    ax1.scatter ( [ xmax ], [ ymax ], label=label, s=130, c="k", marker="v", zorder=5 )
    if type(cmass)==tuple:
        cmass = x[int(len(x)/2)]
    param = "%d GeV" % cmass
    if isSSMPlot():
        param="%.3f" % cmass
    Zmax = Zs[cmass]
    if type(Zmax)==tuple:
        Zmax=Zmax[idx]
    # label = "proto-model\n K(%s)=%.2f" % (param, Zmax )
    label = "proto-model"
    ax1.scatter ( [ cmass ], [ Zmax ], label=label, marker="^", s=130, c="g", zorder=10 )
    # plt.title ( "Test statistic $K=K(%s)$" % pname, fontsize=14 )
    if drawtimestamp:
        plt.text ( .7, -.12, timestamp, c="gray", transform = ax1.transAxes )
    if isSSMPlot():
        plt.xlabel ( "$\\hat\\mu\\kappa(%s)$" % pname, fontsize=17 )
        ax1.set_xlabel ( "$\\hat\\mu\\kappa(%s)$" % pname, fontsize=17 )
        ax1.legend( fontsize = 12, loc= "lower left" )
    else:
        ax1.legend( fontsize = 12 )
        plt.xlabel ( "$m(%s)$ [GeV]" % pname, fontsize=16 )

    # plt.text ( .9*min(x)+.1*(max(x)-min(x)), 1.*max(y), "%d events" % nevents )
    figname = "M%d.png" % pid
    if isSSMPlot():
        figname = "ssm_%d_%d.png" % ( pid, pid2 )
    stdvar =  numpy.std ( y )

    if interactive:
        import IPython
        IPython.embed( using=False )

    if stdvar < 1e-10:
        print ( "[scanner] standard deviation is a %.2f. Not plotting." % stdvar )
        return

    print ( "[scanner] creating %s" % figname )
    plt.tight_layout()
    plt.savefig ( figname )
    plt.close()
    if copy:
        dest = os.path.expanduser ( "~/git/smodels.github.io" )
        cmd = "cp %s/%s %s/protomodels/%s/" % ( rundir,figname, dest, upload )
        o = subprocess.getoutput ( cmd )
        print ( "[scanner] %s: %s" % ( cmd, o ) )

if __name__ == "__main__":
    import argparse
    argparser = argparse.ArgumentParser(
            description='script that takes care of the Z(m) plots' )
    argparser.add_argument ( '-p', '--pid', '--pid1',
            help='pid to consider. If zero, then consider a predefined list [0]',
            type=int, default=0 )
    argparser.add_argument ( '-q', '--pid2',
            help='pid 2. if -1, then scan masses, If not, then scan signal strength multipliers. If zero, then scan all ssms [-1]',
            type=int, default=-1 )
    argparser.add_argument ( '-n', '--nproc',
            help='number of processes, if zero then determine automatically [0]',
            type=int, default=0 )
    argparser.add_argument ( '-f', '--factor',
            help='multiplication factor [None=1.008]',
            type=float, default=None )
    argparser.add_argument ( '-e', '--nevents',
            help='number of events [100000]',
            type=int, default=100000 )
    argparser.add_argument ( '-R', '--rundir',
            help='override the default rundir [None]',
            type=str, default=None )
    argparser.add_argument ( '-P', '--produce',
            help='produce the pickle file',
            action="store_true" )
    argparser.add_argument ( '-D', '--dry_run',
            help='dry_run, dont produce',
            action="store_true" )
    argparser.add_argument ( '-d', '--draw',
            help='produce the plot',
            action="store_true" )
    argparser.add_argument ( '-I', '--interactive',
            help='interactive mode, starts ipython (only works with -d, and not in bulk mode)',
            action="store_true" )
    argparser.add_argument ( '--preserve_xsecs',
            help='when scanning masses, adjust the SSMs so that the xsecs are constant',
            action="store_true" )
    argparser.add_argument ( '-F', '--force_copy',
            help='force copying the hiscore.hi file',
            action="store_true" )
    argparser.add_argument ( '-u', '--upload',
            help='choose upload directory [latest]',
            type=str, default="latest" )
    argparser.add_argument ( '-c', '--copy',
            help='copy plots to ~/git/smodels.github.io/protomodels/<upload>/',
            action="store_true" )
    argparser.add_argument ( '-N', '--notimestamp',
            help='dont put a timestamp on it',
            action="store_true" )
    args = argparser.parse_args()
    drawtimestamp = not args.notimestamp
    rundir = setup( args.rundir )
    nproc = args.nproc
    if nproc < 1:
        nproc = nCPUs() + nproc
    allpids = findPids( rundir )
    pids = args.pid
    if pids == 0:
        pids = allpids
    if args.produce:
        hi = getHiscore( args.force_copy, rundir )
        if args.pid2 > 0:
            produceSSMs( hi, args.pid, args.pid2, args.nevents, args.dry_run, nproc, args.factor, rundir = rundir )
        else:
            produce( hi, pids, args.nevents, args.dry_run, nproc, args.factor, rundir = rundir, preserve_xsecs = args.preserve_xsecs )
    pred = Predictor( 0 )
    rthreshold = pred.rthreshold
    if args.draw:
        if args.pid != 0:
            draw( pids, args.interactive, args.pid2, args.copy, drawtimestamp, rundir, \
                  plotrmax = False, rthreshold = rthreshold, upload = args.upload )
        else:
            for pid in allpids:
                try:
                    draw( pid, args.interactive, args.pid2, args.copy, drawtimestamp,
                          rundir, plotrmax = False, rthreshold = rthreshold, 
                          upload = args.upload )
                except Exception as e:
                    print ( "[scanner] skipping %d: %s" % ( pid, e ) )
