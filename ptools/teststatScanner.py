#!/usr/bin/env python3

__all__ = [ "TeststatScanner" ]

""" draw K/Z/r as a function of a model parameter,
i.e. scan the various test statistics 
"""

import numpy, sys, os, copy, time, subprocess, glob, pickle
from protomodels.csetup import setup
setup()
from smodels.tools.wrapperBase import WrapperBase
WrapperBase.defaulttempdir="./" ## keep the temps in our folder
from builder.manipulator import Manipulator
from builder.protomodel import ProtoModel
from smodels.tools.runtime import nCPUs
from tester.predictor import Predictor
from ptools.sparticleNames import SParticleNames
from ptools.hiscoreTools import fetchHiscoresObj
from walker.hiscores import Hiscores
from typing import Union, List, Dict
from builder.loggerbase import LoggerBase

class TeststatScanner ( LoggerBase ):
    def __init__ ( self, args ):
        super ( TeststatScanner, self ).__init__ ( 0 )
        self.namer = SParticleNames ( susy = False )
        self.args = vars(args)
        # get the standard rthreshold
        pred = Predictor( 0, do_srcombine = True, dbpath = self.args['dbpath'] )
        self.rthreshold = pred.rthreshold
        self.rundir = setup ( self.args["rundir"] )
        self.nproc = args.nproc
        if self.nproc < 1:
            self.nproc = nCPUs() + self.nproc
        self.hi = fetchHiscoresObj ( )
        self.show_plot = self.args["show_plot"] # dont take this one from pickle

    def massesAreTied ( self, pid1, pid2 ):
        """ are the masses of pid1 and pid2 tied originally? """
        if not pid1 in self.origmasses:
            return False
        if not pid2 in self.origmasses:
            return False
        dm = self.origmasses[pid1] - self.origmasses[pid2]
        if abs(dm)<1e-5:
            return True
        return False

    def setMass ( self, model, pid : int, mass : float ):
        """ set mass of <pid> to <mass> """
        partners = [ ( 1000023, 1000024 ) ]
        model.masses[pid]=mass

        for pair in partners:
            if not pid in pair:
                continue
            for p in pair:
                if p in model.masses and self.massesAreTied ( p, pid ):
                    model.masses[p]=mass

    def createCriticResults ( self, model : ProtoModel ) -> List:
        """ create a list of the analysis ids of the critic

        :returns: list of dictionaries of robs, rexp, anaid
        """
        critics = []
        if not hasattr ( model, "tpList" ):
            return critics
        for p in model.tpList:
            d = { "robs": p["robs"], "rexp": p["rexp"], 
                  "anaid": p["tp"].dataset.globalInfo.id }
            critics.append ( d )
        return critics


    def predProcess ( self, args ):
        """ one thread that computes predictions for masses given in mrange
        """
        i = args["i"]
        self.pprint ( f"#{i}: starting thread" )
        model = args["model"]
        model.walkerid = 100000+10000*i + model.walkerid
        model.loggerid = model.walkerid
        pid = args["pid"]
        predictor = args["predictor"]
        nevents = args["nevents"]
        rundir = args["rundir"]
        mrange = args["mrange"]
        ret = {}
        right_xsecs = copy.deepcopy ( model._stored_xsecs )
        for ctr,m in enumerate(mrange):
            fname = model.createNewSLHAFileName ( prefix = f"s{i}p{pid}{m:.2f}" )
            self.setMass ( model, pid, m )
            if self.args['preserve_xsecs']:
                pass
                ## fake a computation of xsecs, update the masses check
                #model._stored_xsecs = copy.deepcopy ( right_xsecs )
                #model._xsecMasses = dict([[pid,m] for pid,m in model.masses.items()])
            self.pprint ( f"#{i}: start with {ctr}/{len(mrange)}, "\
                    f"m({self.namer.asciiName(pid)})={m:.1f} ({self.args['nevents']} events)" )
            ret[m]={}
            if self.args['dry_run']:
                self.pprint ( f"#{i}: dry-run, not doing anything" )
                ret[m]["K"] = 0.
                ret[m]["Z"] = 0.
                ret[m]["r"] = [ 0, 0, 0 ]
                ret[m]["critics"] = []
            else:
                _ = predictor.predict ( model, keep_predictions = True )
                ret[m]["K"] = model.K
                ret[m]["Z"] = model.Z
                ret[m]["r"] = model.rvalues
                ret[m]["critics"] = self.createCriticResults ( model )
            print ()
            def prettyPrint(v):
                if type(v) in [ float, numpy.float64 ]:
                    return f"{v:.2f}"
                if type(v) in [ list ]:
                    return prettyPrint(v[0])
                return v
            self.pprint ( f"m({self.namer.asciiName(pid)})={m:.1f}: { ', '.join ( f'{k}={prettyPrint(v)}' for k,v in ret[m].items() ) }" )
            self.pprint ()
            # model.delCurrentSLHA()
        return ret

    def ssmProcess ( self, args ):
        """ one thread that computes predictions for ssms given in ssmrange
        """
        i = args["i"]
        # time.sleep(5*i) ## would that help??
        self.pprint ( f"#{i}: starting thread" )
        model = args["model"]
        pids = args["pids"]
        predictor = args["predictor"]
        nevents = args["nevents"]
        ssmrange = args["ssmrange"]
        origssm = args["ssm"]
        model.walkerid = 200000+10000*i + model.walkerid
        model.loggerid = model.walkerid
        model.createNewSLHAFileName ( prefix = f"ssm{i}p{pids[0]}{pids[1]}{origssm:.2f}" )
        if not pids in model.ssmultipliers:
            self.pprint ( f"#{i}: error cannot find pids {str(pids)}" )
            return
        ret = {}
        model.delXSecs()
        # model.predict ( nevents = nevents, recycle_xsecs = True )
        predictor.predict ( model )
        self.pprint ( f"#{i}: before we begin, Z is {model.Z:.3f}" )

        for ctr,ssm in enumerate(ssmrange):
            # fname = model.createNewSLHAFileName ( prefix = f"ssm{i}p{pids[0]}{pids[1]}{ssm:.2f}" )
            ssmold = model.ssmultipliers[pids]
            ma = Manipulator ( model )
            ma.changeSSM ( pids, ssm, recursive = True )
            ma.M.computeXSecs ( nevents = 100000 )
            # model.predict ( nevents = nevents, recycle_xsecs = True )
            predictor.predict ( ma.M, keep_predictions = True )
            self.pprint ( f"#{i}: we change the ssm from {ssmold:.3f} to {ssm:.3f}" )
            self.pprint ( f"#{i}: start with {ctr}/{len(ssmrange)}, ssm={ssm:.2f} ({self.args['nevents']} events)" )
            self.pprint ( f"#{i}:   `- K({ssm:.3f})={ma.M.K:.3f} Z={ma.M.Z:.3f}" )
            # import sys, IPython; IPython.embed( colors = "neutral" ); sys.exit()
            mssm = ma.M.muhat*ssm
            ret[ssm]={ "Z": ma.M.Z, "r": ma.M.rvalues, 
                "K": ma.M.K,"muhat": ma.M.muhat }
            ret[ssm]["critics"] = self.createCriticResults ( ma.M )
            #if False : # os.exists ( fname ):
            #    os.unlink ( fname )

        return ret

    def produce( self, pid : int = 1000022 ):
        """ produce pickle files of mass scans for pid
        """
        self.pid1 = pid
        if type(pid) in [ list, tuple, set ]:
            pid = set(map(abs,pid))
            for p in pid:
                produce ( p )
            return
        pid = abs(pid)
        model = self.hi.hiscores[0]
        self.origmasses = copy.deepcopy ( model.masses )
        fac = self.args["factor"]
        if fac == None:
            fac = 1.008
            if model.masses[pid]<150.:
                fac = 1.007
            if model.masses[pid]<80.:
                fac = 1.006
        if self.args['preserve_xsecs'] and not hasattr ( model, "stored_xsecs" ):
            self.pprint ( "preserve_xsec mode, so computing the xsecs now" )
        model.delXSecs()
        model.nevents = 100000
        model.computeXSecs( nevents = 100000, keep_slha = True )
        if model == None:
            self.pprint ( f"cannot find a model in {self.hi.pickleFile}" )
        apid = abs(pid)
        mass = model.masses[pid]
        if mass > 9e5:
            self.pprint ( f"mass {mass} too high. Wont produce." )
            return
        # model.createNewSLHAFileName ( prefix = f"scan{str(pid)}" )
        values = {}
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
        mranges = [ mrangetot[i::self.nproc] for i in range(self.nproc) ]
        self.pprint ( f"start scanning with m({self.namer.asciiName(pid)})={mass:.1f} with {self.nproc} procs, {len(mrangetot)} mass points, {self.args['nevents']} events, select={self.args['select']}, do_srcombine={self.args['do_srcombine']}" )
        expected = False
        predictor =  Predictor( 0, dbpath=self.args['dbpath'], 
                do_srcombine=self.args['do_srcombine'],
                expected=expected, select=self.args['select'] )
        import multiprocessing
        pool = multiprocessing.Pool ( processes = len(mranges) )
        args = [ { "model": model, "rundir": self.args["rundir"], "pid": pid, "nevents": self.args['nevents'], "predictor": predictor, "dry_run": self.args['dry_run'], "i": i, "mrange": x } for i,x in enumerate(mranges) ]
        values={}
        tmp = pool.map ( self.predProcess, args )
        # model.delCurrentSLHA()
        for r in tmp:
            values.update(r)
        if self.args['dry_run']:
            return
        with open ( f"scanM{self.namer.asciiName(pid)}.pcl", "wb" ) as f:
            pickle.dump ( values, f )
            pickle.dump ( mass, f )
            pickle.dump ( time.asctime(), f )
            pickle.dump ( self.args, f )
            f.close()

    def produceSSMs( self, pid1, pid2 ):
        """ produce pickle files for ssm scan, for (pid1,pid2)
        """
        self.pid1, self.pid2 = pid1, pid2
        fac = self.args["factor"]
        if fac == None:
            fac = 1.008
        self.pprint ( f"produceSSMs {pid1} {pid2}" )
        model = self.hi.hiscores[0]
        origpids = ( pid1, pid2 )
        pids = ( pid1, pid2 )
        if pid2 < pid1:
            pids = ( pid2, pid1 )
        if not pids in model.ssmultipliers:
            self.pprint ( f"could not find pids {str(pids)} in multipliers" )
            self.pprint ( f"only {model.ssmultipliers}" )
            return
        ssm = model.ssmultipliers[pids]
        if ssm == 0.:
            self.pprint ( f"WARNING ssm is zero!!" )
            ssm = 1e-3
        # self.pprint ( f"starting with {str(pids)}: {ssm:.2f}"
        values = {}
        ## meaning we our range goes from ssm_hat * fm to ssm_hat / fm
        fm = .4 ## relative lower bound on ssm
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
        nproc = self.nproc
        if nproc > len(ssmrangetot):
            nproc = len(ssmrangetot)
        ssmranges = [ ssmrangetot[i::nproc] for i in range(nproc) ]
        # ssmranges[0].insert(0,ssm) # compute it early!
        self.pprint ( f"start scanning with ssm({self.namer.asciiName(pid1)},{self.namer.asciiName(pid2)})={ssm:.2f} with {nproc} procs, {len(ssmrangetot)} ssm points, {self.args['nevents']} events" )
            
        expected = False
        predictor =  Predictor( 0, dbpath=self.args["dbpath"], 
                do_srcombine=self.args['do_srcombine'],
                expected=expected, select=self.args['select'] )
        args = [ { "model": model, "pids": pids, "nevents": self.args['nevents'], 
                   "ssm": ssm, "predictor": predictor, "rundir": self.args["rundir"], 
                   "dry_run": self.args['dry_run'],
                   "i": i, "ssmrange": x } for i,x in enumerate(ssmranges) ]
        values={}
        if nproc == 1:
            values = self.ssmProcess ( args[0] )
        else:
            import multiprocessing
            pool = multiprocessing.Pool ( processes = len(ssmranges) )
            tmp = pool.map ( self.ssmProcess, args )
            for r in tmp:
                values.update(r)
        if self.args['dry_run']:
            return
        filename = f"ssm{self.namer.asciiName(origpids[0])}{self.namer.asciiName(origpids[1])}.pcl"
        if os.path.exists ( filename ):
            subprocess.getoutput ( f"cp {filename} {filename}old" )
        with open ( filename, "wb" ) as f:
            pickle.dump ( values, f )
            pickle.dump ( ssm, f )
            pickle.dump ( time.asctime(), f )
            pickle.dump ( self.args, f )
            f.close()

    def findPidPairs ( self ):
        """ search for ssm*pcl files, report the corresponding pid pairs.
        :returns: list of tuple of pids
        """
        ret = []
        files = glob.glob("ssm*pcl")
        files += glob.glob( f"{self.rundir}/ssm*pcl" )
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

    def findPids ( self  ):
        """ search for scanM*pcl files, report the corresponding pids.
        :returns: set of pids
        """
        ret = set()
        files = glob.glob("scanM*pcl")
        files += glob.glob( f"{self.rundir}/scanM*pcl" )
        for f in files:
            p = f.find("scanM")
            s = f[p+5:]
            s = s.replace(".pcl","")
            try:
                ret.add ( int(s) )
            except ValueError as e:
                self.pprint  ( f"when trying to find pids to scan: {e}. will ignore {f}." )
                ## if there are additional files: ignore them.
                pass
        return ret


    def getClosest ( self, key, D, maxDist=1e-5 ):
        """ get the entry in dictionary D whose key is closest to k 
        :returns: D[key] of something close. None if nothing was found within maxDist
        """
        if key in D: ## exact match
            return D[key]
        dmax=maxDist+1e-5
        entry = None
        for k,v in D.items():
            if abs(key-k)<dmax:
                dmax=abs(key-k)
                entry = v
        return v

    def draw( self, pid : int = 1000022, interactive : bool = False, 
              pid2 : int = 0  ) -> Union[str,List]:
        """ draw plots

        :returns: filename of plot(s)
        """
        self.pid1, self.pid2 = pid, pid2
        plotrmax = self.args['plot_rmax']
        if pid2 == 0: ## means all
            pidpairs = self.findPidPairs( )
            if len(pidpairs) == 0:
                self.pprint ( "could not find ssm*pcl files. Maybe run ./fetchFromClip.py --ssms" )
                return
            fnames = []
            for pids in pidpairs:
                try:
                    fname = draw ( pids[0], interactive, pids[1] )
                    fnames.append ( fname )
                except Exception as e:
                    self.pprint ( f"{str(e)}" )
            return fnames

        def isSSMPlot():
            ## is this an ssm or a mass plot
            return pid2!=-1

        import matplotlib
        matplotlib.use("Agg")
        from matplotlib import pyplot as plt
        picklefile = f"{self.rundir}/scanM{self.namer.asciiName(pid)}.pcl"
        if isSSMPlot():
            picklefile = f"{self.rundir}/ssm{self.namer.asciiName(pid)}{self.namer.asciiName(pid2)}.pcl"
        with open ( picklefile, "rb" ) as f:
            values = pickle.load( f )
            cmass = pickle.load ( f ) ## cmass is pids
            timestamp = pickle.load ( f )
            self.args = pickle.load ( f )
            nevents = self.args["nevents"]
        x = list(values.keys()) ## the x's correspond to e.g. the particle masses
        x.sort()
        Kvalues, Kdashed = [], []
        # Kdashed is for regions excluded by critic
        rs = []
        rsarea = []
        for i in x:
            y_ = values[i]
            y0=y_
            y0 = y_["K"]
            if y_["r"][0] > self.rthreshold+.05 and plotrmax:
                rsarea.append ( y_["r"][0] )
                y0 = float("nan")
            else:
                rsarea.append ( 0. )
            rs.append ( y_["r"][0] )
            K = values[i]["K"]
            if K == None:
                K = float("nan")
            Kdashed.append ( K )
            Kvalues.append ( y0 )
        pname = self.namer.texName ( pid, addDollars=False )
        if isSSMPlot():
            pname = self.namer.texName ( pid2, addDollars=False, addSign=True )+\
                    self.namer.texName ( pid, addDollars=False, addSign=True )
        fig,ax1 = plt.subplots()
        plt.plot ( x, Kdashed, linewidth=.3, c="tab:blue", zorder=0 )
        plt.plot ( x, Kvalues, linewidth=2., c="tab:blue", zorder=0, label=f"$K({pname})$" )
        ax1.tick_params ( axis="y", labelcolor="tab:blue", labelsize=12, 
                          labelleft=True )
        ax1.tick_params ( axis="x", labelsize=12 )
        ax1.set_ylabel ( "K", c="tab:blue", fontsize=15 )
        # ax1.set_xlabel ( "m [GeV]", fontsize=13 )
        ax1.set_xlabel ( f"$m({pname})$ [GeV]", fontsize=16 )
        maxyr = numpy.nanmax(Kvalues)
        # ax1.set_ylim ( bottom = 0., top=8.4 )
        #ax1.set_ylim ( bottom = 2., top=maxyr*1.03 )
        rsarea[0]=0.
        rsarea[-1]=0.
        if len(rs) == len(x) and plotrmax:
            ax2 = ax1.twinx()
            ax1.plot ([], [], label="$r_\mathrm{max}$", c="tab:red", zorder=1 )
            ax2.plot ( x, rs, label="$r_\mathrm{max}$", c="tab:red", zorder=2 )
            ax2.tick_params ( axis="y", labelcolor="tab:red", labelsize=12 )
            #ax2.set_ylim ( bottom=min(rs)*.7, top = 1.9 )
            bottom = 0.
            top = min ( 2.1, max(rs)+.2 )
            ax2.set_ylim ( bottom=bottom, top = top )
            nbins = int ( numpy.ceil ( ( top - bottom ) / 0.5 ) )
            # tick only every 0.5
            ax2.locator_params(nbins = nbins )
            ax2.set_ylabel ( "$r_\mathrm{max}$", c="tab:red", fontsize=16 )
        if len(rsarea) == len(x) and plotrmax:
            # ax3 = ax1.twinx()
            ax2.fill ( x, rsarea, lw=0, edgecolor="white", alpha=.2, 
                       facecolor="tab:red", zorder=-1 )
        ymax = numpy.nanmax(Kvalues)
        imax = Kvalues.index ( ymax )
        xmax = x[imax]
        param=f"{xmax} GeV"
        if isSSMPlot():
            param=f"{xmax:.3f}"
        label = f"maximum K\n K({xmax:.2f})={ymax:.2f}"
        #label = "maximum K"
        ax1.scatter ( [ xmax ], [ ymax ], label=label, s=130, c="k", marker="v", zorder=5 )
        if type(cmass)==tuple:
            cmass = x[int(len(x)/2)]
        param = f"{cmass} GeV"
        if isSSMPlot():
            param= f"{cmass:.3f}" 
        Zmax = self.getClosest( cmass, values )
        if type(Zmax)==dict:
            Zmax=Zmax["K"]
        # label = f"proto-model\n K({param})={Zmax:.2f}"
        label = f"proto-model\nK({cmass:.2f})={Zmax:.2f}"
        self.pprint ( f"Zmax=Z({cmass:.3f})={Zmax:.3f}" )
        self.pprint ( f"r({cmass})={values[cmass]['r'][0]:.2f}:{values[cmass]['critics'][0]}" )
        ax1.scatter ( [ cmass ], [ Zmax ], label=label, marker="^", s=130, c="g", zorder=10 )
        plotCriticAtMax = True
        if plotCriticAtMax:
            #ratmax_text = f"r({cmass:.2f})=\n{values[cmass]['critics'][0]}"
            # print ( "critics", values[cmass]["r"] )
            # print ( "critics", values[cmass]["critics"] )
            def findStrongestCritic ( critics : List ):
                """ in the list of critics, find the one with highest robs """
                rmax, cname = -1., ""
                for critic in critics:
                    if critic["robs"]>rmax:
                        cname = critic["anaid"]
                        rmax = critic["robs"]
            #  import sys, IPython; IPython.embed( colors = "neutral" ); sys.exit()
            ratmax_text = f"{values[cmass]['critics'][0]['anaid']}"
            ratmax_text = ratmax_text.replace("-agg","")
            p1 = ratmax_text.find(":")
            if p1 > 0:
                ratmax_text = ratmax_text[:p1]
            ax1.scatter ( [], [],
                   label = ratmax_text, s = 100, marker="*", c="tab:red", zorder=5 )
            ax2.scatter ( [ cmass ], [ values[cmass]['r'][0] ],
                   label = ratmax_text, s = 100, marker="*", c="tab:red", zorder=5 )
        # plt.title ( f"Test statistic $K=K({pname})$", fontsize=14 )
        if not self.args["notimestamp"]:
            plt.text ( .72, -.14, timestamp, c="gray", transform = ax1.transAxes )
        if isSSMPlot():
            plt.xlabel ( f"$\\hat\\mu\\kappa({pname})$", fontsize=17 )
            ax1.set_xlabel ( f"$\\hat\\mu\\kappa({pname})$", fontsize=17 )
            ax1.legend( fontsize = 12, loc= "lower left" )
        else:
            ax1.legend( fontsize = 12 )
            plt.xlabel ( f"$m({pname})$ [GeV]", fontsize=16 )

        figname = f"M{self.namer.asciiName(pid)}.png"
        if isSSMPlot():
            figname = f"ssm_{self.namer.asciiName(pid)}_{self.namer.asciiName(pid2)}.png"
        stdvar =  numpy.std ( Kvalues )

        if interactive:
            import IPython
            IPython.embed( using=False )

        if stdvar < 1e-10:
            self.pprint ( f"standard deviation is a {stdvar:.2f}. Not plotting." )
            # return

        self.pprint ( f"creating {figname}" )
        plt.tight_layout()
        plt.savefig ( figname )
        plt.close()
        if self.args["copy"]:
            dest = os.path.expanduser ( "~/git/smodels.github.io" )
            cmd = f"cp {self.rundir}/{figname} {dest}/protomodels/{self.args['uploadTo']}/"
            o = subprocess.getoutput ( cmd )
            self.pprint ( f"{cmd}: {o}" )
        if self.show_plot:
            cmd = f"see {self.rundir}/{figname}"
            o = subprocess.getoutput ( cmd )
            self.pprint ( f"{cmd}: {o}" ) 
        return figname

if __name__ == "__main__":
    import argparse
    argparser = argparse.ArgumentParser(
            description='script that scans/plots various test statistics such as K,Z,r' )
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
            help='multiplication factor [1.008]',
            type=float, default=None )
    argparser.add_argument ( '-e', '--nevents',
            help='number of events [100000]',
            type=int, default=100000 )
    argparser.add_argument ( '-R', '--rundir',
            help='override the default rundir [None]',
            type=str, default=None )
    argparser.add_argument ( '-S', '--select',
            help='select ["all"]',
            type=str, default="all" )
    argparser.add_argument ( '--dbpath',
            help='database path [official]',
            type=str, default="official" )
    argparser.add_argument ( '-P', '--produce',
            help='produce the pickle file',
            action="store_true" )
    argparser.add_argument ( '-D', '--dry_run',
            help='dry_run, dont produce', action="store_true" )
    argparser.add_argument ( '-d', '--draw',
            help='produce the plot',
            action="store_true" )
    argparser.add_argument ( '-r', '--plot_rmax',
            help='plot also rmax', action="store_true" )
    argparser.add_argument ( '-s', '--show_plot',
            help='show the plot, after plotting', action="store_true" )
    argparser.add_argument ( '-I', '--interactive',
            help='interactive mode, starts ipython (only works with -d, and not in bulk mode)',
            action="store_true" )
    argparser.add_argument ( '--preserve_xsecs',
            help='when scanning masses, adjust the SSMs so that the xsecs are constant',
            action="store_true" )
    argparser.add_argument ( '-F', '--force_copy',
            help='force copying the hiscores.cache file',
            action="store_true" )
    argparser.add_argument ( '-u', '--uploadTo',
            help='choose upload directory [latest]',
            type=str, default="latest" )
    argparser.add_argument ( '-c', '--copy',
            help='copy plots to ~/git/smodels.github.io/protomodels/<upload>/',
            action="store_true" )
    argparser.add_argument ( '-N', '--notimestamp',
            help='dont put a timestamp on it',
            action="store_true" )
    argparser.add_argument ( '--do_srcombine',
            help='do srombine',
            action="store_true" )
    args = argparser.parse_args()
    scanner = TeststatScanner( args )
    pids = args.pid
    if pids == 0:
        allpids = scanner.findPids( )
        pids = allpids
    if args.produce:
        if args.pid2 > 0:
            scanner.produceSSMs( args.pid, args.pid2 )
        else:
            scanner.produce( pids )
    if args.draw:
        if args.pid != 0:
            scanner.draw( pids, args.interactive, pid2 = args.pid2 )
        else:
            allpids = scanner.findPids( )
            for pid in allpids:
                try:
                    scanner.draw( pid, args.interactive, args.pid2  )
                except Exception as e:
                    self.pprint ( f"skipping {pid}: {e}" )
