#!/usr/bin/env python3

""" script used to produce the likelihood scans """

__all__ = [ "NLLScanner" ]

import os, sys, multiprocessing, time, numpy, subprocess, copy, glob
import pickle, random, shutil
import numpy as np
from typing import Dict, Tuple, Union, List

try:
    from csetup import setup
    setup()
except ModuleNotFoundError as e:
    pass
from smodels.tools.wrapperBase import WrapperBase
WrapperBase.defaulttempdir="./" ## keep the temps in our folder
from smodels.base.physicsUnits import fb
from smodels.base.runtime import nCPUs
from smodels.matching.theoryPrediction import TheoryPrediction
from smodels.statistics.basicStats import observed, apriori, aposteriori,\
         NllEvalType

from base.loggerbase import LoggerBase
from base.runEnviron import RunEnviron

from tester.combiner import Combiner
from tester.predictor import Predictor
from tester.critic import Critic
from ptools.sparticleNames import SParticleNames
from ptools import moreHelpers, helpers
from ptools.helpers import py_dumps

namer = SParticleNames ( False )
t0 = time.time() ## define t0, to measure how long things took

def findPids ( rundir ):
    """ search for nll*pcl files, report the corresponding pids.
    :returns: set of pids
    """
    ret = set()
    files = glob.glob("nll*pcl")
    files += glob.glob( f"{rundir}/nll*pcl" )
    for f in files:
        p = f.find("nll")
        s = f[p+4:]
        s = s.replace(".pcl","")
        s = s.replace("1000022","")
        s = s.replace("X1Z","")
        ret.add ( int(s) )
    print ( f"[NLLScanner] pids are {ret}" )
    return ret

class NLLThread ( LoggerBase ):
    """ one thread of the sweep """
    def __init__ ( self, threadnr: str, obj ):
        """ the constructor.
        """
        super ( NLLThread, self ).__init__ ( threadnr )
        self.environ = obj.environ
        self.resultsdir = obj.resultsdir
        self.topo = obj.topo
        self.threadnr = threadnr
        self.dict_file = obj.dict_file
        self.picklefile = obj.picklefile
        self.M = copy.deepcopy ( obj.M )
        self.origmasses = copy.deepcopy ( self.M.masses )
        self.origssmultipliers = copy.deepcopy ( self.M.ssmultipliers )
        self.M.createNewSLHAFileName ( prefix=f"lthrd{self.threadnr}_{obj.xvariable}" )
        self.xvariable = obj.xvariable
        self.yvariable = obj.yvariable
        self.mxvariable = obj.mxvariable
        self.myvariable = obj.myvariable
        self.nevents = obj.nevents
        self.predictor = obj.predictor
        self.critic = obj.critic
        helpers.mkdir ( self.resultsdir )

    def getDefaultDictionary ( self ):
        """ initialise the dictionary for the pickle file """
        import time
        d = { "masspoints": [], "mxvariable": self.mxvariable,
              "myvariable": self.myvariable, "nevents": self.nevents,
              "topo": self.topo, "timestamp": time.asctime(),
              "xvariable": self.xvariable, "yvariable": self.yvariable,
              "model": self.M.dict() }
        return d

    def getMetaInformation ( self ) -> dict:
        """ the meta information about this nll scan.
        contains command line arguments, date of production, 
        x, y variables.
        """ 
        from smodels_utils.helper.various import getCommandLine
        meta = { "cmdline": getCommandLine() }
        from datetime import datetime
        from zoneinfo import ZoneInfo
        now = datetime.now(ZoneInfo("Europe/Vienna"))
        meta["created"]=now.isoformat()
        import socket
        hostname = socket.gethostname()
        meta["hostname"]=socket.gethostname()
        meta["dt[h]"]=(time.time()-t0)/60./60. # time it took in hours
        meta["xvariable"]=self.xvariable
        meta["yvariable"]=self.yvariable
        return meta

    def createPickleBackup ( self ):
        if os.path.exists ( self.picklefile ) and \
                os.stat ( self.picklefile ).st_size > 1000:
            try:
                f = open ( self.picklefile, "rb" )
                d = pickle.load(f) ## ok can read. copy, then!
                cmd = f"cp {self.picklefile} {self.picklefile}.old"
                subprocess.getoutput ( cmd )
            except Exception as e:
                pass

    def writePickleFile ( self, d : Dict ):
        """ write the dictionary into the picklefile """
        # self.createPickleBackup()
        f = open ( self.picklefile, "wb" )
        pickle.dump ( d, f )
        f.close()
        writeDictFile = False
        if self.dict_file:
            dictfile = self.picklefile.replace(".pcl",".dict")
            self.pprint ( f"writing to {dictfile}" )
            with open ( dictfile, "wt" ) as f:
                befores = { "mx": "my", "my": "critic", "critic": "oul",
                            "oul": "eul", "eul": "nll" }
                d = py_dumps ( d, level = 0, a_before = befores )
                f.write ( d )
            f.close()

    def lockPickleFile ( self ):
        """ make sure we write sequentially """
        lockfile = f"{self.picklefile}.lock"
        ctr = 0
        while os.path.exists ( lockfile ):
            ctr+=1
            time.sleep ( .1*ctr )
            if ctr > 16:
                self.unlockPickleFile()
                return
        from pathlib import Path
        Path ( lockfile ).touch()

    def unlockPickleFile ( self ):
        lockfile = f"{self.picklefile}.lock"
        if os.path.exists ( lockfile ):
            try:
                os.unlink ( lockfile )
            except FileNotFoundError as e:
                pass

    def isSameMassPoint ( self, point1 : Dict, point2 : Dict ) -> bool:
        """ are the two points identical? """
        dx = point1["mx"]-point2["mx"]
        dy = point1["my"]-point2["my"]
        if abs(dx)<1e-10 and abs(dy)<1e-10:
            return True
        return False

    def writeRunMeta ( self ):
        with open ( f"{self.resultsdir}/run.meta", "wt" ) as f:
            f.write ( f"{{ 'timestamp': '{time.asctime()}', 'ntotal': {self.ntotal} }}\n" )
            f.close()

    def getDictFileName ( self, mx : float, my : float ) -> str:
        """ the dict file name of the point """
        return f"{self.resultsdir}/{mx:.3f}_{my:.3f}.dict"

    def addNewPoint ( self, point : Dict ):
        """ add point to resultsdir. if already in, replace """
        if not "mx" in point: #dont add anything!
            return
        point["mx"]=round( point["mx"], 7 )
        point["my"]=round( point["my"], 7 )
        dictfile = self.getDictFileName ( point["mx"], point["my"] )
        with open ( dictfile, "wt" ) as f:
            # f.write ( f"{point}\n" )
            befores = { "mx": "my", "my": "critic", "critic": "oul",
                        "oul": "eul", "eul": "nll" }
            d = py_dumps ( point, level = 0, a_before = befores )
            f.write ( d + "\n" )
            f.close()
        nfiles = len ( glob.glob ( f"{self.resultsdir}/*.dict" ) )
        if nfiles % 100 == 0: # update with every 20th entry
            self.updatePickleFile()

    def getAllMassPoints ( self ) -> list:
        """ retrieve all mass points from resultsdir """
        files = glob.glob ( f"{self.resultsdir}/*.dict" )
        masspoints = []
        for fname in files:
            with open ( fname, "rt" ) as h:
                try:
                    d = eval(h.read())
                except (SyntaxError,ValueError,TypeError) as e:
                    self.error ( f"could not read: {fname}: {e}" )
                    sys.exit()
                d["mx"] = float ( d["mx"] )
                d["my"] = float ( d["my"] )
                masspoints.append ( d )
        return masspoints

    def updatePickleFile ( self ):
        """ collect all the entries in resultsdir, and compile them
        into one big pickle file """
        self.pprint ( f"updating {self.picklefile}" )
        self.lockPickleFile()
        Dict = self.getDefaultDictionary()
        files = glob.glob ( f"{self.resultsdir}/*.dict" )
        masspoints = self.getAllMassPoints()
        Dict["masspoints"] = masspoints
        Dict["meta"] = self.getMetaInformation()
        self.writePickleFile ( Dict )
        self.unlockPickleFile()

    def unlinkResultsDir ( self ):
        """ clean up in the end """
        if os.path.exists ( self.resultsdir ):
            shutil.rmtree ( self.resultsdir )

    def massesAreTied ( self, xvariable, yvariable ):
        """ are the masses of xvariable and yvariable tied originally? """
        if not xvariable in self.origmasses:
            return False
        if not yvariable in self.origmasses:
            return False
        dm = self.origmasses[xvariable] - self.origmasses[yvariable]
        if abs(dm)<1e-5:
            return True
        return False

    def getPredictions ( self, recycle_xsecs : bool = True ) -> Dict:
        """ get predictions, return likelihoods

        :param recycle_xsecs: if true, then recycle the cross sections, dont
        recompute
        :returns: a diction with likelihoods ("nll"), critics' responses ("critic"),
        observed ("oul") and expected ("eul") upper limits on mu.
        """
        self.debug ( f"asking for predictions for xmy={self.mxvariable:.2f},{self.myvariable:.2g}")
        slhaf = self.M.createSLHAFile( )
        sigmacut=.02*fb
        if max(self.M.masses)>1600:
            sigmacut=.01*fb
        if max(self.M.masses)>1800:
            sigmacut=.003*fb
        if max(self.M.masses)>2000:
            sigmacut=.001*fb
        ## first get rmax
        if hasattr ( self.predictor, "predictions" ):
            del self.predictor.predictions
        from builder.manipulator import Manipulator
        manipulator = Manipulator ( self.M, self.environ )
        worked = self.predictor.predict ( manipulator, keep_predictions = True )
        cr, _ = self.critic.predict_critic ( self.M, keep_predictions = True )
        ret = { "nll": None, "critic": None, "oul": None, "eul": None }

        self.M.delCurrentSLHA()
        critics={ "llhd": None, "ul": self.M.ul_critic }
        # max_allowed, n_excluding = critics["ul"]["max_allowed"], critics["ul"]["n_excluding"]
        # critics["ul"]["passes"] = (max_allowed+1) >= n_excluding ## loosened!
        if hasattr ( self.M, "llhd_critic" ):
            r, rexp = self.M.llhd_critic["robs"], self.M.llhd_critic["rexp"]
            # passes = (r < 1.2) or (r < 1.5 and r/rexp < .5 ) # thats what we really have
            passes = (r < 1.3) or (r < 1.6 and r/rexp < .55 ) ## loosened!!
            self.M.llhd_critic["passes"] = passes
            critics["llhd"] = self.M.llhd_critic
            ret["critic"] = critics

        if not worked:
            self.error( f"worked: {worked}" )
        if not worked:
            return ret
        ## now get the likelihoods
        nlls={}
        ## start with the SM likelihood
        nlls[0.] = self.getNLLs ( self.predictor.predictions, mu=0. )
        ## get for the others FIXME should adapt to ssm?
        for mu in numpy.arange(.4,1.8,.05):
            rounded_mu = float(round(mu,5))
            nlls[rounded_mu] = self.getNLLs ( self.predictor.predictions, mu=mu )
        ret["nll"] = nlls
        ouls = self.getLimits ( self.predictor.predictions, observed )
        ret["oul"] = ouls
        euls = self.getLimits ( self.predictor.predictions, apriori )
        ret["eul"] = euls

        del self.predictor.predictions

        return ret

    def getLimits ( self, predictions : List[TheoryPrediction],
                    evaluationType : NllEvalType ) -> Dict:
        """ get the limits for all predictions

        :param evaluationType: one of: observed, apriori, aposteriori
        """
        limits = {}
        for tp in predictions:
            txname = ','.join ( set( [ i.txName for i in tp.txnames ] ) )
            dId = tp.dataId()
            if dId == "(combined)":
                dId = "combined"
            name = f"{tp.analysisId()}:{dId}:{txname}"
            limits[ name ] = tp.getUpperLimitOnMu (
                    evaluationType = evaluationType )
        return limits

    def getNLLs ( self, predictions, mu = 1. ) -> Dict:
        """ return dictionary with the nlls per analysis """
        nlls = {}
        for tp in predictions:
            txname = ','.join ( set( [ i.txName for i in tp.txnames ] ) )
            dId = tp.dataId()
            if dId == "(combined)":
                dId = "(comb)"
            name = f"{tp.analysisId()}:{dId}:{txname}"
            nlls[ name ] = tp.nll ( mu )
        return nlls

    def clean ( self ):
        """ clean up after the run """
        cmd = f"rm {self.M.currentSLHA}"
        subprocess.getoutput ( cmd )

    def setSSM ( self, pids : tuple, ssm : float ):
        """ set the signal strength multiplier
        """
        self.M.ssmultipliers[pids]=ssm

    def setParameter ( self, pid : Union[int,tuple], mass : float ):
        """ set mass or ssm of <pid> to <mass> """
        if type(pid) == tuple:
            return self.setSSM ( pid, mass )
        partners = [ ( 1000023, 1000024 ) ]
        self.M.masses[pid]=mass
        for pair in partners:
            if not pid in pair:
                continue
            for p in pair:
                if p in self.M.masses and self.massesAreTied ( p, pid ):
                    self.M.masses[p]=mass

    def setSSMultiplier ( self, pids : tuple, ssm : float ):
        """ set the ssm multipliers for pids to ssm.
        for now, set also for all signs
        """
        if pids in self.M.ssmultipliers:
            self.M.ssmultipliers[pids]=ssm
        pids1 = ( -pids[0], pids[0] )
        if pids1 in self.M.ssmultipliers:
            self.M.ssmultipliers[pids1]=ssm

    def hasResultsForPoint ( self, m1 : float, m2 : float ) -> bool:
        """ return true if we have already run point (m1,m2) """
        dictfile = self.getDictFileName ( m1, m2 )
        hasResult = os.path.exists ( dictfile )
        # self.pprint ( f"do we have a result for {m1:.2f},{m2:.2f}? {hasResult}" )
        if hasResult:
            return True
        return False

    def run ( self, rxvariable, ryvariable ):
        """ run for the points given """
        oldmasses = {}
        masspoints=self.getAllMassPoints()
        nxvariables = len(rxvariable)
        ct = 0
        for i1,m1 in enumerate(rxvariable):
            thrnr = 0
            try:
                thrnr = int ( self.threadnr.replace("nll","") )
            except Exception as e:
                pass
            setnr = i1+1 + thrnr * ( nxvariables )
            self.pprint ( f"now starting with point set #{setnr} [of {nxvariables} in this thread]" )
            self.pprint ( f"this point set contains {len(ryvariable)} points" )
            self.setParameter ( self.xvariable, m1 )
            if type(self.mxvariable)==int:
                self.M.masses[self.xvariable]=self.mxvariable ## reset LSP mass
            if type(self.mxvariable)==tuple:
                ## reset LSP mass
                self.setSSMultiplier ( self.xvariable, self.mxvariable )
            if type(self.myvariable)==int:
                self.M.masses[self.yvariable]=self.myvariable ## reset LSP mass
            if type(self.myvariable)==tuple:
                ## reset LSP mass
                self.setSSMultiplier ( self.yvariable, self.myvariable )
            for k,v in oldmasses.items():
                self.pprint ( f"WARNING: setting mass of {namer.asciiName(k)} back to {v}" )
                self.M.masses[k]=v
            oldmasses={}
            self.M.delXSecs() ## make sure we compute
            xsecs = self.M.getXsecs()
            xsectot = 0.*fb
            if len(xsecs)>0:
                for xsec in xsecs[0]:
                    xsectot += xsec.value
            if xsectot.asNumber ( fb ) < 1e-10:
                self.pprint ( "WARNING no xsec??" )
            for i2,m2 in enumerate(ryvariable):
                if m2 > m1: ## we assume yvariable to be the daughter
                    continue
                if m2 < 0.:
                    self.warning ( f"m2({namer.asciiName(self.yvariable)})={m2:.1f}<0. skipping!" )
                    continue
                if self.hasResultsForPoint ( m1, m2 ):
                    continue
                # self.pprint ( f"processing m({m1:.2f},{m2:.2f})" )
                if type(self.yvariable)==int:
                    self.M.masses[self.yvariable]=m2
                if type(self.yvariable)==tuple:
                    self.setSSMultiplier ( self.yvariable, m2 )
                for pid_,m_ in self.M.masses.items():
                    if pid_ != self.yvariable and m_ < m2: ## make sure LSP remains the LSP
                        self.warning ( f"have to raise {namer.asciiName(pid_)} {m_} -> {m2+1.}, so X1Z stays the LSP" )
                        oldmasses[pid_]=m_
                        self.M.masses[pid_]=m2 + 1.
                point = self.getPredictions ( False )
                nlls = point["nll"]
                if not nlls: continue
                nnlls,nnonzeroes=0,0

                for mu,nll in nlls.items():
                    nnlls+=len(nll)

                self.pprint ( f"{i1}/{nxvariables}: m({namer.asciiName(self.xvariable)})={m1:.1f}, m2({namer.asciiName(self.yvariable)})={m2:.1f}, {len(nlls)} mu's, {nnlls} nlls." )
                point["mx"] = float ( m1 )
                point["my"] = float ( m2 )
                masspoints.append ( point )
                self.addNewPoint ( point ) ## add the point
        return masspoints

def runThread ( threadid: int, obj, rxvariable, ryvariable,
        return_dict : Union[Dict,None] = None ):
    """ the method needed for parallelization to work """

    thread = NLLThread ( f"nll{threadid}", obj )
    newpoints = thread.run ( rxvariable, ryvariable )
    if return_dict != None:
        return_dict[threadid]=newpoints
    thread.clean()
    # thread.updatePickleFile()
    return newpoints

class NLLScanner ( LoggerBase ):
    """ class that encapsulates a likelihood sweep """
    def __init__ ( self, protomodel, xvariable, yvariable, nproc,
                   environ : RunEnviron, skip_production : bool = False,
                   dry_run : bool = False, dict_file : bool = False,
                   output : str = "nll" ):
        """
        :param rundir: the rundir
        :param environ: the RunEnviron
        :param skip_production: if possible, skip production, go to plotting
        :param dry_run: dont actually perform the actions
        :param output: prefix for output file [nll]
        """
        super ( NLLScanner, self ).__init__ ( "nll" )
        self.dry_run = dry_run
        self.output = output
        self.dict_file = dict_file
        self.environ = environ
        self.M = protomodel
        self.xvariable = xvariable
        self.yvariable = yvariable
        x_short = namer.asciiName(self.xvariable).replace(',','').replace(' ','')
        y_short = namer.asciiName(self.yvariable).replace(',','').replace(' ','')
        picklefile = f"{self.output}{x_short}{y_short}.pcl"
        self.picklefile = picklefile
        self.nproc = nproc
        self.skip_production = skip_production
        self.predictor = Predictor ( 'nll', environ=self.environ )
        self.critic = Critic ( 'nll', environ=self.environ )
        self.cprint ( "yellow", f"starting with {nproc} threads" )
        self.pprint ( f"self.predictor = Predictor ( 'nll', environ='{self.environ.runDictFile}' )" )
        yname = moreHelpers.shortYVarName( self.yvariable )
        self.resultsdir = f"{self.environ.rundir}/nlls_{namer.asciiName(self.xvariable)}{yname}/"

    def describeRange ( self, r ):
        """ describe range r in a string """
        if len(r)==0:
            return ""
        if len(r)==1:
            return f"{r[0]:.2f}"
        if len(r)==2:
            return f"{r[0]:.2f},{r[1]:.2f}"
        return f"{r[0]:.2f},{r[1]:.2f} ... {r[-1]:.2f} -> {len(r)} points"

    def runForMassPoints ( self, rxvariable, ryvariable ):
        """ run for the given mass points
        :param rxvariable: list of masses for xvariable
        :param ryvariable: list of masses for yvariable
        :returns: masspoints
        """
        if self.dry_run:
            self.pprint ( f"dry_run. stopping here" )
            sys.exit()
            self.pprint ( f"dry_run. would run for xvariable={rxvariable}" )
            self.pprint ( f"yvariable={ryvariable}" )
            sys.exit()
        np.random.shuffle ( rxvariable )
        mask = []
        thread = NLLThread ( "nll0", self )
        for rxv in rxvariable:
            hasMissing = False
            for rxy in ryvariable:
                hasFile = thread.hasResultsForPoint ( rxv, rxy )
                if not hasFile:
                    hasMissing = True
            mask.append ( hasMissing )
        rxvariable = rxvariable[mask]
        if len(rxvariable)==0:
            thread.updatePickleFile()
            return

        if self.nproc == 1:
            return runThread ( 0, self, rxvariable, ryvariable )
        chunkedRxvariable = [ list(rxvariable[i::self.nproc]) for i in range(self.nproc) ]
        processes = []
        manager = multiprocessing.Manager()
        return_dict=manager.dict()
        for ctr,chunk in enumerate(chunkedRxvariable):
            self.M.walkerid = 2000+ctr
            p = multiprocessing.Process ( target = runThread, args = ( ctr, self, chunk, ryvariable, return_dict ) )
            p.start()
            processes.append ( p )

        for p in processes:
            p.join()
        thread.updatePickleFile()

    def scanLikelihoodFor ( self, range1 : Dict, range2 : Dict,
                            nevents : int, topo : str ):
        """ plot the likelihoods as a function of xvariable and yvariable

        :param range1: dictionary for range1 with min, max, dm
        :param range2: dictionary for range1 with min, max, dm
        """
        self.nevents = nevents
        self.topo = topo
        xvariable = self.xvariable
        yvariable = self.yvariable
        if yvariable != self.M.LSP:
            self.pprint ( f"we currently assume yvariable to be the mass of the LSP, but it is {yvariable}" )
        if os.path.exists ( self.picklefile ) and self.skip_production:
            self.pprint ( f"we were asked to skip production: {self.picklefile} exists." )
            return
        import numpy
        c = Combiner()
        anaIds = c.getAnaIdsWithPids ( self.M.bestCombo, [ xvariable, yvariable ] )
        ## mass range for xvariable
        if type(yvariable) == int:
            self.myvariable = self.M.masses[yvariable]
        if type(yvariable) == tuple:
            self.myvariable = self.M.ssmultipliers[yvariable]

        # choose the axis boundaries and step sizes such that the hiscore values
        # are nicely central
        from numpy import ceil
        ndxmin = int ( ceil (( self.mxvariable - range1["min"] ) / range1["dm"]) )
        ndxmax = int ( ceil (( range1["max"] - self.mxvariable ) / range1["dm"]) )
        rxvariable = numpy.arange ( self.mxvariable - ndxmin*range1["dm"],
                       self.mxvariable + ndxmax * range1["dm"] + 1e-5, range1["dm"] )
        # rxvariable = numpy.arange ( range1["min"], range1["max"]+1e-8, range1["dm"] )
        # rxvariable = numpy.insert ( rxvariable, 8, self.mxvariable )
        ndymin = int ( ceil (( self.myvariable - range2["min"] ) / range2["dm"]) )
        ndymay = int ( ceil (( range2["max"] - self.myvariable ) / range2["dm"]) )
        ryvariable = numpy.arange ( self.myvariable - ndymin*range2["dm"],
                       self.myvariable + ndymay * range2["dm"] + 1e-5, range2["dm"] )

        #ryvariable = numpy.arange ( range2["min"], range2["max"]+1e-8, range2["dm"] )
        #ryvariable = numpy.insert ( ryvariable, 8, self.myvariable )
        self.cprint ( "green", f"range for {namer.asciiName(xvariable)}: {self.describeRange( rxvariable )}" )
        self.cprint ( "green", f"range for {namer.asciiName(yvariable)}: {self.describeRange( ryvariable )}" )
        self.cprint ( "green", f"total {len(rxvariable)*len(ryvariable)} points, {nevents} events for {topo}" )
        self.M.createNewSLHAFileName ( prefix=f"nll{xvariable}" )
        #self.M.initializePredictor()
        self.predictor.filterForTopos ( topo )
        self.M.walkerid = 2000

        thread0 = NLLThread ( "nll0", self )
        thread0.ntotal = len(rxvariable)*len(ryvariable)+1
        thread0.writeRunMeta()
        if not thread0.hasResultsForPoint ( self.mxvariable, self.myvariable ):
            point = thread0.getPredictions ( False )
            point["mx"] = float ( self.mxvariable )
            point["my"] = float ( self.myvariable )
            thread0.addNewPoint ( point )
            nlls = point["nll"]
            critics = point["critic"]
            thread0.clean()
            self.pprint ( f"protomodel point: m1({namer.asciiName(self.xvariable)})={self.mxvariable:.2f}, m2({namer.asciiName(self.yvariable)})={self.myvariable:.2f}, {len(nlls)} nlls" )
            # point [ "mx" ] = float ( self.mxvariable )
            # point [ "my" ] = float ( self.myvariable )
            masspoints = [ point ]
        else:
            masspoints = thread0.getAllMassPoints()

        if False:
            ## freeze out all other particles? We shouldnt!
            for pid_,m_ in self.M.masses.items():
                if pid_ not in [ self.xvariable, self.yvariable ]:
                    self.M.masses[pid_]=1e6

        self.runForMassPoints ( rxvariable, ryvariable )
        self.M.delCurrentSLHA()

    def cleanFirst ( self ):
        """ clean results dir and pickle file before running """
        if os.path.exists ( self.picklefile ):
            self.pprint ( f"cleaning out {self.picklefile}" )
            os.unlink ( self.picklefile )
        if self.dict_file:
            dictfile = self.picklefile.replace(".pcl",".dict")
            if os.path.exists ( dictfile ):
                self.pprint ( f"cleaning out {dictfile}" )
                os.unlink ( dictfile )
        if os.path.exists ( self.resultsdir ):
            self.pprint ( f"cleaning out {self.resultsdir}" )
            shutil.rmtree ( self.resultsdir )

    def overrideWithDefaults ( self, args ):
        topo = { 1000005: "T2bb",1000006: "T2tt", 2000006: "T2tt", 1000021: "T1", \
                 1000023: "electroweakinos,stops,TChiZISRqq,TChiISR",
                 1000024: "electroweakinos,stops,TChiZISRqq,TChiISR",
                 (1000024,1000023): "electroweakinos,stops,TChiZISRqq,TChiISR",
                 1000025: "electroweakinos,stops,TChiZISRqq,TChiISR",
                 1000035: "electroweakinos,stops,TChiZISRqq,TChiISR",
                 1000001: "T2",  1000002: "T2", 1000003: "T2", 1000004: "T2" }
        ### make the LSP scan depend on the mother
        if args.topo == None:
            args.topo = topo[args.xvariable]
        # self.mxvariable = self.M.masses[self.xvariable]
        if type(self.xvariable) == int:
            self.mxvariable = self.M.masses[self.xvariable]
        if type(self.xvariable) == tuple:
            if self.xvariable[0]> self.xvariable[1]:
                self.xvariable = ( self.xvariable[1], self.xvariable[0] )
            self.mxvariable = self.M.ssmultipliers[self.xvariable]
        if type(self.yvariable) == int:
            self.myvariable = self.M.masses[self.yvariable]
        if type(self.yvariable) == tuple:
            if self.yvariable[0]> self.yvariable[1]:
                self.yvariable = ( self.yvariable[1], self.yvariable[0] )
            self.myvariable = self.M.ssmultipliers[self.yvariable]
        nbinsx, nbinsy = 20, 20 # how many bins do we want per dimension
        if args.minx == None:
            args.minx = self.mxvariable*.6
        if args.maxx == None:
            args.maxx = self.mxvariable*1.7
        if args.deltamx == None:
            args.deltamx = ( args.maxx - args.minx ) / nbinsx
        if args.miny == None:
            if type(self.yvariable) == int:
                args.miny = max ( self.myvariable*.6 - 10., 1. )
            if type(self.yvariable) == tuple:
                # args.miny = self.myvariable*.2
                args.miny = 0.
        if args.maxy == None:
            if type(self.yvariable) == int:
                args.maxy = self.myvariable*1.9 + 10.
            if type(self.yvariable) == tuple:
                args.maxy = self.myvariable*5.
        if args.deltamy == None:
            args.deltamy = ( args.maxy - args.miny ) / nbinsy
        return args

def main ():
    import argparse
    argparser = argparse.ArgumentParser(
            description='perform likelhood scans')
    argparser.add_argument ( '-n', '--number',
            help='which hiscore to plot [0]',
            type=int, default=0 )
    argparser.add_argument ( '-x', '--xvariable',
            help='variable for x axis, e.g. 1000006 or "Xt" [Xt]',
            type=str, default='Xt' )
    argparser.add_argument ( '-y', '--yvariable',
            help='variable for y axis, e.g. 1000022 or "(Xt,Xt)", [X1Z]',
            type=str, default="X1Z" )
    argparser.add_argument ( '-P', '--nproc',
            help='number of processes to run in parallel. zero is autodetect.'\
                 'Negative numbers are added to autodetect [0]',
            type=int, default=0 )
    argparser.add_argument ( '-mx', '--minx',
            help='minimum mass of xvariable [None]',
            type=float, default=None )
    argparser.add_argument ( '-Mx', '--maxx',
            help='maximum mass of xvariable [None]',
            type=float, default=None )
    argparser.add_argument ( '-dx', '--deltamx',
            help='delta m of xvariable [None]',
            type=float, default=None )
    argparser.add_argument ( '-my', '--miny',
            help='minimum mass of yvariable [None]',
            type=float, default=None )
    argparser.add_argument ( '-My', '--maxy',
            help='maximum mass of yvariable [None]',
            type=float, default=None )
    argparser.add_argument ( '-dy', '--deltamy',
            help='delta m of yvariable [None]',
            type=float, default=None )
    argparser.add_argument ( '-t', '--topo',
            help='topology [None]',
            type=str, default=None )
    argparser.add_argument ( '-e', '--nevents',
            help='number of events [50000]',
            type=int, default=50000 )
    argparser.add_argument ( '-H', '--hiscores',
            help='hiscore file to draw from [<rundir>/hiscores_global.dict]',
            type=str, default="default" )
    argparser.add_argument ( '-D', '--draw',
            help='also perform the plotting, ie call plotnlls',
            action='store_true' )
    argparser.add_argument ( '-K', '--dontkeep',
            help='remove resultsdir after finished',
            action='store_true' )
    argparser.add_argument ( '-c', '--clean_first',
            help='clean pickle files and results dir before running',
            action='store_true' )
    argparser.add_argument ( '--dry_run',
            help='just tell us what you would be doing, dont actually do it',
            action='store_true' )
    argparser.add_argument ( '-v', '--verbosity',
            help='verbosity -- debug, info, warn, err [info]',
            type=str, default="info" )
    argparser.add_argument ( '-o', '--output',
            help="prefix for output file [nll]",
            type=str, default="nll" )
    argparser.add_argument ( '-u', '--uploadTo',
            help="where do we upload to, on smodels.github.io [latest]",
            type=str, default="latest" )
    argparser.add_argument ( '-r', '--run_environ',
            help="path to run environment [./run.dict]",
            type=str, default="./run.dict" )
    argparser.add_argument ( '-S', '--skip_production',
            help='if possible, skip production', action='store_true' )
    argparser.add_argument ( '--dict_file',
            help='write out as dict file as well', action='store_true' )
    args = argparser.parse_args()
    # rundir = setup( args.rundir )
    nproc = args.nproc
    environ = RunEnviron( args.run_environ )
    if nproc < 1:
        nproc = nCPUs() + nproc
    if args.hiscores == "default":
        args.hiscores = f"{environ.rundir}/hiscores_global.dict"
    from ptools.hiscoreTools import fetchHiscoresObj
    hi = fetchHiscoresObj ( args.hiscores, None, environ, walkerid="nll" )
    protomodel = hi.hiscores[0]
    #self.pprint ( f"fetched {protomodel} from {args.hiscores}" )

    xvariables = [ namer.pid ( args.xvariable ) ]
    if args.xvariable == 0:
        xvariables = findPids( environ.rundir )
    for xvariable in xvariables:
        yvariable = namer.pid ( args.yvariable )
        scanner = NLLScanner( protomodel, xvariable, yvariable, nproc,
                environ = environ, skip_production = args.skip_production,
                dry_run = args.dry_run, dict_file = args.dict_file,
                output = args.output )
        args.xvariable = xvariable
        args = scanner.overrideWithDefaults ( args )
        if args.clean_first:
            scanner.cleanFirst()
        range1 = { "min": args.minx, "max": args.maxx, "dm": args.deltamx }
        range2 = { "min": args.miny, "max": args.maxy, "dm": args.deltamy }
        scanner.scanLikelihoodFor ( range1, range2, args.nevents, args.topo )
        if args.dontkeep:
            scanner.unlinkResultsDir()
        if args.draw:
            verbose = args.verbosity
            copy = True
            max_anas = 5
            interactive = False
            drawtimestamp = True
            compress = False
            upload = args.uploadTo
            from plotting import plotNLLs
            plot = plotNLLs.NLLPlot ( xvariable, yvariable, verbose, copy,
                       max_anas, interactive, drawtimestamp, compress, environ,
                       upload )
            plot.writeScriptFile ( )
            plot.plot()

if __name__ == "__main__":
    main()
