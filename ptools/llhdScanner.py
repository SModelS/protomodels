#!/usr/bin/env python3

""" script used to produce the likelihood scans """

__all__ = [ "LlhdScanner" ]

import os, sys, multiprocessing, time, numpy, subprocess, copy, glob
import pickle, random, shutil
try:
    from protomodels.csetup import setup
    setup()
except ModuleNotFoundError as e:
    pass
from smodels.tools.wrapperBase import WrapperBase
WrapperBase.defaulttempdir="./" ## keep the temps in our folder
from smodels.base.physicsUnits import fb
from smodels.base.runtime import nCPUs
from smodels.matching.theoryPrediction import TheoryPrediction
from tester.combiner import Combiner
from tester.predictor import Predictor
from tester.critic import Critic
from plotting import plotLlhds
from typing import Dict, Tuple, Union, List
from ptools.sparticleNames import SParticleNames
from ptools import moreHelpers
from base.loggerbase import LoggerBase

namer = SParticleNames ( False )

def findPids ( rundir ):
    """ search for llhd*pcl files, report the corresponding pids.
    :returns: set of pids
    """
    ret = set()
    files = glob.glob("llhd*pcl")
    files += glob.glob( f"{rundir}/llhd*pcl" )
    for f in files:
        p = f.find("llhd")
        s = f[p+4:]
        s = s.replace(".pcl","")
        s = s.replace("1000022","")
        s = s.replace("X1Z","")
        ret.add ( int(s) )
    print ( f"[llhdScanner] pids are {ret}" )
    return ret

class LlhdThread ( LoggerBase ):
    """ one thread of the sweep """
    def __init__ ( self, threadnr: int, obj ):
        """ the constructor. 
        """
        super ( LlhdThread, self ).__init__ ( threadnr )
        self.rundir = setup( obj.rundir )
        yname = moreHelpers.shortYVarName( obj.yvariable )
        self.resultsdir = f"{self.rundir}/llhds_{namer.asciiName(xvariable)}{yname}/"
        self.topo = obj.topo
        self.threadnr = obj.threadnr
        self.picklefile = obj.picklefile
        self.M = copy.deepcopy ( obj.protomodel )
        self.origmasses = copy.deepcopy ( self.M.masses )
        self.origssmultipliers = copy.deepcopy ( self.M.ssmultipliers )
        self.M.createNewSLHAFileName ( prefix=f"lthrd{threadnr}_{xvariable}" )
        self.xvariable = obj.xvariable
        self.yvariable = obj.yvariable
        self.mxvariable = obj.mxvariable
        self.myvariable = obj.myvariable
        self.nevents = obj.nevents
        self.predictor = obj.predictor
        self.critic = obj.critic
        self.mkResultsDir()

    def getDefaultDictionary ( self ):
        """ initialise the dictionary for the pickle file """
        d = { "masspoints": [], "mxvariable": self.mxvariable,
              "myvariable": self.myvariable, "nevents": self.nevents,
              "topo": self.topo, "timestamp": time.asctime(),
              "xvariable": self.xvariable, "yvariable": self.yvariable,
              "model": self.M.dict() }
        return d

    def writePickleFile ( self, d : Dict ):
        """ write the dictionary into the picklefile """
        #if os.path.exists ( self.picklefile ) and \
        #        os.stat ( self.picklefile ).st_size > 1000:
        #    try:
        #        f = open ( self.picklefile, "rb" )
        #        d = pickle.load(f) ## ok can read. copy, then!
        #        subprocess.getoutput ( f"cp {self.picklefile} {self.picklefile}.old" )
        #    except Exception as e:
        #        pass
        f = open ( self.picklefile, "wb" )
        pickle.dump ( d, f )
        f.close()

    def lockPickleFile ( self ):
        """ make sure we write sequentially """
        lockfile = self.picklefile+".lock"
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
        lockfile = self.picklefile+".lock"
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

    def mkResultsDir ( self ):
        """ make a results dir if it doesnt exist """
        self.pprint ( f"results will go to {self.resultsdir.replace(self.rundir,'.')}" ) 
        if os.path.exists ( self.resultsdir ):
            return
        os.mkdir ( self.resultsdir )


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
            f.write ( f"{point}\n" )
            f.close()
        nfiles = len ( glob.glob ( f"{self.resultsdir}/*.dict" ) )
        if nfiles % 100 == 0: # update with every 20th entry
            self.updatePickleFile()

    def getAllMassPoints ( self ):
        """ retrieve all mass points from resultsdir """
        files = glob.glob ( f"{self.resultsdir}/*.dict" )
        masspoints = []
        for fname in files:
            with open ( fname, "rt" ) as h:
                d = eval(h.read())
                masspoints.append ( d )
        return masspoints

    def updatePickleFile ( self ):
        """ collect all the entries in resultsdir, and compile them
        into one big pickle file """
        self.pprint ( f"update {self.picklefile}" )
        self.lockPickleFile()
        Dict = self.getDefaultDictionary()
        files = glob.glob ( f"{self.resultsdir}/*.dict" )
        masspoints = self.getAllMassPoints()
        Dict["masspoints"] = masspoints
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
        :returns: a diction with likelihoods ("llhd"), critics' responses ("critic"),
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
        worked = self.predictor.predict ( self.M, keep_predictions = True )
        
        ## now get the likelihoods
        llhds={}
        ## start with the SM likelihood
        llhds[0.] = self.getLikelihoods ( self.predictor.predictions, mu=0. )
        ## get for the others FIXME should adapt to ssm?
        for mu in numpy.arange(.4,1.8,.05):
            llhds[float(mu)] = self.getLikelihoods ( self.predictor.predictions, mu=mu )
        ouls = self.getLimits ( self.predictor.predictions, False )
        euls = self.getLimits ( self.predictor.predictions, True )
        del self.predictor.predictions
        self.M.delCurrentSLHA()
        critics={}
        for critic in self.M.critic_description.split(","):
            tokens = critic.split(":")
            if len(tokens)>1:
                critics[tokens[0]]=float(tokens[1])

        return { "llhd": llhds, "critic": critics, "oul": ouls, "eul": euls }

    def getLimits ( self, predictions : List[TheoryPrediction], 
                    expected : bool ) -> Dict:
        """ get the limits for all predictions 

        :param expected: if true, get expected limits on mu, else observed
        """
        limits = {}
        for tp in predictions:
            txname = ','.join ( set( [ i.txName for i in tp.txnames ] ) )
            dId = tp.dataId()
            if dId == "(combined)":
                dId = "(comb)"
            name = f"{tp.analysisId()}:{dId}:{txname}"
            limits[ name ] = tp.getUpperLimitOnMu ( expected = expected )
        return limits

    def getLikelihoods ( self, predictions, mu = 1. ) -> Dict:
        """ return dictionary with the likelihoods per analysis """
        llhds= {}
        for tp in predictions:
            txname = ','.join ( set( [ i.txName for i in tp.txnames ] ) )
            dId = tp.dataId()
            if dId == "(combined)":
                dId = "(comb)"
            name = f"{tp.analysisId()}:{dId}:{txname}"
            llhds[ name ] = tp.likelihood ( mu )
        return llhds

    def clean ( self ):
        """ clean up after the run """
        cmd = f"rm {self.M.currentSLHA}"
        subprocess.getoutput ( cmd )

    def setMass ( self, pid : int, mass : float ):
        """ set mass of <pid> to <mass> """
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
            self.pprint ( f"now starting with {i1}/{nxvariables}" )
            self.setMass ( self.xvariable, m1 )
            if type(self.myvariable)==int:
                self.M.masses[self.yvariable]=self.myvariable ## reset LSP mass
            if type(self.myvariable)==tuple:
                ## reset LSP mass
                self.setSSMultiplier ( self.yvariable, self.myvariable )
            for k,v in oldmasses.items():
                self.pprint ( f"WARNING: setting mass of {k} back to {v}" )
                self.M.masses[k]=v
            oldmasses={}
            self.M.delXSecs() ## make sure we compute
            xsecs = self.M.getXsecs()
            xsectot = 0.*fb
            for xsec in xsecs[0]:
                xsectot += xsec.value
            if xsectot.asNumber ( fb ) < 1e-10:
                self.pprint ( "WARNING no xsec??" )
            for i2,m2 in enumerate(ryvariable):
                if m2 > m1: ## we assume yvariable to be the daughter
                    continue
                if self.hasResultsForPoint ( m1, m2 ):
                    continue
                self.pprint ( f"processing m({m1:.2f},{m2:.2f})" )
                if type(self.yvariable)==int:
                    self.M.masses[self.yvariable]=m2
                if type(self.yvariable)==tuple:
                    self.setSSMultiplier ( self.yvariable, m2 )
                for pid_,m_ in self.M.masses.items():
                    if pid_ != self.yvariable and m_ < m2: ## make sure LSP remains the LSP
                        self.pprint ( f"WARNING: have to raise {pid_} from {m_} to {m2+1.}" )
                        oldmasses[pid_]=m_
                        self.M.masses[pid_]=m2 + 1.
                point = self.getPredictions ( False )
                llhds = point["llhd"]
                nllhds,nnonzeroes=0,0
                for mu,llhd in llhds.items():
                    nllhds+=len(llhd)

                self.pprint ( f"{i1}/{nxvariables}: m({namer.asciiName(self.xvariable)})={m1:.1f}, m2({namer.asciiName(self.yvariable)})={m2:.1g}, {len(llhds)} mu's, {nllhds} llhds." )
                point["mx"] = m1 
                point["my"] = m2
                masspoints.append ( point )
                self.addNewPoint ( point ) ## add the point
        return masspoints

def runThread ( threadid: int, obj, rxvariable, ryvariable, 
        return_dict : Union[Dict,None] = None ):
    """ the method needed for parallelization to work """

    thread = LlhdThread ( threadid, obj )
    newpoints = thread.run ( rxvariable, ryvariable )
    if return_dict != None:
        return_dict[threadid]=newpoints
    thread.clean()
    # thread.updatePickleFile()
    return newpoints

class LlhdScanner ( LoggerBase ):
    """ class that encapsulates a likelihood sweep """
    def __init__ ( self, protomodel, xvariable, yvariable, nproc, rundir : str,
                   dbpath : str = "official", select : str = "all",
                   do_srcombine : bool = True, skip_production : bool = False, 
                   dry_run : bool = False ):
        """
        :param rundir: the rundir
        :param dbpath: the database path
        :param skip_production: if possible, skip production, go to plotting
        :param dry_run: dont actually perform the actions
        """
        super ( LlhdScanner, self ).__init__ ( 0 )
        self.dry_run = dry_run
        self.rundirarg = rundir
        self.rundir = setup( rundir )
        self.M = protomodel
        self.xvariable = xvariable
        self.yvariable = yvariable
        self.nproc = nproc
        self.skip_production = skip_production
        self.predictor = Predictor ( 0, dbpath=dbpath, 
                select=select, do_srcombine = do_srcombine )
        self.critic = Critic ( 0, dbpath=dbpath, 
                select=select, do_srcombine = do_srcombine )
        self.pprint ( f"self.predictor = Predictor ( 0, dbpath='{dbpath}', select='{select}', do_srcombine = {do_srcombine} )" )

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
        random.shuffle ( rxvariable )
        mask = []
        thread = LlhdThread ( 0, self.rundir, self.M, self.xvariable, 
                self.yvariable, self.mxvariable, self.myvariable, self.nevents,
                self.predictor, self.picklefile, self.topo )
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
                            nevents : int, topo : str, output : str ):
        """ plot the likelihoods as a function of xvariable and yvariable

        :param range1: dictionary for range1 with min, max, dm
        :param range2: dictionary for range1 with min, max, dm
        :param output: prefix for output file [mp]
        """
        self.nevents = nevents
        self.topo = topo
        xvariable = self.xvariable
        yvariable = self.yvariable
        if yvariable != self.M.LSP:
            self.pprint ( f"we currently assume yvariable to be the mass of the LSP, but it is {yvariable}" )
        picklefile = f"{output}{namer.asciiName(xvariable)}{namer.asciiName(yvariable).replace(',','').replace(' ','')}.pcl"
        self.picklefile = picklefile
        if os.path.exists ( picklefile ) and self.skip_production:
            self.pprint ( f"we were asked to skip production: {picklefile} exists." )
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
        print ( f"[llhdScanner] range for {namer.asciiName(xvariable)}: {self.describeRange( rxvariable )}" )
        print ( f"[llhdScanner] range for {namer.asciiName(yvariable)}: {self.describeRange( ryvariable )}" )
        print ( f"[llhdScanner] total {len(rxvariable)*len(ryvariable)} points, {nevents} events for {topo}" )
        self.M.createNewSLHAFileName ( prefix=f"llhd{xvariable}" )
        #self.M.initializePredictor()
        self.predictor.filterForTopos ( topo )
        self.M.walkerid = 2000

        thread0 = LlhdThread ( 0, self )
        thread0.ntotal = len(rxvariable)*len(ryvariable)+1
        thread0.writeRunMeta()
        if not thread0.hasResultsForPoint ( self.mxvariable, self.myvariable ):
            point = thread0.getPredictions ( False )
            point["mx"] = self.mxvariable
            point["my"] = self.myvariable
            thread0.addNewPoint ( point )
            llhds = point["llhd"]
            critics = point["critic"]
            thread0.clean()
            self.pprint ( f"protomodel point: m1({namer.asciiName(self.xvariable)})={self.mxvariable:.2f}, m2({namer.asciiName(self.yvariable)})={self.myvariable:.2f}, {len(llhds)} llhds" )
            point [ "mx" ] = self.mxvariable
            point [ "my" ] = self.myvariable
            masspoints = [ point ]
        else:
            masspoints = thread0.getAllMassPoints()

        if False:
            ## freeze out all other particles? We shouldnt!
            for pid_,m_ in self.M.masses.items():
                # print ( f"@@a freezing {pid_}? {pid_ not in [ self.xvariable, self.yvariable ]}" )
                if pid_ not in [ self.xvariable, self.yvariable ]:
                    self.M.masses[pid_]=1e6

        self.runForMassPoints ( rxvariable, ryvariable )
        self.M.delCurrentSLHA()

    def overrideWithDefaults ( self, args ):
        topo = { 1000005: "T2bb",1000006: "T2tt", 2000006: "T2tt", 1000021: "T1", \
                 1000023: "electroweakinos,stops", 
                 1000024: "electroweakinos,stops",
                 1000001: "T2",  1000002: "T2", 1000003: "T2", 1000004: "T2" }
        ### make the LSP scan depend on the mother
        if args.topo == None:
            args.topo = topo[args.xvariable]
        self.mxvariable = self.M.masses[self.xvariable]
        if type(self.yvariable) == int:
            self.myvariable = self.M.masses[self.yvariable]
        if type(self.yvariable) == tuple:
            self.myvariable = self.M.ssmultipliers[self.yvariable]
        nbinsx, nbinsy = 20, 20 # how many bins do we want per dimension
        if args.min1 == None:
            args.min1 = self.mxvariable*.6
        if args.max1 == None:
            args.max1 = self.mxvariable*1.7
        if args.deltam1 == None:
            args.deltam1 = ( args.max1 - args.min1 ) / nbinsx
        if args.min2 == None:
            if type(self.yvariable) == int:
                args.min2 = max ( self.myvariable*.6 - 10., 1. )
            if type(self.yvariable) == tuple:
                # args.min2 = self.myvariable*.2
                args.min2 = 0.
        if args.max2 == None:
            if type(self.yvariable) == int:
                args.max2 = self.myvariable*1.9 + 10.
            if type(self.yvariable) == tuple:
                args.max2 = self.myvariable*5.
        if args.deltam2 == None:
            args.deltam2 = ( args.max2 - args.min2 ) / nbinsy
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
            help='number of process to run in parallel. zero is autodetect. Negative numbers are added to autodetect [0]',
            type=int, default=0 )
    argparser.add_argument ( '-m1', '--min1',
            help='minimum mass of xvariable [None]',
            type=float, default=None )
    argparser.add_argument ( '-M1', '--max1',
            help='maximum mass of xvariable [None]',
            type=float, default=None )
    argparser.add_argument ( '-d1', '--deltam1',
            help='delta m of xvariable [None]',
            type=float, default=None )
    argparser.add_argument ( '-m2', '--min2',
            help='minimum mass of yvariable [None]',
            type=float, default=None )
    argparser.add_argument ( '-M2', '--max2',
            help='maximum mass of yvariable [None]',
            type=float, default=None )
    argparser.add_argument ( '-d2', '--deltam2',
            help='delta m of yvariable [None]',
            type=float, default=None )
    argparser.add_argument ( '-t', '--topo',
            help='topology [None]',
            type=str, default=None )
    argparser.add_argument ( '-R', '--rundir',
            help='override the default rundir [None]',
            type=str, default=None )
    argparser.add_argument ( '-e', '--nevents',
            help='number of events [50000]',
            type=int, default=50000 )
    argparser.add_argument ( '-H', '--hiscores',
            help='hiscore file to draw from [<rundir>/hiscores.dict]',
            type=str, default="default" )
    argparser.add_argument ( '-D', '--draw',
            help='also perform the plotting, ie call plotLlhds',
            action='store_true' )
    argparser.add_argument ( '-K', '--dontkeep',
            help='remove resultsdir after finished',
            action='store_true' )
    argparser.add_argument ( '--dry_run',
            help='just tell us what you would be doing, dont actually do it',
            action='store_true' )
    argparser.add_argument ( '-v', '--verbosity',
            help='verbosity -- debug, info, warn, err [info]',
            type=str, default="info" )
    argparser.add_argument ( '-o', '--output',
            help="prefix for output file [llhd]",
            type=str, default="llhd" )
    argparser.add_argument ( '-u', '--uploadTo',
            help="where do we upload to, on smodels.github.io [latest]",
            type=str, default="latest" )
    argparser.add_argument ( '-s', '--select',
            help="what do we select for [all]",
            type=str, default="all" )
    argparser.add_argument ( '--dbpath',
            help="path to database [official]",
            type=str, default="official" )
    argparser.add_argument ( '-c', '--do_srcombine',
            help='do_srcombine', action='store_true' )
    argparser.add_argument ( '-S', '--skip_production',
            help='if possible, skip production', action='store_true' )
    args = argparser.parse_args()
    rundir = setup( args.rundir )
    nproc = args.nproc
    if nproc < 1:
        nproc = nCPUs() + nproc
    if args.hiscores == "default":
        args.hiscores = f"{rundir}/hiscores.dict"
    from ptools.hiscoreTools import fetchHiscoresObj
    hi = fetchHiscoresObj ( args.hiscores, None, args.dbpath )
    protomodel = hi.hiscores[0]
    #self.pprint ( f"fetched {protomodel} from {args.hiscores}" )

    xvariables = [ namer.pid ( args.xvariable ) ]
    if args.xvariable == 0:
        xvariables = findPids( rundir )
    for xvariable in xvariables:
        yvariable = namer.pid ( args.yvariable )
        scanner = LlhdScanner( protomodel, xvariable, yvariable, nproc, rundir, 
                dbpath = args.dbpath, select = args.select, 
                do_srcombine = args.do_srcombine, 
                skip_production = args.skip_production,
                dry_run = args.dry_run )
        args.xvariable = xvariable
        args = scanner.overrideWithDefaults ( args )
        range1 = { "min": args.min1, "max": args.max1, "dm": args.deltam1 }
        range2 = { "min": args.min2, "max": args.max2, "dm": args.deltam2 }
        scanner.scanLikelihoodFor ( range1, range2, args.nevents, args.topo, 
                args.output )
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
            plot = plotLlhds.LlhdPlot ( xvariable, yvariable, verbose, copy, max_anas,
                  interactive, drawtimestamp, compress, rundir, upload, args.dbpath )
            plot.writeScriptFile ( )
            plot.plot()

if __name__ == "__main__":
    main()
