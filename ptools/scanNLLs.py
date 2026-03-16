#!/usr/bin/env python3

""" script used to produce the likelihood scans """

__all__ = [ "NLLScanner" ]

import os, sys, multiprocessing, time, numpy, subprocess, copy, glob
import pickle, random, shutil
import numpy as np
from typing import Dict, Tuple, Union, List
from base.pbase import prettyFileName

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
from smodels_utils.helper.terminalcolors import RED, GREEN, YELLOW, RESET, CYAN

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
        super ( NLLThread, self ).__init__ ( threadnr, "info" )
        self.obj = obj
        #self.predictor = copy.deepcopy ( obj.predictor )
        # self.critic = copy.deepcopy ( obj.critic )
        self.predictor = obj.predictor
        self.critic = obj.critic
        self.resultsdir = obj.resultsdir
        self.threadnr = threadnr
        self.M = copy.deepcopy ( obj.M )
        self.xvalue = obj.xvalue
        self.yvalue = obj.yvalue
        self.setSLHAFileName( self.xvalue, self.yvalue )
        helpers.mkdir ( self.resultsdir )

    def setSLHAFileName ( self, m1 : float, m2 : float ):
        """ set the slha file name
        :param m1: value of x coordinate
        :param m2: value of y coordinate
        :returns: filename
        """
        # self.pprint ( f"slhadir {self.obj.slhadir}" )
        fname = f"lthrd{self.threadnr}_{self.getBaseName(m1,m2)}.slha"
        self.M.createNewSLHAFileName ( prefix=fname )
        return self.M.currentSLHA

    def getBaseName ( self, m1 : float, m2 : float ):
        """ get our base name for a parameter point
        :param m1: value of x coordinate
        :param m2: value of y coordinate
        :returns: base name
        """
        x_name = namer.asciiName(self.obj.xvariable)
        y_name = namer.asciiName(self.obj.yvariable)
        #x_value = f"{self.xvalue:.3f}"
        #y_value = f"{self.yvalue:.3f}"
        x_value = f"{m1:.3f}"
        y_value = f"{m2:.3f}"
        basename = f"{x_name}_{x_value}-{y_name}_{y_value}"
        basename = basename.replace(" ","").replace("(","_").replace(")","")
        basename = basename.replace(",","")
        return basename

    def getDefaultDictionary ( self ):
        """ initialise the dictionary for the pickle file """
        import time
        d = { "parameterpoints": [], "xvalue": self.xvalue,
              "yvalue": self.yvalue, "nevents": self.obj.nevents,
              "topo": self.obj.topo, "timestamp": time.asctime(),
              "xvariable": self.obj.xvariable, "yvariable": self.obj.yvariable,
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
        meta["y_is_dm"] = self.obj.y_is_dm
        meta["profile_mu"] = self.obj.args["profile_mu"]
        meta["xvariable"]=self.obj.xvariable
        meta["yvariable"]=self.obj.yvariable
        return meta

    def createPickleBackup ( self ):
        if os.path.exists ( self.obj.picklefile ) and \
                os.stat ( self.obj.picklefile ).st_size > 1000:
            try:
                f = open ( self.obj.picklefile, "rb" )
                d = pickle.load(f) ## ok can read. copy, then!
                cmd = f"cp {self.obj.picklefile} {self.obj.picklefile}.old"
                subprocess.getoutput ( cmd )
            except Exception as e:
                pass

    def writeFiles ( self, d : Dict ):
        """ write the dictionary into the picklefile and dict file """
        # self.createPickleBackup()
        ## lets see if that removes the extra entries
        for pm in d["parameterpoints"]:
            for var in [ "xvariable", "yvariable" ]:
                if var in pm:
                    pm.pop ( var )
        f = open ( self.obj.picklefile, "wb" )
        pickle.dump ( d, f )
        f.close()
        writeDictFile = False
        if self.obj.args["dict_file"]:
            dictfile = self.obj.picklefile.replace(".pcl",".dict")
            self.pprint ( f"writing to {CYAN}{dictfile}{RESET}" )
            with open ( dictfile, "wt" ) as f:
                befores = { "x": "xvariable", "xvariable": "y",
                    "y": "yvariable", "yvariable": "critic", "critic": "oul",
                    "oul": "eul", "eul": "nll" }
                d = py_dumps ( d, level = 0, a_before = befores )
                # not needed to write comments, there is a meta
                f.write ( d )
            f.close()

    def lockPickleFile ( self ):
        """ make sure we write sequentially """
        lockfile = f"{self.obj.picklefile}.lock"
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
        lockfile = f"{self.obj.picklefile}.lock"
        if os.path.exists ( lockfile ):
            try:
                os.unlink ( lockfile )
            except FileNotFoundError as e:
                pass

    def isSameMassPoint ( self, point1 : Dict, point2 : Dict ) -> bool:
        """ are the two points identical? """
        dx = point1["x"]-point2["x"]
        dy = point1["y"]-point2["y"]
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
        if not "x" in point: #dont add anything!
            return
        point["x"]=round( point["x"], 7 )
        point["y"]=round( point["y"], 7 )
        # this information, xvariable, yvariable is not needed,
        # but helpful for debugging
        point["xvariable"]=self.obj.xvariable
        point["yvariable"]=self.obj.yvariable
        dictfile = self.getDictFileName ( point["x"], point["y"] )
        with open ( dictfile, "wt" ) as f:
            # f.write ( f"{point}\n" )
            befores = { "x": "xvariable", "xvariable": "y", "y": "yvariable",
                        "yvariable": "critic", "critic": "oul",
                        "oul": "eul", "eul": "nll" }
            d = py_dumps ( point, level = 0, a_before = befores )
            f.write ( d + "\n" )
            f.close()
        nfiles = len ( glob.glob ( f"{self.resultsdir}/*.dict" ) )
        if nfiles % 100 == 0: # update with every 20th entry
            self.updatePickleFile()

    def getAllParameterPoints ( self ) -> list:
        """ retrieve all mass points from resultsdir """
        files = glob.glob ( f"{self.resultsdir}/*.dict" )
        parameterpoints = []
        for fname in files:
            with open ( fname, "rt" ) as h:
                try:
                    d = eval(h.read())
                except (SyntaxError,ValueError,TypeError) as e:
                    self.error ( f"could not read: {fname}: {e}" )
                    sys.exit()
                for var in [ "x", "y", "mx", "my" ]:
                    if var in d:
                        d[var]=float(d[var])
                if not "x" in d and "mx" in d:
                    d["x"]=d["mx"]
                    d.pop("mx")
                if not "y" in d and "my" in d:
                    d["y"]=d["my"]
                    d.pop("my")
                parameterpoints.append ( d )
        return parameterpoints

    def updatePickleFile ( self ):
        """ collect all the entries in resultsdir, and compile them
        into one big pickle file """
        self.debug ( f"updating {self.obj.picklefile}" )
        self.lockPickleFile()
        Dict = self.getDefaultDictionary()
        files = glob.glob ( f"{self.resultsdir}/*.dict" )
        parameterpoints = self.getAllParameterPoints()
        Dict["parameterpoints"] = parameterpoints
        Dict["meta"] = self.getMetaInformation()
        self.writeFiles ( Dict )
        self.unlockPickleFile()

    def unlinkResultsDir ( self ):
        """ clean up in the end """
        if os.path.exists ( self.resultsdir ):
            shutil.rmtree ( self.resultsdir )

    def massesAreTied ( self, xvariable, yvariable ):
        """ are the masses of xvariable and yvariable tied originally? """
        if not xvariable in self.obj.origmasses:
            return False
        if not yvariable in self.obj.origmasses:
            return False
        dm = self.obj.origmasses[xvariable] - self.obj.origmasses[yvariable]
        if abs(dm)<1e-5:
            return True
        return False

    def getPredictions ( self, recycle_xsecs : bool,
           m1 : float, m2 : float ) -> Dict:
        """ get predictions, return likelihoods

        :param recycle_xsecs: if true, then recycle the cross sections, dont
        recompute
        :returns: a diction with likelihoods ("nll"), critics' responses ("critic"),
        observed ("oul") and expected ("eul") upper limits on mu.
        """
        self.debug ( f"asking for predictions for x,y={self.xvalue:.2f},{self.yvalue:.2g}")
        slhaf = self.M.createSLHAFile( )
        ## first get rmax
        if os.path.exists ( slhaf ) and self.obj.slhadir is not None:
            newf = f"{self.obj.slhadir}/{self.getBaseName(m1,m2)}.slha"
            shutil.copy ( slhaf, newf )
            snewf = newf.replace ( os.getcwd(), "." )
            self.debug ( f"created {snewf}" )
        sigmacut=0.*fb
        if hasattr ( self.predictor, "predictions" ):
            del self.predictor.predictions
        from builder.manipulator import Manipulator
        manipulator = Manipulator ( self.M, self.obj.environ )
        force_K = self.obj.args["profile_mu"]
        worked, explanation = self.predictor.predict ( manipulator,
            sigmacut = sigmacut, keep_predictions = True,
            force_computation_K = force_K,
            give_explanation = True )

        cr, _ = self.critic.predict_critic ( self.M, keep_predictions = True )
        ret = { "nll": None, "critic": None, "oul": None, "eul": None }

        if self.obj.slhadir == None:
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
            self.error( f"predictor for {m1},{m2} failed: {explanation}." )
        if not worked:
            return ret
        ## now get the likelihoods
        nlls={}
        ## start with the SM likelihood
        # for debugging only!
        # hasCMSSUS20004 = self.checkPredictions ()
        # ret["hasCMSSUS20004"] = hasCMSSUS20004
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
        if hasattr ( self.critic, "predictions" ):
            del self.critic.predictions

        return ret

    def checkPredictions ( self, anaId : str = "CMS-SUS-20-004" ):
        """ a debug function to find out why sometimes CMS-SUS-20-004 gets dropped
        """
        # return
        hasCMSSUS20004 = False
        for tp in self.predictor.predictions:
            my_anaId = tp.analysisId()
            if my_anaId == anaId:
                hasCMSSUS20004 = True
        #if not hasCMSSUS20004:
        #    self.error ( f"x {self.xvalue} y {self.yvalue} has no CMS-SUS-20-004" )
        return hasCMSSUS20004

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
        if self.obj.slhadir != None:
            return
        cmd = f"rm {self.M.currentSLHA}"
        subprocess.getoutput ( cmd )

    def setMass ( self, pid : Union[int,tuple], mass : float,
                  heed_partners : bool = False ):
        """ set mass or ssm of <pid> to <mass>
        :param heed_partners: if true, then also set partner particles masses
        """
        partners = [ ( 1000023, 1000024 ) ]
        self.M.masses[pid]=mass
        if not heed_partners:
            return
        """
        for pair in partners:
            if not pid in pair:
                continue
            for p in pair:
                if p in self.M.masses and self.massesAreTied ( p, pid ):
                    self.M.masses[p]=float(mass)
        """

    def setParameter ( self, pid : Union[int,tuple], value : float,
                       coord : str ):
        """
        :param coord: "x" or "y"
        """
        value = float(value)
        assert type(pid) in [ int, tuple ], "pid is neither int nor tuple"
        if type(pid)==int:
            self.setMass ( pid, value )
            return
        if coord == "y" and self.obj.y_is_dm:
            self.setDM ( pid, value, coord )
        self.setSSMultiplier ( pid, value )

    def setDM ( self, pid : tuple, dm : float, coord : str ):
        """ set the delta_m
        :param pid: e.g. ( 1000022, 1000023 )
        :param dm: delta_m
        """
        assert len(pid)==2, "setDM needs two pids"
        assert coord == "y", "y coordinates only"
        if self.obj.xvariable == pid[1]:
            self.M.masses[ pid[0] ] = self.M.masses[ pid[1] ]- dm
        else:
            self.M.masses[ pid[1] ] = self.M.masses[ pid[0] ]+ dm

    def setSSMultiplier ( self, pids : tuple, ssm : float ):
        """ set the ssm multipliers for pids to ssm.
        for now, set also for all signs
        """
        if pids in self.M.ssmultipliers:
            self.M.ssmultipliers[pids]=ssm
        if pids[1] == 1000024:
            pids1 = ( -pids[1], pids[0] )
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
        parameterpoints=self.getAllParameterPoints()
        nxvariables = len(rxvariable)
        ct = 0
        lspmass = self.M.masses[self.M.LSP]
        self.Morig = copy.deepcopy ( self.M )
        for i1,m1 in enumerate(rxvariable):
            sx = self.obj.getFullVariableName ( self.obj.xvariable, "x" )
            if type(self.obj.xvariable)==int:
                ## heed the LSP mass limit
                if m1 < lspmass and self.obj.xvariable != self.M.LSP:
                    self.pprint ( f"skipping {sx}={m1:.1f} < m(LSP)={lspmass:.1f}" )
                    continue
            thrnr = 0
            try:
                thrnr = int ( self.threadnr.replace("nll","") )
            except Exception as e:
                pass
            setnr = i1+1 + thrnr * ( nxvariables )
            self.pprint ( f"now starting with point set #{setnr} [of {nxvariables} in this thread]" )
            self.pprint ( f"this point set contains {len(ryvariable)} points" )
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
                self.M = copy.deepcopy( self.Morig )
                self.setParameter ( self.obj.xvariable, m1, "x" )
                self.setParameter ( self.obj.yvariable, m2, "y" )
                sy = self.obj.getFullVariableName ( self.obj.yvariable, "y" )
                if type(self.obj.yvariable)==int:
                    ## heed the LSP mass limit
                    if m2 < lspmass:
                        self.pprint ( f"skipping {sy}={m2:.1f} < {lspmass:.1f}" )
                        continue
                if m2 > m1 and type(self.obj.xvariable)==float and \
                        type(self.obj.yvariable) == float:
                    ## for masses we assume yvariable to be the daughter
                    self.warning ( f"{sx}={m1} < {sy}={m2}. skipping!" )
                    continue
                if m2 < 0.:
                    self.warning ( f"{sy}={m2:.1f}<0. skipping!" )
                    continue
                hasResult = self.hasResultsForPoint ( m1, m2 )
                if hasResult and not self.obj.args["redo"]:
                    self.pprint ( f"loading from cache: {sx}={m1:.2f} {sy}={m2:.2f}" )
                    continue
                point = self.getPredictions ( False, m1, m2 )
                nlls = point["nll"]
                if not nlls: continue
                nnlls,nnonzeroes=0,0

                for mu,nll in nlls.items():
                    nnlls+=len(nll)

                sx = self.obj.getFullVariableName ( self.obj.xvariable, "x" )
                sy = self.obj.getFullVariableName ( self.obj.yvariable, "y" )
                self.pprint ( f"{GREEN}{i1}/{nxvariables}{RESET}: {sx}={GREEN}{m1:.1f}{RESET}, {sy}={GREEN}{m2:.1f}{RESET}, {len(nlls)} mu's, {nnlls} nlls." )
                point["x"] = float ( m1 )
                point["y"] = float ( m2 )
                parameterpoints.append ( point )
                self.addNewPoint ( point ) ## add the point
        return parameterpoints

def runThread ( threadid: int, obj, rxvariable, ryvariable,
        return_dict : Union[Dict,None] = None ):
    """ the method needed for parallelization to work """
    col, res = f"\033[48;5;{threadid+235}m", "\033[0m"

    thread = NLLThread ( f"{col}nll{threadid}{res}", obj )
    newpoints = thread.run ( rxvariable, ryvariable )
    if return_dict != None:
        return_dict[threadid]=newpoints
    thread.clean()
    # thread.updatePickleFile()
    return newpoints

class NLLScanner ( LoggerBase ):
    """ class that encapsulates a likelihood sweep """
    def __init__ ( self, protomodel, xvariable, yvariable, nproc,
                   environ : RunEnviron, args : dict ):
        """
        :param rundir: the rundir
        :param environ: the RunEnviron
        :param skip_production: if possible, skip production, go to plotting
        :param dry_run: dont actually perform the actions
        :param output: prefix for output file [nll]
        :param y_is_dm: y variable is delta_m, not ssm
        """
        super ( NLLScanner, self ).__init__ ( "nll", "info" )
        self.args = args
        self.y_is_dm = args["yvariable"].startswith("dm")
        self.M = protomodel
        self.origmasses = copy.deepcopy ( self.M.masses )
        self.origssmultipliers = copy.deepcopy ( self.M.ssmultipliers )
        self.environ = environ
        self.xvariable = xvariable
        self.yvariable = yvariable
        x_short = self.getFullVariableName(self.xvariable,"x",True)
        y_short = self.getFullVariableName(self.yvariable,"y",True)
        picklefile = f"{self.args['output']}{x_short}_{y_short}.pcl"
        self.picklefile = picklefile
        self.nproc = nproc
        self.predictor = Predictor ( 'nll', environ=self.environ )
        self.critic = Critic ( 'nll', environ=self.environ )
        self.cprint ( "yellow", f"starting with {nproc} threads" )
        self.pprint ( f"self.predictor = Predictor ( 'nll', environ='{self.environ.runDictFile}' )" )
        yname = self.getFullVariableName ( self.yvariable, "y", True )
        xname = self.getFullVariableName ( self.xvariable, "x", True )
        # dmy = ""
        #if self.y_is_dm:
        #    yname = yname.replace("_","_dm")
        self.resultsdir = f"{self.environ.rundir}/nlls{xname}{yname}/"
        self.slhadir = None
        if args["keep_slha"]:
            self.slhadir = f"{self.environ.rundir}/slha{xname}{yname}/"
            helpers.mkdir  ( self.slhadir )

    def getVariableName ( self, variable, var_type : str ) -> str:
        """
        :returns: ssm, dm, or m
        """
        if var_type == "y":
            if type( variable )==int:
                return "m"
            if self.y_is_dm:
                return "dm"
            return "ssm"
        if type( variable )==int:
            return "m"
        return "ssm"

    def getFullVariableName ( self, variable : Union[str,int],
            var_type : str, for_filename : bool = False ) -> str:
        """
        :param var: e.g. 1000006
        :param for_filename: ssmXtXt instead of ssm(Xt,Xt)
        :returns: e.g. ssm(Xt,Xt)
        """
        prefix = self.getVariableName ( variable, var_type )
        name = namer.asciiName(variable)
        if for_filename:
            name = name.replace(",","").replace(" ","")
            name = name.replace("(","").replace(")","")
            return f"{prefix}{name}"
        ret = f"{prefix}({name})"
        return ret

    def describeRange ( self, r ):
        """ describe range r in a string """
        if len(r)==0:
            return ""
        if len(r)==1:
            return f"{r[0]:.2f}"
        if len(r)==2:
            return f"{r[0]:.2f},{r[1]:.2f}"
        return f"{r[0]:.2f},{r[1]:.2f} ... {r[-1]:.2f} -> {len(r)} points"

    def runForParameterPoints ( self, rxvariable, ryvariable ):
        """ run for the given mass points
        :param rxvariable: list of masses for xvariable
        :param ryvariable: list of masses for yvariable
        :returns: parameterpoints
        """
        if self.args["dry_run"]:
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
            ret = runThread ( 0, self, rxvariable, ryvariable )
            thread.updatePickleFile()

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
        #if yvariable != self.M.LSP:
        #    self.pprint ( f"we currently assume yvariable to be the mass of the LSP, but it is {yvariable}" )
        if os.path.exists ( self.picklefile ) and self.args["skip_production"]:
            self.pprint ( f"we were asked to skip production: {self.picklefile} exists." )
            return
        import numpy
        c = Combiner()
        anaIds = c.getAnaIdsWithPids ( self.M.bestCombo, [ xvariable, yvariable ] )
        ## mass range for xvariable
        if type(yvariable) == int:
            self.yvalue = self.M.masses[yvariable]
        if type(yvariable) == tuple:
            if self.y_is_dm:
                self.yvalue = self.M.masses[ yvariable[1] ] - self.M.masses[ yvariable[0] ]
            else:
                self.yvalue = self.M.ssmultipliers[yvariable]

        # choose the axis boundaries and step sizes such that the hiscore values
        # are nicely central
        from numpy import ceil
        ndxmin = int ( ceil (( self.xvalue - range1["min"] ) / range1["dm"]) )
        ndxmax = int ( ceil (( range1["max"] - self.xvalue ) / range1["dm"]) )
        rxvariable = numpy.arange ( self.xvalue - ndxmin*range1["dm"],
                       self.xvalue + ndxmax * range1["dm"] + 1e-5, range1["dm"] )
        rxvariable = rxvariable[rxvariable>=0.]
        ndymin = int ( ceil (( self.yvalue - range2["min"] ) / range2["dm"]) )
        ndymay = int ( ceil (( range2["max"] - self.yvalue ) / range2["dm"]) )
        ryvariable = numpy.arange ( self.yvalue - ndymin*range2["dm"],
                       self.yvalue + ndymay * range2["dm"] + 1e-5, range2["dm"] )
        ryvariable = ryvariable[ryvariable>=0.]

        # FIXME make a method that returns m(X1Z), dm(X2Z,X1Z), etc
        self.cprint ( "green", f"range for {namer.asciiName(xvariable)}: {self.describeRange( rxvariable )}" )
        self.cprint ( "green", f"range for {namer.asciiName(yvariable)}: {self.describeRange( ryvariable )}" )
        self.cprint ( "green", f"total {len(rxvariable)*len(ryvariable)} points, {nevents} events for {topo}" )
        # self.M.createNewSLHAFileName ( prefix=f"nll{xvariable}" )
        #self.M.initializePredictor()
        self.predictor.filterForTopos ( topo )
        self.M.walkerid = 2000

        thread0 = NLLThread ( "nll0", self )
        thread0.ntotal = len(rxvariable)*len(ryvariable)+1
        thread0.writeRunMeta()
        if not thread0.hasResultsForPoint ( self.xvalue, self.yvalue ):
            point = thread0.getPredictions ( False, self.xvalue, self.yvalue )
            point["x"] = float ( self.xvalue )
            point["y"] = float ( self.yvalue )
            thread0.addNewPoint ( point )
            nlls = point["nll"]
            critics = point["critic"]
            thread0.clean()
            self.pprint ( f"protomodel point: {self.getFullVariableName(self.xvariable,'x')}={self.xvalue:.2f}, {self.getFullVariableName(self.yvariable,'y')}={self.yvalue:.2f}, {len(nlls)} nlls" )
            # point [ "x" ] = float ( self.xvalue )
            # point [ "y" ] = float ( self.yvalue )
            parameterpoints = [ point ]
        else:
            parameterpoints = thread0.getAllParameterPoints()

        if False:
            ## freeze out all other particles? We shouldnt!
            for pid_,m_ in self.M.masses.items():
                if pid_ not in [ self.xvariable, self.yvariable ]:
                    self.M.masses[pid_]=1e6

        self.runForParameterPoints ( rxvariable, ryvariable )
        self.removeSLHAFile()

    def removeSLHAFile ( self ):
        if self.slhadir == None:
            self.M.delCurrentSLHA()

    def cleanFirst ( self ):
        """ clean results dir and pickle file before running """
        if os.path.exists ( self.picklefile ):
            self.pprint ( f"cleaning out {prettyFileName(self.picklefile)}" )
            os.unlink ( self.picklefile )
        if self.args["dict_file"]:
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
                 1000022: "electroweakinos,stops,TChiZISRqq,TChiISR",
                 1000024: "electroweakinos,stops,TChiZISRqq,TChiISR",
                 (1000024,1000023): "electroweakinos,stops,TChiZISRqq,TChiISR",
                 1000025: "electroweakinos,stops,TChiZISRqq,TChiISR",
                 1000035: "electroweakinos,stops,TChiZISRqq,TChiISR",
                 1000001: "T2",  1000002: "T2", 1000003: "T2", 1000004: "T2" }
        ### make the LSP scan depend on the mother
        if args.topo == None:
            args.topo = topo[args.xvariable]
        # self.xvalue = self.M.masses[self.xvariable]
        if type(self.xvariable) == int:
            self.xvalue = self.M.masses[self.xvariable]
        if type(self.xvariable) == tuple:
            if self.xvariable[0]> self.xvariable[1]:
                self.xvariable = ( self.xvariable[1], self.xvariable[0] )
            self.xvalue = self.M.ssmultipliers[self.xvariable]
        if type(self.yvariable) == int:
            self.yvalue = self.M.masses[self.yvariable]
        if type(self.yvariable) == tuple:
            if self.yvariable[0]> self.yvariable[1]:
                self.yvariable = ( self.yvariable[1], self.yvariable[0] )
            self.yvalue = self.M.ssmultipliers[self.yvariable]
        nbinsx, nbinsy = 20, 20 # how many bins do we want per dimension
        if args.minx == None:
            args.minx = self.xvalue*.6
        if args.maxx == None:
            args.maxx = self.xvalue*1.7
        if args.deltamx == None:
            args.deltamx = ( args.maxx - args.minx ) / nbinsx
        if args.miny == None:
            if type(self.yvariable) == int:
                args.miny = max ( self.yvalue*.6 - 10., 1. )
            if type(self.yvariable) == tuple:
                # args.miny = self.yvalue*.2
                args.miny = 0.
        if args.maxy == None:
            if type(self.yvariable) == int:
                args.maxy = self.yvalue*1.9 + 10.
            if type(self.yvariable) == tuple:
                args.maxy = self.yvalue*5.
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
            help='variable for y axis, e.g. 1000022 or "ssm(Xt,Xt)" or "dm(X2Z,X1Z)", [X1Z]',
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
            help='clean pickle files, results dir before running',
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
    argparser.add_argument ( '--redo',
            help="ignore cache, redo all points",
            action="store_true" )
    argparser.add_argument ( '--keep_slha',
            help="keep the SLHA files",
            action="store_true" )
    argparser.add_argument ( '--profile_mu',
            help="profile the overall signal strength parameter",
            action="store_true" )
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
    environ.pprint ( f"fetched {protomodel}\n                 from {prettyFileName(args.hiscores)}" )
    # environ.pprint ( f"ssms={protomodel.ssmultipliers}" )

    xvariables = [ namer.pid ( args.xvariable ) ]
    args.xvariable = args.xvariable.replace("ssm","")
    args.yvariable = args.yvariable.replace("ssm","")

    if args.xvariable == 0:
        xvariables = findPids( environ.rundir )
    for xvariable in xvariables:
        s_yvariable = args.yvariable.replace("dm","")
        yvariable = namer.pid ( s_yvariable )
        scanner = NLLScanner( protomodel, xvariable, yvariable, nproc,
                environ = environ, args = vars(args) )
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
