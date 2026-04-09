#!/usr/bin/env python3

""" a simple class that performs global modifications on a list of results.
Used to ``take out potential signals'' i.e. put all observations to values
expected from background, by sampling the background model. """

__all__ = [ "readDatabaseDictFile", "ExpResModifier" ]

# https://link.springer.com/content/pdf/10.1007/JHEP02(2015)004.pdf

import copy, os, sys, time, subprocess, math, numpy, shutil
import scipy.spatial
sys.path.insert( 0, "../" )
sys.path.append('../smodels')
from csetup import setup
setup()
from scipy import stats
from builder.protomodel import ProtoModel
from builder.manipulator import Manipulator
from ptools.helpers import computeP, computeZFromP, computePForDataSet, \
         computePSLv2, py_dumps
from smodels.base import runtime
if False:
    runtime._experimental["truncatedgaussians"] = True
from smodels.base.model import Model
from smodels.share.models.SMparticles import SMList
from share.model_spec import BSMList
from smodels.matching.theoryPrediction import theoryPredictionsFor
from smodels.statistics.simplifiedLikelihoods import Data, UpperLimitComputer, \
         LikelihoodComputer
from smodels.statistics.basicStats import observed, apriori, \
         aposteriori, NllEvalType
from smodels.base.physicsUnits import fb, GeV, TeV
from smodels.decomposition import decomposer
from smodels.base.smodelsLogging import logger
from smodels.experiment.databaseObj import Database
from base.loggerbase import LoggerBase
from base.runEnviron import RunEnviron
from tester.combinationsmatrix import getYamlMatrix
from typing import Dict, List, Text, Callable, Union
# from icecream import ic
from smodels_utils.helper.terminalcolors import *

from smodels.base.runtime import _deltas_rel_default
from smodels.statistics.statsTools import StatsComputer
from smodels.experiment.datasetObj import CombinedDataSet
import pyhf

logger.setLevel("ERROR")

hasWarned = { "noupperlimits": 0 }

def readDatabaseDictFile ( filename : str = "default.dict",
       filterWith : Union[None,Callable] = None ) -> Dict:
    """ read in content of a database dict file <filename>
    :param filename: the filename of the database dictionary.
    often it is <dbversion>.dict or *_database.dict or
    signal_database.dict.
    :param filterWith: optionally supply a filter function
    that takes the analysis name and the analysis dictionary as
    arguments, and is supposed to return a boolean, with true
    meaning a pass, false meaning drop

    :returns: a dictionary with 'meta', 'data', and 'basename'
    """

    with open( filename,"rt") as f:
        tmp=f.readlines()
    lines = []
    firstcommentline=None
    for i,line in enumerate(tmp):
        if line.startswith("#"):
            if firstcommentline == None:
                firstcommentline = i
            continue
        lines.append ( line )
    basename = os.path.basename ( filename ).replace(".dict","")
    meta = eval('\n'.join(lines[:firstcommentline]))
    nan=float("nan")
    inf=float("inf")
    data = eval("\n".join(lines[firstcommentline:]))
    newdata = {}
    for i,v in data.items():
        keepThis = True
        if filterWith is not None:
            keepThis = filterWith ( i, v )
        if keepThis == False:
            continue
        if "expectedBG" in v and v["expectedBG"]>=0.:
            newdata[i]=v
        else:
            if i.endswith ( ":ul" ):
                pass
            else:
                eBG=None
                if "expectedBG" in v:
                    eBG = v["expectedBG"]
    basename = os.path.basename ( filename ).replace(".dict","")
    from ptools.helpers import countObsAtExp
    meta["countObsAtExp"]=countObsAtExp
    return { "meta": meta, "data": newdata, "basename": basename }

class ExpResModifier ( LoggerBase ):
    epilog="""
Examples:
=========

Fake SM-only database:
----------------------
./expResModifier.py -R $RUNDIR -d original.pcl -s fake1

Database with a fake signal:
----------------------------
./expResModifier.py -R $RUNDIR -d original.pcl -s signal1 -P signal_model.dict

No data synthesis, just create the _database.dict file:
-------------------------------------------------------
./expResModifier.py -R $RUNDIR -d original.pcl -s original --no_synthesis

Playback the modifications described in playback file "db.dict":
----------------------------------------------------------------
WARNING this functionality has not yet been tested!
./expResModifier.py -R $RUNDIR -d original.pcl -p db.dict -o playedback.pcl

Build a database:
-----------------
./expResModifier.py -B -d ../../../smodels-database

Just filter the database:
-------------------------
./expResModifier.py -d ./original.pcl --remove_orig --nofastlim --onlyvalidated --nosuperseded --dontsample --remove_nonagg -o test.pcl

"""

    def __init__ ( self, args : dict ):
        """ constructor.

        args ( dict ):
          - database: path to database
          - max: upper limit on an individual excess
          - suffix: suffix to use, e.g. fake, signal, etc
          - lognormal: if True, use lognormal for nuisances, else Gaussian
          - fixedsignals: if True, then use the central value of theory prediction
          as the signal yield, dont draw from Poissonian
          - fixedbackgrounds: if True, then use the central value of theory
          prediction as the background yield, dont draw from Poissonian
          - seed: if int and not None, set random number seed
          - ulmassscale: maximum distance (in GeV) for the euclidean space in
          masses, for a signal to populate an UL map
        """
        super ( ExpResModifier, self ).__init__ ( "erm" )
        self.superseded = set() ## take note of everything superseded
        self.fastlim = set() # take note of everything fastlim
        self.sigNTotal = { "total": 0 } # total numbers of injected signals
        self.nMCMC_min =    500000
        self.nMCMC_max = 100000000
        self.defaults()
        if "max" in args and args["max"] == None:
            args["max"] = 100
        for a,value in args.items():
            setattr ( self, a, value )
        if "rundir" in args:
            self.rundir = setup( args["rundir"] )
        if self.outfile in [ "None", "NONE", "none" ]:
            self.outfile = None
        self.logfile = "modifier.log"
        self.startLogger()
        self.logCall()
        #self.database = setup(args["dbpath"])
        if "seed" in args:
            self.setSeed ( args["seed"] )
        self.run()

    def createMyTruth ( self ):
        """ create the truth.dict file, with the true BSM model,
        and how it would score with the new signal database """
        from multiverse.mhelpers import createMyTruthFile
        self.truth = createMyTruthFile ( signal_model = self.pmodel,
            dbpath = self.outfile, outfile = "truth.dict",
            interactive = False )

    def defaults( self ):
        """ define the defaults """
        self.db = None
        self.comments = {} ## comments on entries in dict
        self.hasFiltered = False
        self.timestamps = False
        self.no_synthesis = True
        self.protomodel = None
        self.stats = {}
        self.dbpath = "../../smodels-database"
        self.outfile = ""
        self.suffix = "db"
        self.rundir = None
        self.fudge = 1.
        self.nofastlim = False
        self.onlyvalidated = False
        self.nosuperseded = False
        self.noupperlimits = False
        self.remove_orig = False
        self.remove_nonagg = False
        self.dontsample = False
        self.lognormal = False
        self.fixedsignals = False
        self.fixedbackgrounds = False
        self.max = 100
        self.ulmassscale = 300
        self.seed = None
        self.nproc = 1
        self.pmodel = ""
        self.allowN1N1Prod = True
        self.playback = ""
        self.verbose = 0
        self.interactive = False
        self.build = False
        self.check = False
        self.compute_ps = True
        self.extract_stats = False
        self.upload = False
        self.symlink = False
        self.keep = False

    def setSeed( self, seed ):
        if seed is None:
            return
        from ptools import helpers
        helpers.seedRandomNumbers( seed )
        self.pprint ( f"setting random seed to {seed}" )

    def interact ( self, listOfExpRes ):
        import IPython
        IPython.embed( using=False )

    def extractStats ( self ):
        """ dont produce a new fake database, extract a stats dict
            from an existing database. """
        self.info ( "extracting stats" )
        picklefile = self.dbpath
        if not "/" in self.dbpath and not self.dbpath in [ "official" ]:
            picklefile = os.path.abspath ( f"{self.rundir}/{self.dbpath}" )
        if self.rundir in self.dbpath:
            picklefile = self.dbpath
        self.pprint ( f"Extracting stats from {picklefile}" )
        if self.db == None:
            combinationsmatrix, status = getYamlMatrix()
            if not combinationsmatrix or status != 0:
                logger.error("Combination matrix not loaded correctly.")
            self.pprint ( f"loading database {self.dbpath}" )
            self.db = Database ( picklefile, combinationsmatrix=combinationsmatrix )
            self.pprint ( f"loaded db v{self.db.databaseVersion}" )
        self.dbversion = self.db.databaseVersion
        listOfExpRes = self.db.expResultList
        self.stats = {}
        for expRes in listOfExpRes:
            for i,dataset in enumerate(expRes.datasets):
                dId = dataset.dataInfo.dataId
                if dId == None:
                    dId = "ul"
                label = f"{dataset.globalInfo.id}:{dId}"
                D={}
                info = dataset.dataInfo
                dt = info.dataType
                if dt == "upperLimit":
                    for txname in dataset.txnameList:
                        D[txname.txName]=list ( map ( float, txname.txnameData.y_values ) )

                for i in [ "observedN", "origN", "expectedBG", "lmbda", "bgError",
                           "origUpperLimit", "origExpectedUpperLimit", "upperLimit",
                           "expectedUpperLimit", "thirdMoment" ]:
                    if hasattr ( info, i ):
                        D[i] = getattr ( info, i )
                        if i in [ "expectedUpperLimit", "upperLimit" ]:
                            D[i]=float ( D[i].asNumber(fb) )
                if self.timestamps:
                    D["timestamp"]=dataset.globalInfo.lastUpdate
                self.addToStats ( label, D, dataset.globalInfo )

    def drawNuisance ( self, mu = 0., sigma = 1. ):
        """ draw from the nuisance model.
        :param mu: the expecation value of the distribution
        :param sigma: the standard deviation of the distribution
        :returns: a single fake observation
        """
        ## if lognormal is not selected, or mu is very small
        ## (as the lognormal cannot really produce mu=0)
        if not self.lognormal or abs(mu)<(sigma/4.):
            return stats.norm.rvs ( mu, sigma )
        loc = mu**2 / numpy.sqrt ( mu**2 + sigma**2 )
        stderr = numpy.sqrt ( numpy.log ( 1 + sigma**2 / mu**2 ) )
        ret = stats.lognorm.rvs ( s=stderr, scale=loc )
        return ret

    def computeNewObserved ( self, txname, globalInfo, x_ = None ):
        """ given expected upper limit, compute a fake observed limit
            by sampling the non-truncated Gaussian likelihood
        :param x: if not None, use it
        """
        expected = txname.txnameDataExp
        observed = txname.txnameData
        ## we only draw once for the entire UL map, equivalent to assuming
        ## that we are dealing with only one signal region
        ## second basic assumption: sigma_obs approx sigma_exp
        allpositive = False
        ctr = 0
        x = float("inf")
        ## stop when all values are positive
        while not allpositive and ctr < 10:
            ret = copy.deepcopy ( expected )
            ctr += 1
            x = float("inf")
            if x_ != None:
                x = x_
            D = {}
            ctr2 = 0
            while x > self.max and ctr2 < 10:
                x = 0.
                ctr2 += 1
                if not self.fixedbackgrounds:
                    x = self.drawNuisance() * self.fudge # draw but once from standard-normal
                # x = stats.norm.rvs() * self.fudge # draw but once from standard-normal
            D["x"] = float(x)
            D["lumi"] =float ( globalInfo.lumi * fb)
            allpositive = True
            for i,y in enumerate( expected.y_values ):
                sigma_exp = y / 1.96 ## the sigma of the Gaussian
                #D["yexp"]= y
                #D["yobs"]= float("nan")
                #if len(expected.y_values) == len(observed.y_values):
                #    D["yobs"]=observed.y_values[i]
                #D["sigma_exp"]= sigma_exp
                ## now lets shift, observed limit = expected limit + dx
                obs = y + sigma_exp * x ## shift the expected by the random fake signal
                if i == 0:
                    D["y0old"] = float ( observed.y_values[i] )
                    D["y0exp"] = float ( y )
                    D["sigma_exp0"] = float ( sigma_exp )
                    D["y0new"] = float ( obs )
                    self.comments["y0old"]="the old observed y value for first entry in UL map"
                    self.comments["y0exp"]="the expected y value for first entry in UL map"
                    self.comments["y0new"]="the fake new observed bg y value for first entry in UL map"
                    self.comments["sigma_exp0"]="the computed sigma for first entry in UL map"
                #D["y"]= obs ## we keep only last entry, but thats ok
                if obs <= 0.:
                    ## try again
                    allpositive = False
                ret.y_values[i] = obs ## now we simply shift
            if ctr > 2:
                self.log ( f"WARNING seems like I am having a hard time getting all "\
                        "values of {globalInfo.id} positive." )

        label = f"{globalInfo.id}:ul:{txname.txName}"
        D["fudge"]=self.fudge
        if self.timestamps:
            D["timestamp"]=globalInfo.lastUpdate
        self.addToStats ( label, D, globalInfo )
        self.log ( f"computed new UL result {globalInfo.id}:{txname.txName}, x={x:.2f}" )
        if x > 3.5:
            self.log ( f"WARNING high UL x={x:.2f}!!!" )
        return ret

    def bgUpperLimit ( self, dataset ):
        """ fix the upper limits, use expected (if exists) as observed """
        ## FIXME wherever possible, we should sample from the non-truncated likelihood, take that as the signal strength and re-computed a likelihood with it.
        for i,txname in enumerate(dataset.txnameList):
            if hasattr ( txname, "txnameDataExp" ) and txname.txnameDataExp != None:
                txnd = self.computeNewObserved ( txname, dataset.globalInfo )
                dataset.txnameList[i].txnameData = txnd
        return dataset

    def logCall ( self ):
        """ log how expResModifier got called """
        f=open("expResModifier.log","at")
        args = ""
        for i in sys.argv:
            if " " in i or "," in i:
                i = f'"{i}"'
            args += f"{i} "
        f.write ( f"[expResModifier.py-{time.asctime()}]\n{args.strip()}\n")
        f.close()

    def startLogger ( self ):
        subprocess.getoutput ( f"mv {self.logfile} modifier.old" )
        self.log ( f"starting at {time.asctime()} with zmax of {self.max}" )
        self.log ( f"arguments were {' '.join ( sys.argv )}" )

    def finalize ( self ):
        """ finalize, delete files, create my.truth """
        # print ( "[expResModifier] finalize" )
        if self.keep:
            pass
        elif hasattr ( self, "protomodel" ) and self.protomodel is not None and \
                type(self.protomodel) != str:
            self.protomodel.delCurrentSLHA()
        self.createNewRunDict()
        self.createMyTruth()

    def createNewRunDict ( self ):
        """ create the new run.dict file, referencing the signal database now """
        self.environ.moveRunDict ( "run_creation.dict" )
        args = vars(self.environ)
        self.dbversion = self.db.databaseVersion
        dbpath = self.outfile
        if dbpath in [ "none", None ]:
            ## if no outfile is given, we revert to the previous dbpath
            dbpath = self.dbpath
        newargs = { "dbpath": dbpath, "dbversion": self.dbversion }
        for i in [ "allowN1N1Prod", ]:
            newargs[i]=args[i]
        ## this line triggers creation of the new run.dict
        self.environ = RunEnviron.new ( **newargs )

    def produceProtoModel ( self, filename : str, dbversion : str,
           allowN1N1Prod : bool = True ) -> Union[None,ProtoModel]:
        """ try to produce a protomodel from pmodel
        :param filename: filename of pmodel dictionary
        :param dbversion: version of database, for tracking
        :param allowN1N1Prod: if bool, then have also N1N1 production
        :returns: none if not succesful, else protomodel object
        """
        self.environ = RunEnviron.new( allowN1N1Prod = allowN1N1Prod,
                dbversion = dbversion, dbpath = self.dbpath )
        if filename == "":
            return None
        if not os.path.exists ( filename ):
            self.pprint ( f"When trying to construct protomodel, {filename} does not exist" )
            return None
        shutil.copyfile ( filename, f"{self.rundir}/my.signal" )
        walkerid = "erm"
        expected = False
        select = "all"
        keep_meta = True
        # M = ProtoModel ( walkerid, self.dbpath, expected, select, keep_meta )
        ## create a new environment, possibly overwriting old run.dicts
        M = ProtoModel ( walkerid, keep_meta, environ = self.environ )
        M.createNewSLHAFileName ( prefix="erm" )
        ma = Manipulator ( M, walkerid = walkerid, environ = self.environ )
        with open ( filename, "rt" ) as f:
            try:
                m = eval ( f.read() )
            except (SyntaxError,TypeError) as e:
                print ( f"[expResModifier] error parsing {filename}: {e}" )
                print ( f"[expResModifier] is this a protomodel with 'masses', etc defined?" )
                sys.exit()
        if allowN1N1Prod and not (1000022,1000022) in m['ssmultipliers']:
            self.warn ( f"we allow N1N1 production but not N1N1 production in signal model" )
        if not allowN1N1Prod and (1000022,1000022) in m['ssmultipliers']:
            self.warn ( f"we disallow N1N1 production but N1N1 production in signal model" )
        ma.initFromDict ( m, initTestStats=True )
        ma.M.computeXSecs( keep_slha = True )
        self.log ( f"xsecs produced {ma.M.currentSLHA}" )
        self.log ( f" `- does currentslha exist? {os.path.exists ( ma.M.currentSLHA )}" )
        self.pprint ( f"BSM model's xsecs:" )
        ma.printXSecs( useParticleNames = True )
        self.protomodel = ma.M
        return self.protomodel

    def removeEmpty ( self, listOfExpRes ):
        ret = []
        for er in listOfExpRes:
            hasEntry = False
            dses = []
            for dataset in er.datasets:
                txnames = [ tx.txName for tx in dataset.txnameList ]
                if len(txnames)>0:
                    hasEntry = True
                    dses.append ( dataset )
                else:
                    self.info( f"{er.globalInfo.id}:{dataset.dataInfo.dataId} has only empty txnames. will remove." )
                er.datasets = dses
            if hasEntry:
                ret.append ( er )
            else:
                self.info ( f"{er.globalInfo.id} has only empty datasets, will remove." )
        return ret

    def removeMLModels ( self, updatedListOfExpRes : list, listOfExpRes : list ):
        """ remove the ML models from updatedListOfExpRes.
        In a more refined version, we can actually compare against listOfExpRes,
        to retain situations where nobs has not changed """
        if self.no_synthesis:
            return # no need to remove
        for er in updatedListOfExpRes:
            if hasattr ( er.globalInfo, "mlModels" ):
                del er.globalInfo.mlModels

    def modifyDatabase ( self ):
        """ modify the database, possibly write out to a pickle file
        :param outfile: if not empty, write the database into file
        :param suffix: suffix to append to database version
        :param pmodel: if not empty, then this is the file name of the signal
                       model. in this case fake a signal
        :returns: the database
        """
        spmodel = f"protomodel is '{self.pmodel}'"
        if self.pmodel == "":
            spmodel = "no protomodel given"
        self.info ( f"starting to create {self.outfile} from {self.dbpath}. suffix is '{self.suffix}', {spmodel}." )
        if self.db == None:
            combinationsmatrix, status = getYamlMatrix()
            if not combinationsmatrix or status != 0:
                logger.error("Combination matrix not loaded correctly.")
            self.pprint ( f"loading database {self.dbpath}" )
            self.db = Database ( self.dbpath, combinationsmatrix=combinationsmatrix)
            self.pprint ( f"loaded db v{self.db.databaseVersion}" )
        self.dbversion = self.db.databaseVersion
        self.orig_dbversion = self.dbversion
        listOfExpRes = self.removeEmpty ( self.db.expResultList ) ## seems to be the safest bet?
        self.produceProtoModel ( self.pmodel, self.db.databaseVersion,
                                 self.allowN1N1Prod )
        self.log ( f"{len(listOfExpRes)} results before faking bgs" )
        updatedListOfExpRes = self.fakeBackgrounds ( listOfExpRes )
        self.log ( f"{len(updatedListOfExpRes)} results after faking bgs" )
        updatedListOfExpRes = self.addSignals ( updatedListOfExpRes )
        self.log ( f"{len(updatedListOfExpRes)} results after adding signals" )
        self.removeMLModels ( updatedListOfExpRes, listOfExpRes )
        if hasattr ( self.db, "subs" ): ## for smodels 2.1
            self.db.subs[0].expResultList = updatedListOfExpRes
            self.db.subs = [ self.db.subs[0] ]
        else:
            self.db.expResultList = updatedListOfExpRes
        newver = self.db.databaseVersion
        if self.suffix is not None:
            newver = self.db.databaseVersion + "_" + self.suffix
        self.db.txt_meta.databaseVersion = newver
        self.db.pcl_meta.databaseVersion = newver
        self.pprint ( f"Constructed fake database with {len(updatedListOfExpRes)} (of {len(listOfExpRes)}) results" )
        self.createBinaryFile()
        return self.db

    def computeP ( self, obsN : float, bgExp : float, bgErr : float,
            thirdMoment : Union[float,None] ):
        """ a convenience function for computation of p values """
        if thirdMoment is None:
            p = computeP ( obsN, bgExp, bgErr, nmin = self.nMCMC_min,
                           nmax = self.nMCMC_max )
            return p
        p = computePSLv2 ( obsN, bgExp, bgErr, thirdMoment,
                nmin = self.nMCMC_min, nmax = self.nMCMC_max )
        return p

    def computePForDataSet ( self, dataset,
            obsN : Union[int,None] = None )-> float:
        """ convenience function to compute p for dataset, with right
        nmin and nmax """
        return computePForDataSet ( dataset, obsN, nmin = self.nMCMC_min,
                                    nmax = self.nMCMC_max )

    def sampleEfficiencyMap ( self, dataset ):
        """ for the given dataset,
        sample from background and put the value as observed """
        orig = dataset.dataInfo.observedN
        exp = dataset.dataInfo.expectedBG
        thirdMoment = None
        if hasattr ( dataset.dataInfo, "thirdMoment" ):
            thirdMoment = dataset.dataInfo.thirdMoment * self.fudge**3
        err = 0.
        if not self.fixedbackgrounds:
            err = dataset.dataInfo.bgError * self.fudge
        D = { "origN": int(orig), "expectedBG": exp, "bgError": err,
              "fudge": self.fudge, "lumi": float(dataset.globalInfo.lumi * fb) }
        porig = self.computePForDataSet ( dataset )
        self.checkIfZero ( porig, dataset )
        D["orig_p"]=porig
        self.comments["orig_p"]="p-value (Gaussian nuisance) of original observation (no fudge factor applied)"
        origZ = computeZFromP ( porig )
        D["orig_Z"]=origZ
        self.comments["orig_Z"]="the significance Z of the original observation (no fudge factor applied)"
        label = f"{dataset.globalInfo.id}:{dataset.dataInfo.dataId}"
        txnames = [ tx.txName for tx in dataset.txnameList ]
        txnames.sort()
        if len ( txnames ) == 0:
            self.warning ( f"no txnames for {label}." )
        D["txns"]=tuple(txnames )
        self.comments["txns"]="tuple of txnames that populate this signal region / analysis"
        if self.timestamps:
            D["timestamp"]=dataset.globalInfo.lastUpdate
        constraints = set()
        for txni in dataset.txnameList:
            constraints.add ( txni.constraint )
        D["constraints"]=tuple( constraints )
        self.comments["constraints"]="tuple of the sms constraints"

        if thirdMoment is not None:
            D["thirdMoment"]=thirdMoment
            self.comments["thirdMoment"]="third moment for SLv2 likelihoods"
        if self.no_synthesis:
            D["newObs"] = D["origN"]
            D["new_p"]=D["orig_p"]
            D["new_Z"]=D["orig_Z"]
            self.comments["newObs"]="the new fake observation (signal + background) -- in our case same as 'origN'"
            self.comments["new_p"]="p-value (Gaussian nuisance) of newObs -- in our case same as 'orig_p'"
            self.comments["new_Z"]="significance (Gaussian nuisance) of newObs -- in our case same as 'orig_Z'"
            self.addToStats ( label, D, dataset.globalInfo )
            return dataset
        if self.compute_ps:
            p = self.computeP ( orig, exp, err, thirdMoment )
            self.checkIfZero ( p, dataset )
            if abs(self.fudge-1.)>1e-10:
                self.comments["orig_p_fudged"]="p-value (Gaussian nuisance) of original observation (with fudge factor applied -- is this useful?)"
                D["orig_p_fudged"]=p
                origZ = computeZFromP ( p )
                D["orig_Z_fudged"]=origZ
                self.comments["orig_Z_fudged"]="the significance Z of the original observation (with fudge factor applied -- is this useful?)"
        Z = float("inf")
        ct = 0
        while Z > self.max and ct < 10:
            ct += 1
            # lmbda = stats.norm.rvs ( exp, err )
            lmbda = exp
            if not self.fixedbackgrounds:
                lmbda = self.drawNuisance ( exp, err )
            dataset.dataInfo.lmbda = lmbda
            if hasattr ( dataset.dataInfo, "thirdMoment" ):
                thirdMoment = float ( dataset.dataInfo.thirdMoment ) * self.fudge**3
                lmbda += thirdMoment * (lmbda-exp)**2 / err**4
            if lmbda < 0.:
                lmbda = 0.
                # lmbda = [self.drawNuisance ( exp, err ) for _ in range(1000)]
                # lmbda = numpy.mean([v if v > 0 else 0 for v in lmbda])
            obs = lmbda
            toterr = 0.
            if not self.fixedbackgrounds:
                obs = stats.poisson.rvs ( lmbda )
                toterr = math.sqrt ( err**2 + exp )
            if True: # toterr > 0.:
                pnew = self.computeP ( orig, exp, err, thirdMoment )
                if thirdMoment is not None:
                    D["thirdMoment"]=thirdMoment
                    self.comments["thirdMoment"]="third moment for SLv2 likelihoods"
                self.checkIfZero ( pnew, dataset )
                Z = - scipy.stats.norm.ppf ( pnew )
                # Z = ( obs - exp ) / toterr
                # origZ = ( orig - exp ) / toterr
            if Z < self.max or ct > 9:
                self.log ( f"effmap replacing old nobs={orig} (bg={exp:.2f}+/-{err:.2f}, lmbda={lmbda:.2f}, Z={Z:.2f}) with nobs={obs} for {dataset.globalInfo.id}:{dataset.dataInfo.dataId}" )
                dataset.dataInfo.observedN = obs
        if Z > 3.5:
            self.log ( f"WARNING!!! high em Z={Z:.2f}!!!!" )
        D["Zbg"]=float(Z)
        self.comments["Zbg"]="the significance of the observation, bg only"
        #D["Z"]=float(Z)
        #self.comments["Z"]="the significance of the observation, taking into account the signal"
        self.comments["lmbda"]="Poissonian lambda of the fake background"
        D["lmbda"]=float(lmbda)
        D["newObs"]=obs
        self.comments["newObs"]="the new fake observation (signal + background)"
        if self.compute_ps:
            p = self.computeP ( obs, exp, err, thirdMoment )
            self.checkIfZero( p, dataset )
            self.comments["new_p"]="p-value (Gaussian nuisance) of newObs"
            self.comments["new_Z"]="significance (Gaussian nuisance) of newObs"
            D["new_p"]=p
            newZ = computeZFromP ( p )
            D["new_Z"]=newZ
        D["obsBg"]=obs
        self.comments["obsBg"]="the new fake observation, background component"
        D["toterr"]=toterr
        ## origN stores the n_observed of the original database
        dataset.dataInfo.origN = orig
        self.addToStats ( label, D, dataset.globalInfo )
        return dataset

    def createEMStatsDict ( self, dataset ) -> Dict:
        """ given the dataset, create a stats dictionary, for SL and pyhf
        datasets """
        orig = dataset.dataInfo.observedN
        exp = dataset.dataInfo.expectedBG
        thirdMoment = None
        if hasattr ( dataset.dataInfo, "thirdMoment" ):
            thirdMoment = dataset.dataInfo.thirdMoment * self.fudge**3
        err = 0.
        if not self.fixedbackgrounds:
            err = dataset.dataInfo.bgError * self.fudge
        D = { "origN": int(orig), "expectedBG": exp, "bgError": err, "fudge": self.fudge,
              "lumi": float(dataset.globalInfo.lumi * fb) }
        if thirdMoment is not None:
            D["thirdMoment"]=thirdMoment
            self.comments["thirdMoment"]="third moment for SLv2 likelihoods"
        if self.compute_ps:
            p = self.computePForDataSet ( dataset )
            D["orig_p"]=p
            self.comments["orig_p"]="p-value (Gaussian nuisance) of original observation (no fudge factor applied)"
            origZ = computeZFromP ( p )
            self.comments["orig_Z"]="the significance Z of the original observation (no factor applied)"
            D["orig_Z"]=origZ
            p = self.computeP ( orig, exp, err, thirdMoment )
            self.comments["orig_p_fudged"]="p-value (Gaussian nuisance) of original observation (fudge factor applied)"
            self.checkIfZero ( p, dataset )
            D["orig_p_fudged"]=p
            origZ = computeZFromP ( p )
            D["orig_Z_fudged"]=origZ
            self.comments["orig_Z_fudged"]="the significance Z of the original observation (fudge factor applied)"
        txnames = [ tx.txName for tx in dataset.txnameList ]
        txnames.sort()
        if len ( txnames ) == 0:

            self.warning ( f"no txnames for {label}." )
        D["txns"]=tuple(txnames )
        self.comments["txns"]="tuple of txnames that populate this signal region / analysis"
        constraints = set()
        for txni in dataset.txnameList:
            constraints.add ( txni.constraint )
        D["constraints"]=tuple( constraints )
        self.comments["constraints"]="tuple of the sms constraints"
        # D["txns"]=",".join(txnames )
        if self.timestamps:
            D["timestamp"]=dataset.globalInfo.lastUpdate
        return D

    def getPyhfname ( self, dataset ):
        for jsonfile, SRs in dataset.globalInfo.jsonFiles.items():
            for sr in SRs:
                if sr["smodels"] == dataset.dataInfo.dataId:
                    return sr["pyhf"]
        return None

    def addSignalForPyhf ( self, dataset, sigN ):
        """ add sigN to the json file in dataset.globalInfo.jsons """
        pyhfname = self.getPyhfname ( dataset )
        if pyhfname == None:
            print ( f"ERROR no pyhfname!!!!" )
            import sys, IPython; IPython.embed( colors = "neutral" ); sys.exit()
        for i,json in enumerate(dataset.globalInfo.jsons):
            for obs in json["observations"]:
                if obs["name"] == pyhfname:
                    oldBG = obs["data"][0]
                    # print ( f"[expResModifier] adding {sigN} to {oldBG} in {pyhfname}" )
                    if len(obs["data"])>1:
                        self.error ( f"more than one bin!!! {obs['data']}" )
                        import sys, IPython; IPython.embed( colors = "neutral" ); sys.exit()
                    obs["data"][0]+=sigN

    def checkIfZero ( self, p, dataset ):
        """ give a warning if a p-value is zero zero """
        if p > 1e-100:
            return
        self.warn ( f"{dataset.globalInfo.id}:{dataset.dataInfo.dataId} has p={p} -- maybe you injected too strong a signal?" )

    def addSignalForEfficiencyMap ( self, dataset, tpred, lumi ):
        """ add a signal to this efficiency map. background sampling is
            already taken care of """
        txns = list ( map ( str, tpred.txnames ) )
        txns.sort()
        self.log ( f"add EM matching tpred {tpred.analysisId()}/{tpred.dataId()} {','.join(txns)}: {tpred.xsection.asNumber(fb):.2g} fb" )
        label = f"{dataset.globalInfo.id}:{dataset.dataInfo.dataId}"
        if not label in self.stats:
            if not "CR" in label:
                print ( f"[expResModifier] {label} not found in stats! (dunno if that is a problem)" )
            return dataset
        orig = dataset.dataInfo.observedN
        sigLambda = float ( tpred.xsection * lumi )
        D={}
        ## FIXME sigLumi and sigXSec, we should take it back out later
        D["sigLumi"] = float ( lumi.asNumber(1./fb) )
        self.comments["sigLumi"]="the lumi assumed for the signal in fb FIXME remove"
        D["sigXSec"] = float ( tpred.xsection.asNumber(fb) )
        self.comments["sigXSec"]="the fiducial xsec computed for the signal in fb FIXME remove"
        D["sigLambda"]=sigLambda
        self.comments["sigLambda"]="the lambda for the signal"
        sigN = sigLambda
        if self.fixedsignals:
            self.comments["sigN"]="the number of events from the added signal (using central value)"
        else:
            sigN = stats.poisson.rvs ( sigLambda )
        D["sigN"]=0
        if "sigN" in self.stats[label]:
            ## sigN is the total number of added signals
            ## they may be from multiple topologies
            D["sigN"]=self.stats[label]["sigN"]
        D["sigN"]=D["sigN"]+sigN
        self.sigNTotal["total"] += sigN
        self.comments["sigN"]="the number of events from the added signal"
        D["sigTxns"] = txns
        self.comments["sigTxns"]="the txnames that contributed to sigN"
        txnsc = "_".join( txns )
        ## sigNT<x> denotes the contributions from the individual theory preds
        sigNt = f"sigN{txnsc}"
        D[ sigNt ] = sigN
        if not sigNt in self.sigNTotal:
            self.sigNTotal[sigNt]=0
        self.sigNTotal[sigNt]+=sigN
        D["obsBg"]=self.stats[label]["newObs"]
        err = dataset.dataInfo.bgError * self.fudge
        dataset.dataInfo.sigN = sigN ## keep track of signal
        if sigN == 0:
                self.log ( f" `- signal sigN={sigN} re obsN={orig} too small. skip.")
                dataset.dataInfo.origUpperLimit = dataset.dataInfo.upperLimit
                dataset.dataInfo.origExpectedUpperLimit = dataset.dataInfo.expectedUpperLimit
                D["newObs"]=orig
                self.addToStats ( label, D, dataset.globalInfo )
                return dataset
        self.log ( f" `- effmap adding sigN={sigN} to obsN={orig} -> newObs={orig+sigN}" )
        dataset.dataInfo.trueBG = orig ## keep track of true bg
        dataset.dataInfo.observedN = orig + sigN
        if hasattr ( dataset.globalInfo, "jsons" ):
            self.addSignalForPyhf ( dataset, sigN )

        D["newObs"]=dataset.dataInfo.observedN
        exp = dataset.dataInfo.expectedBG
        err = dataset.dataInfo.bgError * self.fudge
        #toterr = math.sqrt ( err**2 + exp )
        #Z = 0.
        #if toterr > 0.:
        #    Z = ( dataset.dataInfo.observedN - exp ) / toterr
        #D["Z"]=Z
        #self.comments["Z"]="the significance of the observation, taking into account the signal"
        thirdMoment = None
        if hasattr ( dataset.dataInfo, "thirdMoment" ):
            thirdMoment = dataset.dataInfo.thirdMoment
        new_p = self.computeP ( dataset.dataInfo.observedN, exp, err,
                                thirdMoment )
        self.checkIfZero( new_p, dataset )
        new_Z = computeZFromP ( new_p )
        D["new_p"] = new_p
        D["new_Z"] = new_Z
        ## now recompute the limits!!
        if orig == 0.0:
            orig = 0.00001
        m = Data( orig+sigN, orig, err**2, nsignal = 1. )
        computer = UpperLimitComputer( LikelihoodComputer ( m ) )
        lumi = dataset.globalInfo.lumi# .asNumber(1./fb)
        maxSignalXsec = computer.getUpperLimitOnMu ( ) / lumi
        dataset.dataInfo.origUpperLimit = dataset.dataInfo.upperLimit
        dataset.dataInfo.origExpectedUpperLimit = dataset.dataInfo.expectedUpperLimit
        dataset.dataInfo.upperLimit = maxSignalXsec
        maxSignalXsec = computer.getUpperLimitOnMu( evaluationType=apriori ) / lumi #  NllEvalType.apriori ) #/ lumi
        dataset.dataInfo.expectedUpperLimit = maxSignalXsec
        self.addToStats ( label, D, dataset.globalInfo )
        return dataset

    def txNameIsIn ( self, txname, tpred ):
        """ check if txname is in tpred
        :param txname: a txName object
        :param tpred: a theoryPred object
        """
        for txn in tpred.txnames:
            if txn.txName == txname.txName:
                return True
        return False

    def addToStats ( self, label, Dict, globalInfo ):
        """ add the content of dictionary Dict to the stats,
            under the label "label" """
        if hasattr ( globalInfo, "supersedes" ):
            self.superseded.add ( globalInfo.supersedes )
        if hasattr ( globalInfo, "supersededBy" ):
            self.superseded.add ( globalInfo.id )
        if hasattr ( globalInfo, "contact" ) and "fastlim" in globalInfo.contact:
            self.fastlim.add ( globalInfo.id )
        if hasattr ( globalInfo, "comment" ) and "fastlim" in globalInfo.comment:
            self.fastlim.add ( globalInfo.id )
        if not label in self.stats:
            # we dont yet have an entry, so lets start
            self.stats[label]=Dict
            return
        # we have an entry, so we add
        for k,v in Dict.items():
            self.stats[label][k]=v

    def addSupersededFlags ( self ):
        """ at the end, add superseded flag to everything superseded. """
        for k,v in self.stats.items():
            p1 = k.find(":")
            name = k[:p1]
            if name in self.superseded:
                self.stats[k]["superseded"]=True
            if name in self.fastlim:
                self.stats[k]["fastlim"]=True

    def distance ( self, v1, v2 ):
        """ compute distance between v1 and v2 """
        ret = 0.
        nmin = min ( len(v1), len(v2) )
        nmax = max ( len(v1), len(v2) )
        div = nmax / nmin
        v1,v2 = list(v1)[:nmin],list(v2)[:nmin]
        #if len(v1)*2 == len(v2):
        #    v1 = v1*2
        sums = []
        for _1,_2 in zip ( v1, v2 ):
            sums.append ( ( _1 - _2 )**2 )
        ret = math.sqrt ( sum(sums) / div )
        return ret

    def addSignalFromDict ( self, txname, dataset, values ):
        """ add a signal to this UL result. background sampling is
            already taken care of """
        self.pprint ( f"warning, signal playback not yet testededed for ULs" )
        from smodels.base.physicsUnits import fb
        from ptools import helpers
        txns = values["txns"]
        ## so we simply add the theory predicted cross section to the limit
        sigmaN = values["sigmaN"] # tpred.xsection.asNumber(fb)
        label = f"{dataset.globalInfo.id}:ul:{txns}"
        D={}
        D["sigmaN"]=sigmaN
        D["pids"]=list(set(values["pids"]))
        D["masses"]=values["masses"]
        D["txns"]=tuple(txns)
        self.comments["txns"]="tuple of txnames that populate this signal region / analysis"
        self.comments["sigmaN"]="the added theory prediction (in fb), for UL maps"
        ## sigmaN is the predicted production cross section of the signal,
        ## in fb
        if not txname.txName in txns:
            return txname.txnameData
        hasAdded = 0
        txnd = txname.txnameData
        etxnd = txname.txnameDataExp
        masses = values["masses"]
        # coordsTpred = txnd.PCAtransf( masses ) # , txnd._V, txnd.delta_x ) ## coordinates of tpred
        minDist = float("inf") ## for the closest point we store the numbers
        gaussAtMass = scipy.stats.norm.pdf ( 0, 0, self.ulmassscale ) * 1.2
        for yi,y in enumerate(txnd.y_values):
            pt = txnd.tri.points[yi] ## the point in the rotated coords
            pt_masses = txnd.inversePCAtransf ( pt )
            dist = self.distance ( masses, pt_masses )
            # dist = self.distance ( pt, coordsTpred )
            if dist > self.ulmassscale: ## change y_values only in vicinity of protomodel
                continue
            oldv = txnd.y_values[yi]
            oldo = txnd.y_values[yi]
            hasExpected=False
            sqrts = float ( dataset.globalInfo.sqrts.asNumber(TeV))
            mysignal = self.computeXSecForMass ( sigmaN, masses, pt_masses,
                                                 D["pids"], sqrts )
            mysignal = mysignal * scipy.stats.norm.pdf ( dist, 0, self.ulmassscale ) / gaussAtMass
            if etxnd != None and len(txnd.y_values) == len(etxnd.y_values):
                dt = ( ( txnd.delta_x - etxnd.delta_x )**2 ).sum()
                if dt < 1e-2:
                    hasExpected=True
                    oldv = etxnd.y_values[yi] ## FIXME more checks pls
            if dist < minDist:
                ## remember the candidate
                minDist = dist
                D["yold"]=oldo
                D["dist"]=dist
                self.comments["dist"]="distance of closest point to protomodel"
                if hasExpected:
                    D["yexp"]=oldv
                    self.comments["yexp"]="expected y value (fb) closest to signal protomodel for UL map"
                self.comments["yold"]="old y value (fb) closest to signal protomodel for UL map"
                self.comments["ynew"]="new y value (fb) closest to signal protomodel for UL map"
                D["ynew"]=oldv+mysignal
            txnd.y_values[yi]=oldv + mysignal
            hasAdded += 1
            if hasAdded == 0:
                self.pprint ( "warning: signal was not added in {tpred.analysisId()}:{txname.txName}" )
            D[f"signalpoints{txname.txName}"]=hasAdded
            D[f"totalpoints{txname.txName}"]=len(txnd.y_values)
            self.comments["signalpointsTx"]="number of grid points that got the signal injected"
            self.comments["totalpointsTx"]="total number of grid points in that map"
            self.addToStats ( label, D, dataset.globalInfo )
        return txnd

    def getMassVector ( self, tpred ):
        """ get the particle masses of a theory prediction """
        ## FIXME this will have to be smarter
        ret = []
        if tpred is None:
            return ret
        sms=tpred.smsList[0]
        point = tpred.txnames[0].getDataFromSMS ( sms )
        return point

    def getPIDVector ( self, tpred ):
        """ get the particle pdgs of a theory prediction """
        ## FIXME this will have to be smarter
        ret = set()
        if tpred is None:
            return ret
        for sms in tpred.smsList:
            for node in sms.nodes:
                if type(node.particle.pdg) == int and node.particle.pdg > 99:
                    ret.add ( node.particle.pdg )
        return list(ret)

    def computeXSecForMass ( self, sigmaN, oldmasses, newmasses, pids, sqrts ):
        """ given the cross section at mass oldmass,
        compute a rough equivalent cross section at mass point newmass,
        for pids, according to how reference cross sections scale """
        from ptools.xsecFit import XSecFitter
        mypids = pids[:2]
        mypids.sort()
        mypids = tuple ( mypids )
        func = XSecFitter(mypids, sqrts)
        oldmass = oldmasses[0]*GeV ## FIXME what about asymmetric?
        oldxs = func.getValueFromFit ( oldmass, inverse=False )
        newmass = newmasses[0]*GeV
        newxs = func.getValueFromFit ( newmass, inverse=False )
        if False: # newmass == 350*GeV:
            self.pprint ( f"computeXSecForMass {oldmasses} {newmasses} {pids}" )
            self.pprint ( f"           -- oldxs {oldxs} newxs {newxs}" )
        if oldxs is None or newxs is None:
            return float ( sigmaN )
        return float ( sigmaN * newxs / oldxs )

    def addSignalForULMap ( self, dataset, tpred, lumi ):
        """ add a signal to this UL result. background sampling is
            already taken care of """
        if tpred is None:
            return dataset
        from smodels.base.physicsUnits import fb
        from ptools import helpers
        txns = list ( map ( str, tpred.txnames ) )
        txns.sort()
        self.log ( f"add UL matching tpred {tpred.analysisId()}: <{tpred.xsection.asNumber(fb):.2g}> fb {tpred.smsList} {','.join(txns)}" )
        ## so we simply add the theory predicted cross section to the limit
        sigmaN = tpred.xsection.asNumber(fb)
        label = f"{tpred.analysisId()}:ul:{','.join(txns)}"
        D={}
        D["sigmaN"]=sigmaN
        D["pids"]=self.getPIDVector ( tpred )
        # D["smsList"]=tpred.smsList
        D["masses"]=self.getMassVector ( tpred )
        D["txns"]=tuple(txns)
        # D["txns"]=",".join(txns)
        self.comments["txns"]="tuple of txnames that populate this signal region / analysis"
        self.comments["sigmaN"]="the added theory prediction (in fb), for UL maps"
        ## sigmaN is the predicted production cross section of the signal,
        ## in fb
        gaussAtMass = scipy.stats.norm.pdf ( 0, 0, self.ulmassscale ) * 1.2
        for i,txname in enumerate(dataset.txnameList):
            if not self.txNameIsIn ( txname, tpred ):
                continue
            hasAdded = 0
            txnd = txname.txnameData
            etxnd = txname.txnameDataExp
            masses = self.getMassVector ( tpred )
            # coordsTpred = txnd.PCAtransf ( masses ) # , txnd._V, txnd.delta_x ) ## coordinates of tpred
            minDist, minPt = float("inf"),None ## for the closest point we store the numbers
            for yi,y in enumerate(txnd.y_values):
                pt = txnd.tri.points[yi] ## the point in the rotated coords
                ## the masses in the original frame
                pt_masses = txnd.inversePCAtransf ( pt )
                dist = self.distance ( masses, pt_masses )
                # rot_dist = self.distance ( pt, coordsTpred )
                if dist < minDist: ## just so we know how far away we are
                    minDist = dist
                    minPt = txnd.inversePCAtransf ( pt ) # , txnd._V, txnd.delta_x )
                # if we want to add it "binarically"
                #if dist > self.ulmassscale: ## change y_values only in vicinity of protomodel
                #    continue
                oldv = txnd.y_values[yi]
                oldo = txnd.y_values[yi]
                mysignal = self.computeXSecForMass  ( sigmaN, masses, pt_masses,
                   D["pids"], float(tpred.dataset.globalInfo.sqrts.asNumber(TeV)) )
                # lets make it gaussian
                mysignal = mysignal * scipy.stats.norm.pdf ( dist, 0, self.ulmassscale ) / gaussAtMass
                hasExpected=False
                if etxnd != None and len(txnd.y_values) == len(etxnd.y_values):
                    dt = ( ( txnd.delta_x - etxnd.delta_x )**2 ).sum()
                    if dt < 1e-2:
                        hasExpected=True
                        oldv = etxnd.y_values[yi] ## FIXME more checks pls
                if dist < minDist:
                    ## remember the candidate
                    D["yold"]=oldo
                    D["dist"]=dist
                    self.comments["dist"]="distance of closest point to protomodel"
                    if hasExpected:
                        D["yexp"]=oldv
                        self.comments["yexp"]="expected y value (fb) closest to signal protomodel for UL map"
                    self.comments["yold"]="old y value (fb) closest to signal protomodel for UL map"
                    self.comments["ynew"]="new y value (fb) closest to signal protomodel for UL map"
                    D["ynew"]=oldv+mysignal
                # self.pprint ( f"adding {mysignal}*fb to {oldv}*fb in {pt_masses} for {txname}::{tpred.analysisId()}" )
                txnd.y_values[yi]=oldv + mysignal
                hasAdded += 1
            if hasAdded == 0:
                self.pprint ( f"warning: no signal was added in {tpred.analysisId()}:{txname.txName}, closest was point {minPt} at d={minDist:.2f}" )
            D[f"signalpoints{txname.txName}"]=hasAdded
            D[f"totalpoints{txname.txName}"]=len(txnd.y_values)
            self.comments["signalpointsTx"]="number of grid points that got the signal injected"
            self.comments["totalpointsTx"]="total number of grid points in that map"
            self.addToStats ( label, D, dataset.globalInfo )
            dataset.txnameList[i].txnameData = txnd
            dataset.txnameList[i].sigmaN = sigmaN
        return dataset

    def writeStats ( self, statsname : os.PathLike = None ):
        """ write out the collected stats
        :param statsname: sth like *_database.dict
        """
        if self.suffix in [ None, "None", "", "none" ]:
            return
        filename = f"{self.rundir}/{self.suffix}.dict".replace("//","/")
        if statsname != None:
            filename = statsname
        meta = { "orig_database": self.orig_dbversion,
                 "orig_dbpath": self.dbpath,
                 "Zmax": self.max,
                 "fake_database": self.dbversion,
                 "fudge": self.fudge,
                 "timestamp": time.asctime(),
                 "allowN1N1Prod": self.allowN1N1Prod,
                 "lognormal": self.lognormal,
                 "ulmassscale": self.ulmassscale,
                 "fixedsignals": self.fixedsignals,
                 "bsm_file": self.pmodel,
                 "fixedbackgrounds": self.fixedbackgrounds
        }
        if "K" in self.truth:
            meta["K_true"]=self.truth["K"]
        if "TL" in self.truth:
            meta["TL_true"]=self.truth["TL"]
        if self.protomodel != None:
            meta["signal_model"]=self.protomodel.dict()
        #meta["protomodel"]=None
        #if self.protomodel!= None:
        #    meta["protomodel"] = f'{str(self.protomodel)}'
        if hasattr ( runtime, "_drmax" ):
            meta["_drmax"]=runtime._drmax
        if hasattr ( runtime, "_experimental" ):
            meta["_experimental"]=runtime._experimental
        self.addSupersededFlags()
        meta["sigNTotal"] = self.sigNTotal
        self.pprint ( f"saving stats to {filename}" )
        self.log ( f"saving stats to {filename}" )
        with open ( filename, "wt" ) as f:
            ds = py_dumps ( meta, indent = 4 )
            ds = ds.replace( "false", "False" ).replace ( "true", "True" )
            # ds = ds.replace( "inf", "float('inf')" )
            ds = ds.replace( r'"\"None\""', 'None')
            f.write ( ds + "\n"  )
            # f.write ( f"{meta!s}\n" )
            f.write ( f"# this file was created with {' '.join(sys.argv)}\n" )
            if len(self.comments)>0:
                f.write ( "# explanations on the used variables:\n" )
                f.write ( "# =====================================\n" )
            else:
                f.write ( "# no explanations for variables have been given\n" )
            for k,v in self.comments.items():
                f.write ( f"# {k}: {v}\n" )
            ds = py_dumps ( self.stats, indent=4 )
            # ds = ds.replace( "inf", "float('inf')" )
            f.write ( ds+ "\n" )
            f.close()

    def produceTopoList ( self ):
        """ create smstopolist """
        from smodels.base.physicsUnits import fb, GeV
        model = Model ( BSMList, SMList )
        model.updateParticles ( inputFile=self.protomodel.currentSLHA )
        mingap=10*GeV
        sigmacut = 0.*fb
        self.topos = decomposer.decompose ( model, sigmacut, minmassgap=mingap )

    def addSignalsSingleProc ( self, listOfExpRes : list ) -> list:
        """ thats the method that adds a typical signal """
        if self.protomodel == None:
            return listOfExpRes
        self.produceTopoList()
        ctr,els = 0, ""
        for topo in self.topos:
            for sms in self.topos[topo]:
                ctr+= 1
                els += f"{sms!s}, "
            #for el in topo.elementList:
            #    ctr+=1
            #    els += str(el) + ", "
        els=els[:-2]
        self.log ( f"now add the signals from {self.getPModelName()}, {ctr} topologies: {els}" )
        addedUL, addedEM = 0, 0
        print ( f"[expResModifier.addSignalsSingleProc] {len(listOfExpRes)} results get signal added: ", end="" )
        for l,expRes in enumerate(listOfExpRes):
            self.db.selectExpResults ( analysisIDs = [ expRes.globalInfo.id ] )
            print ( ".", flush=True, end="" )
            tpreds = theoryPredictionsFor ( self.db, self.topos,
                    useBestDataset=False, combinedResults=False )
            if tpreds == None:
                continue
            lumi = expRes.globalInfo.lumi
            for i,dataset in enumerate(expRes.datasets):
                dt = dataset.dataInfo.dataType
                dsname = dataset.dataInfo.dataId
                if dt == "upperLimit":
                    for tpred in tpreds:
                        if tpred.dataId() == None:
                            # IPython.embed()
                            addedUL += 1
                            listOfExpRes[l].datasets[i] = self.addSignalForULMap ( dataset, tpred, lumi )
                else:
                    for tpred in tpreds:
                        if tpred.dataId() == dsname:
                            addedEM += 1
                            listOfExpRes[l].datasets[i] = self.addSignalForEfficiencyMap ( dataset, tpred, lumi )
                    ## expRes.datasets[i] = self.fixUpperLimit ( dataset )
        self.db.selectExpResults ( analysisIDs = [ "all" ] )
        print ( )
        self.log ( f"added {addedUL} UL signals and {addedEM} EM signals" )
        return listOfExpRes

    def signalAdder ( self, listOfExpRes ):
        # print ( "[expResModifier] signalAdder", len(listOfExpRes), len(self.topos) )
        addedUL, addedEM = 0, 0
        for l,expRes in enumerate(listOfExpRes):
            print ( ".", flush=True, end="" )
            tpreds = theoryPredictionsFor ( expRes, self.topos, useBestDataset=False,
                                            combinedResults=False )
            if tpreds == None:
                # ret.append ( expRes )
                continue
            lumi = expRes.globalInfo.lumi
            for i,dataset in enumerate(expRes.datasets):
                dt = dataset.dataInfo.dataType
                dsname = dataset.dataInfo.dataId
                if dt == "upperLimit":
                    for tpred in tpreds:
                        if tpred.dataId() == None:
                            # IPython.embed()
                            addedUL += 1
                            listOfExpRes[l].datasets[i] = self.addSignalForULMap ( dataset, tpred, lumi )
                else:
                    for tpred in tpreds:
                        if tpred.dataId() == dsname:
                            addedEM += 1
                            listOfExpRes[l].datasets[i] = self.addSignalForEfficiencyMap ( dataset, tpred, lumi )
                    ## expRes.datasets[i] = self.fixUpperLimit ( dataset )
        self.log ( f"added {addedUL} UL signals and {addedEM} EM signals" )
        return listOfExpRes

    def getPModelName ( self ):
        """ name of protomodel """
        pmodelname = str(self.protomodel)
        for i in [ "<sub>", "<sup>", "</sub>", "</sup>" ]:
            pmodelname = pmodelname.replace( i, "" )
        return pmodelname

    def addSignals ( self, listOfExpRes ):
        """ thats the method that adds a typical signal, parallel version
        :param nproc: number of processes
        """
        if self.protomodel == None:
            return listOfExpRes
        if self.nproc == 1:
            return self.addSignalsSingleProc ( listOfExpRes )
        # print ( "adding signals", os.path.exists ( self.protomodel.currentSLHA ) )
        ret = []
        self.produceTopoList()
        self.log ( f"now add the signals from {self.getPModelName()}, K='?', {len(self.topos)} topos, {self.nproc} procs" )
        import multiprocessing
        ## listOfExpRes=listOfExpRes[:10]
        chunks = [ listOfExpRes[i::self.nproc] for i in range(self.nproc) ]
        pool = multiprocessing.Pool ( processes = self.nproc )
        tmp = pool.map ( self.signalAdder, chunks )
        print ( "done! now collect." )
        ret = []
        for t in tmp:
            for x in t:
                ret.append ( x )
        ## print ( "ret=", ret )
        return ret

    def fakeBackgroundsForSL ( self, expRes ):
        """ synthesize fake observations by sampling a simplified likelihood
        model
        :param expRes: the experimental result to do this for
        """
        # self.error ( f"FIXME fake SL backgrounds for {expRes.globalInfo.id}" )
        import numpy as np
        import scipy.stats
        covm = np.array ( expRes.globalInfo.covariance )
        if abs ( self.fudge - 1 ) > 1e-8:
            covm = self.fudge**2 * covm
        # diag = np.array ([expRes.globalInfo.covariance[i][i] for i in range(len(covm))])
        expectedBGs = np.array([ x.dataInfo.expectedBG for x in expRes.datasets ] )
        observed = np.array([ x.dataInfo.observedN for x in expRes.datasets ] )
        zeroes = np.array ( [0.]*len(covm) )
        rvs = scipy.stats.multivariate_normal.rvs ( zeroes, covm )
        ## FIXME SLv2 needed also!
        thirdMoments = None
        tpe = "SLv1"
        if hasattr ( expRes.datasets[0].dataInfo, "thirdMoment" ):
            tpe = "SLv2"
            thirdMoments = [ x.dataInfo.thirdMoment * self.fudge**3 for x in expRes.datasets ]
        self.comments["type"]="result type (None, SLv1, SLv2, pyhf)"
        from smodels.statistics.simplifiedLikelihoods import Data
        data = Data ( observed, expectedBGs, covm, thirdMoments )
        for i,dataset in enumerate(expRes.datasets):
            newObs = dataset.dataInfo.observedN
            if not self.no_synthesis:
                lmbda = dataset.dataInfo.expectedBG + rvs[i]
                if tpe == "SLv2":
                    lmbda = data.A[i] + rvs[i] + data.C[i] * rvs[i]**2 / data.B[i]**2
                    # lmbda += dataset.dataInfo.thirdMoments / diag[i]**2 * rvs[i]**2
                lmbda = max ( 0., lmbda )
                newObs = scipy.stats.poisson.rvs ( lmbda )
            D = self.createEMStatsDict ( dataset )
            D["type"]=tpe
            if self.fixedbackgrounds:
                D["newObs"]=dataset.dataInfo.expectedBG
            else:
                D["newObs"]=newObs
                if not self.no_synthesis:
                    D["lmbda"]=float(lmbda)
            if self.compute_ps:
                if self.no_synthesis:
                    D["new_p"] = D["orig_p"]
                    D["new_Z"] = D["orig_Z"]
                    self.comments["new_p"]="p-value (Gaussian nuisance) of newObs -- in our case same as 'orig_p'"
                    self.comments["new_Z"]="significance (Gaussian nuisance) of newObs -- in our case same as 'orig_Z'"
                else:
                    p = self.computePForDataSet ( dataset, newObs )
                    self.checkIfZero ( p, dataset )
                    self.comments["new_p"]="p-value (Gaussian nuisance) of newObs"
                    D["new_p"]=p
                    newZ = computeZFromP ( p )
                    self.comments["new_Z"]="significance (Gaussian nuisance) of newObs"
                    D["new_Z"]=newZ
                    if p == 0 or newZ == float("inf"):
                        self.error ( f"we got p={p} Z={newZ}. exit" )
                        sys.exit()
            expRes.datasets[i].dataInfo.observedN = newObs
            label = f"{dataset.globalInfo.id}:{dataset.dataInfo.dataId}"
            self.addToStats ( label, D, dataset.globalInfo )


    def getChannelNames ( self, channels : List[Text] ) -> List:
        """ get the names of channels from the pyhf entry """
        ret = []
        for channel in channels:
            if ";" in channel:
                [ ret.append(x) for x in channel.split(";") ]
            else:
                ret.append ( channel )
        return ret

    def fudgePyhfModel ( self, expRes, computer ):
        """ fudge the pyhf model, ie multiply all errors with self.fudge """
        anaId = expRes.globalInfo.id
        # self.error ( f"FIXME fudge factors not yet implemented for pyhf ({anaId})" )
        # import sys, IPython; IPython.embed( colors = "neutral" ); sys.exit()
        ## FIXME this needs more thought: which errors get rescaled, which dont, etc
        for iws,ws in enumerate(computer.likelihoodComputer.workspaces):
            for ich,channel in enumerate(ws["channels"]):
                for ism,sample in enumerate(channel["samples"]):
                    if not "modifiers" in sample:
                        continue
                    for im,modifier in enumerate(sample["modifiers"] ):
                        if not "data" in modifier:
                            continue
                        if modifier["data"] is None:
                            continue
                        data = modifier["data"]
                        if "hi" in data:
                            center = (data["hi"]+data["lo"])/2.
                            delta = data["hi"] - center
                            data["hi"] = center + delta * self.fudge
                            data["lo"] = center - delta * self.fudge
                        if "hi_data" in data:
                            for idd,(hi,lo) in enumerate ( zip ( data["hi_data"],data["lo_data"] ) ):
                                center = (hi+lo)/2.
                                delta = hi - center
                                data["hi_data"][idd] = center + delta * self.fudge
                                data["lo_data"][idd] = center - delta * self.fudge

    def replaceObservation ( self, expRes, sr, newObs, ws_i ):
        jsonEntries = expRes.globalInfo.jsons[ ws_i ]["observations"]
        for i,jsonEntry in enumerate(jsonEntries):
            idx = 0
            pyhfbasename = sr["pyhf"]
            data = jsonEntry["data"]
            if len(data)>1:
                p1 = sr["pyhf"].find("[")
                pyhfbasename = sr["pyhf"][:p1]
                idx = int ( sr["pyhf"][p1+1:-1] )
            oldE = data[idx]
            if jsonEntry["name"]==pyhfbasename:
                # print ( f"[expResModifier] replacing {oldE} with {newObs} in {jsonEntry} {expRes.globalInfo.id}" )
                expRes.globalInfo.jsons[ ws_i ]["observations"][i]["data"][idx]=newObs
                continue

    def fakeBackgroundsForPyhf ( self, expRes ):
        """ synthesize fake observations by sampling a pyhf model
        :param expRes: the experimental result to do this for
        """
        datasetDict= { ds.getID(): ds for ds in expRes.origdatasets }
        ## store original values
        # origN = { k : v.dataInfo.observedN for k,v in datasetDict.items() }
        srNsigDict = {ds.getID() : 0.0 for ds in expRes.origdatasets}
        # Update with theory predictions
        #srNsigDict.update({pred.dataset.getID() :
        #              (pred.xsection*pred.dataset.getLumi()).asNumber()
        #              for pred in self.datasetPredictions})

        #create combined dataset for pyhf pred
        cdataset = CombinedDataSet ( expRes )
        computer = StatsComputer.forPyhf( cdataset, srNsigDict,
                _deltas_rel_default )
        if abs ( self.fudge - 1. ) > 1e-5:
            self.fudgePyhfModel ( expRes, computer )
        srs_in_workspaces = list(expRes.globalInfo.jsonFiles.values())
        anaId = expRes.globalInfo.id

        for ws_i, (ws, srs) in enumerate(zip(
                    computer.likelihoodComputer.workspaces, srs_in_workspaces) ):
            ## srs are the names of the signal regions
            try:
                model = ws.model()
            except pyhf.exceptions.InvalidModel as e:
                print ( f"[expResModifier] pyhf.InvalidModel for {anaId} [{list(expRes.globalInfo.jsonFiles.keys())[ws_i]}][{ws_i}]: {e}" )
                continue
                # sys.exit(-1)
            channelnames = self.getChannelNames ( model.config.channels )
            pars_bkg = model.config.suggested_init()
            pars_bkg[model.config.poi_index] = 0.0 ## background
            pdf_bkg = model.main_model.make_pdf(pyhf.tensorlib.astensor(pars_bkg))
            # pdf_bkg = model.make_pdf(pyhf.tensorlib.astensor(pars_bkg))
            sample = pdf_bkg.sample()
            sampleDictPyhf = {}
            ##first we need a dictionary to translate the SR names
            #ic ( channelnames )
            sampleidx=0
            for channelname in channelnames:
                if model.config.channel_nbins [ channelname] == 1:
                    sampleDictPyhf[channelname]= float ( sample[sampleidx] )
                    sampleidx+=1
                else:
                    for i in range ( model.config.channel_nbins[channelname]):
                        fullname = f"{channelname}[{i}]"
                        sampleDictPyhf[fullname]= float ( sample[sampleidx] )
                        sampleidx+=1
            sampleDictSModelS = {}

            for sr in srs:
                if sr["pyhf"] in sampleDictPyhf:
                    sampleDictSModelS[ sr["smodels"] ] = sampleDictPyhf[ sr["pyhf"] ]
            for sr in srs:
                srname = sr["smodels"]
                if sr["type"] != "SR":
                    continue
                dataset = datasetDict[srname]
                D = self.createEMStatsDict ( dataset )
                newObs = int( sampleDictSModelS[ srname ] )
                if self.fixedbackgrounds:
                    D["newObs"]=dataset.dataInfo.expectedBG
                else:
                    D["newObs"]=newObs
                if self.no_synthesis:
                    D["newObs"]=D["origN"]
                D["type"]="pyhf"
                if self.compute_ps:
                    p = self.computePForDataSet ( dataset, newObs )
                    self.comments["new_p"]="p-value (Gaussian nuisance) of newObs"
                    D["new_p"]=p
                    newZ = computeZFromP ( p )
                    self.comments["new_Z"]="significance (Gaussian nuisance) of newObs"
                    D["new_Z"]=newZ
                    if p == 0 or newZ == float("inf"):
                        self.error ( f"we got p={p} Z={newZ}. exit" )
                        sys.exit()
                    if self.no_synthesis:
                        D["new_p"] = D["orig_p"]
                        D["new_Z"] = D["orig_Z"]
                        self.comments["new_p"]="p-value (Gaussian nuisance) of newObs -- in our case same as 'orig_p'"
                        self.comments["new_Z"]="significance (Gaussian nuisance) of newObs -- in our case same as 'orig_Z'"
                label = f"{anaId}:{dataset.dataInfo.dataId}"
                self.addToStats ( label, D, dataset.globalInfo )
                ## as the very last measure, we replace the observation with
                ## the fake observation
                dataset.dataInfo.observedN = newObs
                self.replaceObservation ( expRes, sr, newObs, ws_i )
                # this replaces the observation in the json with the new bg



    def fakeBackgrounds ( self, listOfExpRes ):
        """ thats the method that samples the backgrounds """
        ret = []
        self.log ( "now fake backgrounds" )
        for expRes in listOfExpRes:
            t0 = time.time()
            if hasattr ( expRes.globalInfo, "covariance" ):
                self.fakeBackgroundsForSL ( expRes )
            elif hasattr ( expRes.globalInfo, "jsonFiles" ):
                self.fakeBackgroundsForPyhf ( expRes )
            else:
                for i,dataset in enumerate(expRes.datasets):
                    dt = dataset.dataInfo.dataType
                    if dt == "upperLimit":
                        expRes.datasets[i] = self.bgUpperLimit ( dataset )
                    elif dt == "efficiencyMap":
                        expRes.datasets[i] = self.sampleEfficiencyMap ( dataset )
                    else:
                        print ( f"[expResModifier] dataset type {dt} unknown" )
            ret.append ( expRes )
            t1 = time.time()
            dt = t1 - t0
            if dt > 3:
                self.pprint ( f"{expRes.globalInfo.id} took {t1-t0:.2f}s" )
        self.log ( "done faking the backgrounds" )
        return ret

    def cleanTxNameData ( self, txnd ):
        txnd.y_values=numpy.array ( txnd.y_values, dtype=numpy.float32 )

        if txnd.dimensionality == 1:
            return txnd
        txnd.tri._points = numpy.array ( txnd.tri._points, dtype=numpy.float32 )
        return txnd

    def filter ( self, exclude_anas: list = [] ):
        """ filter the list fo experimental results.
        :param outfile: store result in outfile (a pickle file)
        :param nofastlim: remove fastlim results
        :param onlyvalidated: remove non-validated results
        :param nosuperseded: remove superseded results
        :param remove_orig: remove original values
        :param remove_nonagg: remove non-aggregated results
        """
        if not ( self.nofastlim or self.onlyvalidated or self.nosuperseded or self.remove_orig or self.remove_nonagg or self.noupperlimits ):
            if exclude_anas == []:
                return
        self.log ( f"starting to filter {self.outfile}. suffix is {self.suffix}." )
        if self.db == None:
            combinationsmatrix, status = getYamlMatrix()
            if not combinationsmatrix or status != 0:
                logger.error("Combination matrix not loaded correctly.")
            self.db =Database ( self.dbpath, combinationsmatrix=combinationsmatrix)
        listOfExpRes = self.db.expResultList ## seems to be the safest bet?
        if self.remove_nonagg:
            from smodels_utils.helper.databaseManipulations import filterNonAggregatedFromList
            n = len(listOfExpRes )
            listOfExpRes = filterNonAggregatedFromList ( listOfExpRes, verbose=self.verbose )
            print ( f"[expResModifier] nonaggregated filter: from {n} to {len(listOfExpRes)}" )
            if len(listOfExpRes ) <  n:
                self.hasFiltered = True
        newList = []
        for er in listOfExpRes:
            anaId = er.globalInfo.id
            addThisOne = True
            if self.nofastlim:
                if hasattr ( er.globalInfo, "contact" ) and "fastlim" in er.globalInfo.contact:
                    self.pprint ( f" `- skipping fastlim {anaId}" )
                    addThisOne = False
                    self.hasFiltered = True
            if self.nosuperseded:
                if hasattr ( er.globalInfo, "supersededBy" ):
                    self.pprint ( f" `- skipping superseded {anaId}" )
                    addThisOne = False
                    self.hasFiltered = True
            if self.noupperlimits:
                if er.datasets[0].getID() == None:
                    hasWarned["noupperlimits"]+=1
                    if hasWarned["noupperlimits"]<4:
                        self.pprint ( f" `- skipping UL-type {anaId}" )
                    if hasWarned["noupperlimits"]==4:
                        self.pprint ( f" (quenching more of the msgs given above)" )
                    addThisOne = False
                    self.hasFiltered = True
            if hasattr ( er.globalInfo, "private" ) and er.globalInfo.private in [ "True", True ]:
                    self.pprint ( f" `- skipping private {anaId}" )
                    addThisOne = False
                    self.hasFiltered = True
            import fnmatch
            for exclude in exclude_anas:
                if fnmatch.fnmatch ( anaId, exclude ):
                    if exclude == anaId:
                      print ( f"dropping {anaId}" )
                    else:
                        print ( f"dropping {anaId}: matches {exclude}" )
                    addThisOne = False
                    self.hasFiltered = True
                    break
            if not addThisOne:
                self.hasFiltered = True
                continue
            if self.onlyvalidated:
                newDs = []
                hasIssued = 0
                for ds in er.datasets:
                    txnew = []
                    for txn in ds.txnameList:
                        if txn.validated == False:
                            if hasIssued == 0:
                                self.pprint ( f" `- skipping non-validated {txn.txName}/{ds.dataInfo.dataId}/{anaId}" )
                            if hasIssued == 1:
                                self.pprint ( " `- (suppressed more, similar messages)" )
                            hasIssued += 1
                            self.hasFiltered = True
                        else:
                            txnew.append ( txn )
                    ds.txnameList = txnew
                    if len(txnew)>0:
                        newDs.append ( ds )
                er.datasets = newDs
                if len(newDs) == 0:
                    addThisOne = False
            if self.remove_orig:
                from smodels.experiment.txnameObj import TxNameData
                TxNameData._keep_values = False
                for label in [ "prettyName", "arxiv", "publication", "implementedBy",\
                               "lastUpdate", "contact" ]:
                    if hasattr ( er.globalInfo, label ):
                        delattr ( er.globalInfo, label )
                        self.hasFiltered = True
                for iD,ds in enumerate(er.datasets):
                    for it,txn in enumerate(ds.txnameList):
                        #txn.txnameData = self.cleanTxNameData ( txn.txnameData )
                        for label in [ "figureUrl", "dataUrl" ]:
                            if hasattr ( txn, label ):
                                self.hasFiltered = True
                                delattr ( txn, label )
                        if hasattr ( txn.txnameData, "origdata" ):
                            del er.datasets[iD].txnameList[it].txnameData.origdata
                            self.hasFiltered = True
                        if txn.txnameDataExp != None:
                            #txn.txnameDataExp = self.cleanTxNameData ( txn.txnameDataExp )
                            if hasattr ( txn.txnameDataExp, "origdata" ):
                                del er.datasets[iD].txnameList[it].txnameDataExp.origdata
                                self.hasFiltered = True
            if not addThisOne:
                continue
            newList.append ( er )
        self.db.subs[0].expResultList = newList
        self.createBinaryFile()

    def createBinaryFile ( self ):
        """ write binary pickle database to self.outfile """
        if self.outfile == "":
            return
        if self.suffix in [ "None", "none", "", None ]:
            return
        if self.outfile in [ None ]:
            self.pprint ( f"creation of pickle file was suppressed" )
            return
        self.pprint ( f"writing to {self.outfile}" )
        self.db.createBinaryFile( self.outfile )

    def playback ( self, playbackdict ):
        """ playback the mods described in playbackdict """
        self.pprint ( "WARNING: playback functionality has not yet been validated!!" )
        with open ( playbackdict, "rt" ) as h:
            lines = h.readlines()
            h.close()
        ## first line goes directly into database
        line = lines.pop(0)
        D = eval ( line )
        if self.db == None:
            combinationsmatrix, status = getYamlMatrix()
            if not combinationsmatrix or status != 0:
                logger.error("Combination matrix not loaded correctly.")
            self.db =Database ( self.dbpath, combinationsmatrix=combinationsmatrix)
        for k,v in D.items():
            if k in [ "dbpath", "database" ]:
                continue
            setattr ( self, k, v )
        ## now the remaining lines
        cleaned = []
        for line in lines:
            if line.startswith("#"):
                continue
            cleaned.append ( line )
        D = eval ( "\n".join ( cleaned ) )

        self.dbversion = self.db.databaseVersion
        self.lExpRes = self.db.expResultList ## seems to be the safest bet?
        # self.lExpRes = db.getExpResults ( [ "CMS-SUS-19-006" ] ) ## for debugging
        for anaids,values in D.items():
            #if not "CMS-SUS-19-006" in anaids: # for debugging
            #    continue
            self.playbackOneItem ( anaids, values )
        self.db.expResultList = self.lExpRes
        self.db.dbpath = self.outfile
        self.dbversion = f"{self.dbversion}.playedback"
        self.db.txt_meta.databaseVersion = f"{self.db.databaseVersion}.playedback"
        self.db.pcl_meta.databaseVersion = f"{self.db.databaseVersion}.playedback"
        self.createBinaryFile()

    def playbackOneItem ( self, anaids : str, values : dict ):
        """ play back a single item
        :param anaids: e.g. "CMS-SUS-14-021:ul:T2bbWWoff"
        :param values: e.g. xxx
        """
        mytxname = anaids.split(":")[-1]
        # print ( f"playing back {anaids} for ", mytxname, values )
        self.log ( f"playing back {anaids}" )
        tokens = anaids.split(":")
        anaid = tokens[0]
        isEffMap = False
        if len(tokens)==3:
            # dataType = tokens[1]
            txname = tokens[2]
        if len(tokens)==2:
            isEffMap = True
            sr=tokens[1]
        for ier,er in enumerate(self.lExpRes):
            tanaid = er.globalInfo.id
            if tanaid != anaid:
                continue
            tdatasets = er.datasets
            for ids,tds in enumerate(tdatasets):
                if tds.getType() == "upperLimit" and not isEffMap:
                    ### update an UL dataset
                    for itx,txnd in enumerate(tds.txnameList):
                        if txnd.txName != mytxname:
                            continue
                        if hasattr ( txnd, "txnameDataExp" ) and txnd.txnameDataExp != None:
                            self.pprint ( "updating UL map", tds.globalInfo.id )
                            if not "x" in values:
                                self.pprint ( f"error, cannot find x value in {tanaid}:{txnd.txName}: {values}" )
                                continue
                            ## print ( "playing back", tds.globalInfo.id, values["x"], txnd.txName, "txns" in values )
                            ntxnd = self.computeNewObserved ( txnd, tds.globalInfo, values["x"] )
                            txnd.txnameData=ntxnd
                            self.lExpRes[ier].datasets[ids].txnameList[itx].txnameData=ntxnd
                            if "sigmaN" in values and "masses" in values:
                                ntxnd = self.addSignalFromDict ( txnd, tds, values )
                                self.lExpRes[ier].datasets[ids].txnameList[itx].txnameData=ntxnd


                if tds.getType() == "efficiencyMap" and isEffMap and sr == tds.getID():
                    ### update an EM dataset
                    self.pprint ( "found EM to update", tds.getID() )
                    self.lExpRes[ier].datasets[ids].dataInfo.observedN = values["newObs"]

    def upload( self ):
        import filecmp
        # cmd = f"cp {args.outfile} ./modifier.log {self.rundir}"
        for f in [ args.outfile, self.logfile ]:
            if not filecmp.cmp ( f, f"{self.rundir}/{os.path.basename(f)}" ):
                cmd = f"cp {f} {self.rundir}"
                a = subprocess.getoutput ( cmd )
                print ( "[expResModifier]", cmd, a )
        fname = f"{self.rundir}/default.pcl"
        if os.path.exists ( fname ):
            cmd = f"rm {fname}"
            a = subprocess.getoutput ( cmd )
            print ( "[expResModifier]", cmd, a )
        cmd = f"ln -s {self.rundir}/{args.outfile} {self.rundir}/default.pcl"
        a = subprocess.getoutput ( cmd )
        print ( "[expResModifier]", cmd, a )

    def symlink ( self, outfile ):
        """ create a symlink to rundir/default.pcl """
        dest = f"{self.rundir}/default.pcl"
        if os.path.exists ( dest ):
            cmd = f"rm {dest}"
            subprocess.getoutput ( cmd )
        cmd = f"ln -s {outfile} {dest}"
        subprocess.getoutput ( cmd )

    def check ( self ):
        """ check the picklefile """
        picklefile = self.outfile
        print ( "now checking the modified database" )
        self.db = Database ( picklefile )
        # listOfExpRes = db.getExpResults()
        listOfExpRes = self.db.expResultList ## seems to be the safest bet?
        for er in listOfExpRes:
            datasets = er.datasets
            for ds in datasets:
                txnl = ds.txnameList
                for txn in txnl:
                    x = txn.txnameData.dataType
        print ( "we're good", self.db.databaseVersion )

    def run ( self ):
        if self.fixedbackgrounds and not self.fixedsignals:
            self.pprint ( "WARNING fixing backgrounds but not signals. Sounds weird" )
        if self.fixedbackgrounds and self.fudge > 1e-2:
            self.pprint ( "WARNING fixing backgrounds but fudge factor is not zero. Sounds weird" )
        if self.build:
            from smodels.experiment.txnameObj import TxNameData
            TxNameData._keep_values = True
            from smodels.experiment.databaseObj import Database
            #self.database = "official"
            # self.database = "../../smodels-database"
            self.pprint ( f"starting to build database at {self.dbpath}." )
            combinationsmatrix, status = getYamlMatrix()
            if not combinationsmatrix or status != 0:
                logger.error("Combination matrix not loaded correctly.")
            self.pprint ( f"loading database {self.dbpath}" )
            db = Database ( self.dbpath, combinationsmatrix=combinationsmatrix )
            self.orig_dbversion = db.databaseVersion
            self.pprint ( f"built database at {self.dbpath}." )
            # sys.exit()
        if self.rundir == None:
            self.rundir = os.getcwd()
        self.pprint ( f"{GREEN}rundir is {os.getcwd()}{RESET}" )
        if type(self.rundir)==str and not "/" in self.rundir and \
                not self.rundir.startswith("."):
            self.rundir = f"{os.environ['HOME']}/{self.rundir}"
        statsname = f"{self.suffix}_database.dict"
        if self.outfile is not None:
            if self.outfile == "":
                self.outfile = f"{self.suffix}.pcl"
            if self.playback not in [ None, "" ]:
                self.playback ( self.playback, self.outfile )
                statsname = "playback.dict"

            if not self.outfile.endswith(".pcl") and self.outfile != None:
                self.pprint ( f"warning, shouldnt the name of your outputfile ``{self.outfile}'' end with .pcl?" )
        #else: # outfile is None
        #    statsname = None
        self.filter ( )
        if self.dontsample:
            self.pprint ( "we were asked to not sample, so we exit now." )
            sys.exit()

        if self.extract_stats:
            er = self.extractStats()
        else:
            if not self.playback:
                er = self.modifyDatabase ( )

        if self.check:
            self.check ( )

        if self.interactive:
            self.interact ( er )

        if self.upload:
            self.upload()

        if self.symlink:
            self.symlink ( )

        self.finalize()

        if statsname is not None:
            self.writeStats( statsname )

if __name__ == "__main__":
    import argparse
    from argparse import RawTextHelpFormatter
    argparser = argparse.ArgumentParser(
                        description='Experimental results modifier. Used to synthesize fake data by setting all observations to values sampled from the background models. Can insert signals, too.', formatter_class = RawTextHelpFormatter,
                        epilog=ExpResModifier.epilog )
    argparser.add_argument ( '-d', '--dbpath',
            help='database to use [../../smodels-database]',
            type=str, default="../../smodels-database" )
    argparser.add_argument ( '-o', '--outfile',
            help='file to write out database pickle. If left empty, then outfile is <suffix>.pcl. if "none", then dont create pickle file [""]',
            type=str, default="" )
    argparser.add_argument ( '-s', '--suffix',
            help='suffix for database version, if None or "" then do not write out ["fake1"]',
            type=str, default="fake1" )
    argparser.add_argument ( '-R', '--rundir',
            help='override rundir [None]',
            type=str, default=None )
    argparser.add_argument ( '-f', '--fudge',
            help='fudge factor. all systematic errors will be multiplied by that [1.0]',
            type=float, default=1.0 )
    argparser.add_argument ( '--nofastlim',
            help='remove fastlim results',
            action='store_true' )
    argparser.add_argument ( '--onlyvalidated',
            help='remove non-validated results',
            action='store_true' )
    argparser.add_argument ( '--nosuperseded',
            help='remove superseded results',
            action='store_true' )
    argparser.add_argument ( '--no_synthesis',
            help='no data synthesis, just report observations',
            action='store_true' )
    argparser.add_argument ( '--noupperlimits',
            help='remove upper limit results',
            action='store_true' )
    argparser.add_argument ( '--remove_orig',
            help='remove original values',
            action='store_true' )
    argparser.add_argument ( '--remove_nonagg',
            help='remove nonaggregated results',
            action='store_true' )
    argparser.add_argument ( '--dontsample',
            help='do not sample at all, only filter',
            action='store_true' )
    argparser.add_argument ( '--disallowN1N1Prod',
            help='turn off N1N1 production',
            action='store_true' )
    argparser.add_argument ( '--allowN1N1Prod',
            help='explicitly turn on N1N1 production (on per default)',
            action='store_true' )
    argparser.add_argument ( '-l', '--lognormal',
            help='use lognormal, not Gaussian for nuisances (1d regions only)',
            action='store_true' )
    argparser.add_argument ( '--fixedsignals',
            help='fix the contributions from the signals, dont draw from Poissonian',
            action='store_true' )
    argparser.add_argument ( '--fixedbackgrounds',
            help='fix the contributions from the backgrounds, use central values for all.',
            action='store_true' )
    argparser.add_argument ( '-M', '--Zmax',
            help='upper limit on significance of individual excess [None]',
            type=float, default=None )
    argparser.add_argument ( '--ulmassscale',
            help='mass scale (GeV) for adding the signal in the UL maps [300.]',
            type=float, default=300. )
    argparser.add_argument ( '--seed',
            help='set a random number seed [None]',
            type=int, default=None )
    argparser.add_argument ( '-N', '--nproc',
            help='number of parallel processes, for signal adding [1]',
            type=int, default=1 )
    argparser.add_argument ( '-P', '--pmodel',
            help='supply filename of a pmodel, in which case create a signal-infused database [""]',
            type=str, default="" )
    argparser.add_argument ( '-p', '--playback',
            help='playback the modifications described in given dictionary file [""]',
            type=str, default="" )
    argparser.add_argument ( '-v', '--verbose',
            help='print results to stdout', action='store_true' )
    argparser.add_argument ( '-I', '--interactive',
            help='interactive mode', action='store_true' )
    argparser.add_argument ( '-B', '--build',
            help='build the original pickle file with all relevant info, then exit (use --dbpath to specify path)', action='store_true' )
    argparser.add_argument ( '-c', '--check',
            help='check the pickle file <outfile>', action='store_true' )
    ## turn this on always
    #argparser.add_argument ( '-C', '--compute_ps',
    #        help='compute p-values for all SRs', action='store_true' )
    argparser.add_argument ( '-t', '--timestamps',
            help='add time-stamps (only to be used with -C)', action='store_true' )
    argparser.add_argument ( '-x', '--extract_stats',
            help='dont create new database, extract stats from existing database',
            action='store_true' )
    argparser.add_argument ( '-u', '--upload',
            help='upload to $RUNDIR', action='store_true' )
    argparser.add_argument ( '-S', '--symlink',
            help='symlink default.pcl to <outfile> (in rundir)', action='store_true' )
    argparser.add_argument ( '-k', '--keep',
            help='keep temporary files (for debugging)', action='store_true' )
    args = argparser.parse_args()
    vargs = vars(args)
    if vargs["allowN1N1Prod"] and vargs["disallowN1N1Prod"]:
        print ( f"[ExpResModifier:...] you at the same time allow and disallow N1N1 prod. fix this." )
        sys.exit()
    vargs["allowN1N1Prod"]= not vargs["disallowN1N1Prod" ]
    vargs.pop ( "disallowN1N1Prod" )
    modifier = ExpResModifier( vargs )
