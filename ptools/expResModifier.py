#!/usr/bin/env python3

""" a simple class that performs global modifications on a list of results.
Used to ``take out potential signals'' i.e. put all observations to values
expected from background, by sampling the background model. """

__all__ = [ "readDictFile", "ExpResModifier" ]

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
from ptools.helpers import computeP, computeZFromP, computePForDataSet, computePSLv2
from smodels.base import runtime
if False:
    runtime._experimental = True
from smodels.base.model import Model
from smodels.share.models.SMparticles import SMList
from smodels.share.models.mssm import BSMList
from smodels.matching.theoryPrediction import theoryPredictionsFor
from smodels.statistics.simplifiedLikelihoods import Data, UpperLimitComputer
from smodels.base.physicsUnits import fb
from smodels.decomposition import decomposer
from smodels.base.smodelsLogging import logger
from smodels.experiment.databaseObj import Database
from base.loggerbase import LoggerBase
from tester.combinationsmatrix import getYamlMatrix
from typing import Dict, List, Text
from icecream import ic

logger.setLevel("ERROR")

hasWarned = { "noupperlimits": 0 }

def readDictFile ( filename : str = "default.dict" ) -> Dict:
    """ read in content of filename
    :param filename: the filename of the database dictionary.
    often it is <dbversion>.dict

    :returns: a dictionary with 'meta' and 'data' as keys
    """

    with open( filename,"rt") as f:
        tmp=f.readlines()
    lines = []
    for line in tmp:
        if line.startswith("#"):
            continue
        lines.append ( line )
    basename = os.path.basename ( filename ).replace(".dict","")
    meta = eval(lines[0])
    nan=float("nan")
    inf=float("inf")
    data = eval("\n".join(lines[1:]))
    newdata = {}
    for i,v in data.items():
        if "expectedBG" in v and v["expectedBG"]>=0.:
            newdata[i]=v
        else:
            if i.endswith ( ":ul" ):
                pass
            else:
                eBG=None
                if "expectedBG" in v:
                    eBG = v["expectedBG"]
    return { "meta": meta, "data": newdata }

class ExpResModifier ( LoggerBase ):
    epilog="""
Examples:
=========

Fake SM-only database:
----------------------
./expResModifier.py -R $RUNDIR -d original.pcl -s fake1

Database with a fake signal:
----------------------------
./expResModifier.py -R $RUNDIR -d original.pcl -s signal1 -P pmodel9.dict

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

    def __init__ ( self, args ):
        """ args is a dictionary here,
        :param database: path to database
        :param max: upper limit on an individual excess
        :param suffix: suffix to use, e.g. fake, signal, etc
        :param lognormal: if True, use lognormal for nuisances, else Gaussian
        :param fixedsignals: if True, then use the central value of theory prediction
                             as the signal yield, dont draw from Poissonian
        :param fixedbackgrounds: if True, then use the central value of theory prediction
                             as the background yield, dont draw from Poissonian
        :param seed: if int and not None, set random number seed
        :param maxmassdist: maximum distance (in GeV) for the euclidean space in masses,
                            for a signal to populate an UL map
        :param compute_ps: compute p-values for all SRs
        """
        super ( ExpResModifier, self ).__init__ ( 0 )
        self.superseded = set() ## take note of everything superseded
        self.fastlim = set() # take note of everything fastlim
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
        if "seed" in args:
            self.setSeed ( args["seed"] )
        self.run()

    def defaults( self ):
        """ define the defaults """
        self.db = None
        self.comments = {} ## comments on entries in dict
        self.hasFiltered = False
        self.timestamps = False
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
        self.maxmassdist = 400
        self.seed = None
        self.nproc = 1
        self.pmodel = ""
        self.playback = ""
        self.verbose = 0
        self.interactive = False
        self.build = False
        self.check = False
        self.compute_ps = False
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
        picklefile = os.path.abspath ( self.rundir + "/" + self.dbpath )
        if self.rundir in self.dbpath:
            picklefile = self.dbpath
        self.pprint ( f"Extracting stats from {picklefile}" )
        if self.db == None:
            combinationsmatrix, status = getYamlMatrix()
            if not combinationsmatrix or status != 0:
                logger.error("Combination matrix not loaded correctly.")
            print ( f"[expResModifier] loading database {self.database}" )
            self.db = Database ( picklefile, combinationsmatrix=combinationsmatrix )
            print ( f"[expResModifier] loaded" )
        self.dbversion = self.db.databaseVersion
        listOfExpRes = self.db.expResultList
        self.stats = {}
        for expRes in listOfExpRes:
            for i,dataset in enumerate(expRes.datasets):
                dId = dataset.dataInfo.dataId
                if dId == None:
                    dId = "ul"
                label = dataset.globalInfo.id + ":" + dId
                D={}
                info = dataset.dataInfo
                dt = info.dataType
                if dt == "upperLimit":
                    for txname in dataset.txnameList:
                        D[txname.txName]=list ( txname.txnameData.y_values )

                for i in [ "observedN", "origN", "expectedBG", "lmbda", "bgError",
                           "origUpperLimit", "origExpectedUpperLimit", "upperLimit",
                           "expectedUpperLimit", "thirdMoment" ]:
                    if hasattr ( info, i ):
                        D[i] = getattr ( info, i )
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
            D["x"] = x
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
                    D["y0old"] = observed.y_values[i]
                    D["y0exp"] = y
                    D["sigma_exp0"] = sigma_exp
                    D["y0new"] = obs
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

        label = globalInfo.id + ":ul:" + txname.txName
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
            args += i + " "
        f.write ( f"[expResModifier.py-{time.strftime('%H:%M:%S')}] {args.strip()}\n")
        f.close()

    def startLogger ( self ):
        subprocess.getoutput ( f"mv {self.logfile} modifier.old" )
        self.log ( f"starting at {time.asctime()} with zmax of {self.max}" )
        self.log ( f"arguments were {' '.join ( sys.argv )}" )

    def finalize ( self ):
        """ finalize, for the moment its just deleting slha files """
        # print ( "[expResModifier] finalize" )
        if self.keep:
            return
        if hasattr ( self, "protomodel" ) and self.protomodel is not None and \
                type(self.protomodel) != str:
            self.protomodel.delCurrentSLHA()

    def produceProtoModel ( self, filename, dbversion ):
        """ try to produce a protomodel from pmodel
        :param filename: filename of pmodel dictionary
        :param dbversion: version of database, for tracking
        :returns: none if not succesful, else protomodel object
        """
        if filename == "":
            return None
        if not os.path.exists ( filename ):
            self.pprint ( f"When trying to construct protomodel, {filename} does not exist" )
            return None
        shutil.copyfile ( filename, self.rundir+"/my.signal" )
        walkerid = 0
        expected = False
        select = "all"
        keep_meta = True
        # M = ProtoModel ( walkerid, self.dbpath, expected, select, keep_meta )
        M = ProtoModel ( walkerid, keep_meta, dbversion = dbversion )
        M.createNewSLHAFileName ( prefix="erm" )
        ma = Manipulator ( M )
        with open ( filename, "rt" ) as f:
            m = eval ( f.read() )
        ma.initFromDict ( m, initTestStats=True )
        ma.M.computeXSecs( keep_slha = True )
        self.log ( f"xsecs produced {ma.M.currentSLHA}" )
        self.log ( f" `- does currentslha exist? {os.path.exists ( ma.M.currentSLHA )}" )
        ma.printXSecs()
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
            print ( f"[expResModifier] loading database {self.dbpath} [0]" )
            self.db = Database ( self.dbpath, combinationsmatrix=combinationsmatrix)
            print ( f"[expResModifier] loaded" )
        self.dbversion = self.db.databaseVersion
        listOfExpRes = self.removeEmpty ( self.db.expResultList ) ## seems to be the safest bet?
        self.produceProtoModel ( self.pmodel, self.db.databaseVersion )
        # print ( "pm produced", os.path.exists ( self.protomodel.currentSLHA ) )
        self.log ( f"{len(listOfExpRes)} results before faking bgs" )
        updatedListOfExpRes = self.fakeBackgrounds ( listOfExpRes )
        # print ( "fb produced", os.path.exists ( self.protomodel.currentSLHA ) )
        self.log ( f"{len(updatedListOfExpRes)} results after faking bgs" )
        updatedListOfExpRes = self.addSignals ( updatedListOfExpRes )
        self.log ( f"{len(updatedListOfExpRes)} results after adding signals" )
        if hasattr ( self.db, "subs" ): ## for smodels 2.1
            self.db.subs[0].expResultList = updatedListOfExpRes
            self.db.subs = [ self.db.subs[0] ]
        else:
            self.db.expResultList = updatedListOfExpRes
        newver = self.db.databaseVersion
        if self.suffix is not None:
            newver = self.db.databaseVersion + self.suffix
        self.db.txt_meta.databaseVersion = newver
        self.db.pcl_meta.databaseVersion = newver
        self.pprint ( f"Constructed fake database with {len(updatedListOfExpRes)} (of {len(listOfExpRes)}) results" )
        self.createBinaryFile()
        return self.db

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
        D = { "origN": int(orig), "expectedBG": exp, "bgError": err, "fudge": self.fudge,
              "lumi": float(dataset.globalInfo.lumi * fb) }
        if self.compute_ps:
            if thirdMoment is None:
                p = computeP ( orig, exp, err )
            else:
                D["thirdMoment"]=thirdMoment
                self.comments["thirdMoment"]="third moment for SLv2 likelihoods"
                p = computePSLv2 ( orig, exp, err, thirdMoment )
            self.comments["orig_p"]="p-value (Gaussian nuisance) of original observation"
            D["orig_p"]=p
            origZ = computeZFromP ( p )
            D["orig_Z"]=origZ
            self.comments["orig_Z"]="the significance Z of the original observation"
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
                if thirdMoment is None:
                    pnew = computeP ( orig, exp, err )
                else:
                    D["thirdMoment"]=thirdMoment
                    self.comments["thirdMoment"]="third moment for SLv2 likelihoods"
                    pnew = computePSLv2 ( orig, exp, err, thirdMoment )
                Z = - scipy.stats.norm.ppf ( pnew )
                # Z = ( obs - exp ) / toterr
                # origZ = ( orig - exp ) / toterr
            if Z < self.max or ct > 9:
                self.log ( f"effmap replacing old nobs={orig} (bg={exp:.2f}+/-{err:.2f}, lmbda={lmbda:.2f}, Z={Z:.2f}) with nobs={obs} for {dataset.globalInfo.id}:{dataset.dataInfo.dataId}" )
                dataset.dataInfo.observedN = obs
        if Z > 3.5:
            self.log ( f"WARNING!!! high em Z={Z:.2f}!!!!" )
        D["Zbg"]=Z
        self.comments["Zbg"]="the significance of the observation, bg only"
        D["Z"]=Z
        self.comments["Z"]="the significance of the observation, taking into account the signal"
        self.comments["lmbda"]="Poissonian lambda of the fake background"
        D["lmbda"]=lmbda
        D["newObs"]=obs
        self.comments["newObs"]="the new fake observation"
        if self.compute_ps:
            if thirdMoment is None:
                p = computeP ( obs, exp, err )
            else:
                D["thirdMoment"]=thirdMoment
                self.comments["thirdMoment"]="third moment for SLv2 likelihoods"
                p = computePSLv2 ( obs, exp, err, thirdMoment )
            self.comments["new_p"]="p-value (Gaussian nuisance) of newObs"
            D["new_p"]=p
        D["obsBg"]=obs
        self.comments["obsBg"]="the new fake observation, background component"
        D["toterr"]=toterr
        ## origN stores the n_observed of the original database
        dataset.dataInfo.origN = orig
        label = dataset.globalInfo.id + ":" + dataset.dataInfo.dataId
        txnames = [ tx.txName for tx in dataset.txnameList ]
        txnames.sort()
        if len ( txnames ) == 0:
            self.warning ( f"no txnames for {label}." )
        D["txns"]=",".join(txnames )
        self.comments["txns"]="list of txnames that populate this signal region / analysis"
        if self.timestamps:
            D["timestamp"]=dataset.globalInfo.lastUpdate
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
        if self.compute_ps:
            if thirdMoment is None:
                p = computeP ( orig, exp, err )
            else:
                D["thirdMoment"]=thirdMoment
                self.comments["thirdMoment"]="third moment for SLv2 likelihoods"
                p = computePSLv2 ( orig, exp, err, thirdMoment )
            self.comments["orig_p"]="p-value (Gaussian nuisance) of original observation"
            D["orig_p"]=p
            origZ = computeZFromP ( p )
            D["orig_Z"]=origZ
            self.comments["orig_Z"]="the significance Z of the original observation"
        txnames = [ tx.txName for tx in dataset.txnameList ]
        txnames.sort()
        if len ( txnames ) == 0:

            self.warning ( f"no txnames for {label}." )
        D["txns"]=",".join(txnames )
        self.comments["txns"]="list of txnames that populate this signal region / analysis"
        if self.timestamps:
            D["timestamp"]=dataset.globalInfo.lastUpdate
        return D

    def addSignalForEfficiencyMap ( self, dataset, tpred, lumi ):
        """ add a signal to this efficiency map. background sampling is
            already taken care of """
        txns = list ( map ( str, tpred.txnames ) )
        txns.sort()
        self.log ( f"add EM matching tpred {tpred.analysisId()}/{tpred.dataId()} {','.join(txns)}: {tpred.xsection.value}" )
        label = dataset.globalInfo.id + ":" + dataset.dataInfo.dataId
        orig = dataset.dataInfo.observedN
        sigLambda = float ( tpred.xsection.value * lumi )
        D={}
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
        self.comments["sigN"]="the number of events from the added signal"
        txnsc = "_".join( txns )
        ## sigNT<x> denotes the contributions from the individual theory preds
        D[ f"sigN{txnsc}" ] = sigN
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
        ## the signal is less than permille of bg?
        if orig > 0. and sigN / orig < 1e-3:
                self.log ( f" `- signal sigN={sigN} re obsN={orig} too small. skip.")
                dataset.dataInfo.origUpperLimit = dataset.dataInfo.upperLimit
                dataset.dataInfo.origExpectedUpperLimit = dataset.dataInfo.expectedUpperLimit
                D["newObs"]=orig
                self.addToStats ( label, D, dataset.globalInfo )
                return dataset
        self.log ( f" `- effmap adding sigN={sigN} to obsN={orig} -> newObs={orig+sigN}" )
        dataset.dataInfo.trueBG = orig ## keep track of true bg
        dataset.dataInfo.observedN = orig + sigN
        D["newObs"]=dataset.dataInfo.observedN
        exp = dataset.dataInfo.expectedBG
        err = dataset.dataInfo.bgError * self.fudge
        toterr = math.sqrt ( err**2 + exp )
        Z = 0.
        if toterr > 0.:
            Z = ( dataset.dataInfo.observedN - exp ) / toterr
        D["Z"]=Z
        self.comments["Z"]="the significance of the observation, taking into account the signal"
        ## now recompute the limits!!
        alpha = .05
        if orig == 0.0:
            orig = 0.00001
        computer = UpperLimitComputer(cl=1.-alpha )
        m = Data( orig+sigN, orig, err**2, nsignal = 1. )
        lumi = dataset.globalInfo.lumi# .asNumber(1./fb)
        maxSignalXsec = computer.ulSigma(m, marginalize=True ) / lumi
        dataset.dataInfo.origUpperLimit = dataset.dataInfo.upperLimit
        dataset.dataInfo.origExpectedUpperLimit = dataset.dataInfo.expectedUpperLimit
        dataset.dataInfo.upperLimit = maxSignalXsec
        maxSignalXsec = computer.ulSigma(m, marginalize=True, expected=True ) / lumi
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
        sigmaN = values["sigmaN"] # tpred.xsection.value.asNumber(fb)
        label = dataset.globalInfo.id + ":ul:" + txns
        D={}
        D["sigmaN"]=sigmaN
        D["pids"]=values["pids"]
        D["masses"]=values["masses"]
        D["txns"]=txns
        self.comments["txns"]="list of txnames that populate this signal region / analysis"
        self.comments["sigmaN"]="the added theory prediction (in fb), for UL maps"
        ## sigmaN is the predicted production cross section of the signal,
        ## in fb
        if not txname.txName in txns:
            return txname.txnameData
        hasAdded = 0
        txnd = txname.txnameData
        etxnd = txname.txnameDataExp
        coordsTpred = txnd.dataToCoordinates ( values["masses"], txnd._V, txnd.delta_x ) ## coordinates of tpred
        minDist = float("inf") ## for the closest point we store the numbers
        for yi,y in enumerate(txnd.y_values):
            pt = txnd.tri.points[yi] ## the point in the rotated coords
            dist = self.distance ( pt, coordsTpred )
            if dist > self.maxmassdist: ## change y_values only in vicinity of protomodel
                continue
            oldv = txnd.y_values[yi]
            oldo = txnd.y_values[yi]
            hasExpected=False
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
                D["ynew"]=oldv+sigmaN
            txnd.y_values[yi]=oldv + sigmaN
            hasAdded += 1
            if hasAdded == 0:
                self.pprint ( "warning: signal was not added in {tpred.analysisId()}:{txname.txName}" )
            D[f"signalpoints{txname.txName}"]=hasAdded
            D[f"totalpoints{txname.txName}"]=len(txnd.y_values)
            self.comments["signalpointsTx"]="number of grid points that got the signal injected"
            self.comments["totalpointsTx"]="total number of grid points in that map"
            self.addToStats ( label, D, dataset.globalInfo )
        return txnd

    def addSignalForULMap ( self, dataset, tpred, lumi ):
        """ add a signal to this UL result. background sampling is
            already taken care of """
        from smodels.base.physicsUnits import fb
        from ptools import helpers
        txns = list ( map ( str, tpred.txnames ) )
        txns.sort()
        self.log ( f"add UL matching tpred {tpred.analysisId()}: <{tpred.xsection.value}> {tpred.PIDs} {','.join(txns)}" )
        ## so we simply add the theory predicted cross section to the limit
        sigmaN = tpred.xsection.value.asNumber(fb)
        label = tpred.analysisId() + ":ul:" + ",".join(txns)
        D={}
        D["sigmaN"]=sigmaN
        D["pids"]=tpred.PIDs
        D["masses"]=helpers.stripUnits ( tpred.mass )
        D["txns"]=",".join(txns)
        self.comments["txns"]="list of txnames that populate this signal region / analysis"
        self.comments["sigmaN"]="the added theory prediction (in fb), for UL maps"
        ## sigmaN is the predicted production cross section of the signal,
        ## in fb
        for i,txname in enumerate(dataset.txnameList):
            if not self.txNameIsIn ( txname, tpred ):
                continue
            hasAdded = 0
            txnd = txname.txnameData
            etxnd = txname.txnameDataExp
            coordsTpred = txnd.dataToCoordinates ( tpred.mass, txnd._V, txnd.delta_x ) ## coordinates of tpred
            minDist, minPt = float("inf"),None ## for the closest point we store the numbers
            for yi,y in enumerate(txnd.y_values):
                pt = txnd.tri.points[yi] ## the point in the rotated coords
                dist = self.distance ( pt, coordsTpred )
                if dist < minDist: ## just so we know how far away we are
                    minDist = dist
                    minPt = txnd.coordinatesToData ( pt, txnd._V, txnd.delta_x )
                if dist > self.maxmassdist: ## change y_values only in vicinity of protomodel
                    continue
                oldv = txnd.y_values[yi]
                oldo = txnd.y_values[yi]
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
                    D["ynew"]=oldv+sigmaN
                txnd.y_values[yi]=oldv + sigmaN
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

    def saveStats ( self, statsname = None ):
        """ write out the collected stats, so we can discuss experimentalists'
            conservativeness """
        if self.suffix in [ None, "None", "", "none" ]:
            return
        filename = f"{self.rundir}/{self.suffix}.dict".replace("//","/")
        if statsname != None:
            filename = statsname
        self.log ( f"saving stats to {filename}" )
        meta = { "dbpath": self.dbpath, "Zmax": self.max,
                 "database": self.dbversion, "fudge": self.fudge,
                 "protomodel": f'"{self.protomodel}"', "timestamp": time.asctime(),
                 "lognormal": self.lognormal, "maxmassdist": self.maxmassdist,
                 "fixedsignals": self.fixedsignals,
                 "fixedbackgrounds": self.fixedbackgrounds }
        if hasattr ( runtime, "_drmax" ):
            meta["_drmax"]=runtime._drmax
        if hasattr ( runtime, "_experimental" ):
            meta["_experimental"]=runtime._experimental
        self.pprint ( f"saving stats to {filename}" )
        self.addSupersededFlags()
        with open ( filename,"wt" ) as f:
            f.write ( str(meta)+"\n" )
            if len(self.comments)>0:
                f.write ( "# explanations on the used variables:\n" )
                f.write ( "# =====================================\n" )
            else:
                f.write ( "# no explanations for variables have been given\n" )
            for k,v in self.comments.items():
                f.write ( f"# {k}: {v}\n" )
            f.write ( '{' )
            for ctr,(k,v) in enumerate(self.stats.items()):
                f.write ( f"'{k}': {v}" )
                if ctr != len(self.stats)-1:
                    f.write ( ",\n" )
            f.write ( '}\n' )
            f.close()

    def produceTopoList ( self ):
        """ create smstopolist """
        from smodels.base.physicsUnits import fb, GeV
        model = Model ( BSMList, SMList )
        model.updateParticles ( inputFile=self.protomodel.currentSLHA )
        mingap=10*GeV
        sigmacut = 0.02*fb
        self.topos = decomposer.decompose ( model, sigmacut, minmassgap=mingap )

    def addSignalsSingleProc ( self, listOfExpRes ):
        """ thats the method that adds a typical signal """
        if self.protomodel == None:
            return listOfExpRes
        ret = []
        self.produceTopoList()
        ctr,els = 0, ""
        for topo in self.topos:
            for el in topo.elementList:
                ctr+=1
                els += str(el) + ", "
        els=els[:-2]
        self.log ( f"now add the signals from {self.getPModelName()}, {ctr} topologies: {els}" )
        addedUL, addedEM = 0, 0
        self.pprint ( f"{len(listOfExpRes)} results: ", end="" )
        for l,expRes in enumerate(listOfExpRes):
            print ( ".", flush=True, end="" )
            tpreds = theoryPredictionsFor ( expRes, self.topos, useBestDataset=False,
                                            combinedResults=False )
            if tpreds == None:
                ret.append ( expRes )
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
        from smodels.statistics.simplifiedLikelihoods import Data
        data = Data ( observed, expectedBGs, covm, thirdMoments )
        # print ( f"@@3 for {expRes.globalInfo.id} rvs {rvs[:3]} cnt {centers[:3]} diag {diag[:3]}" )
        for i,dataset in enumerate(expRes.datasets):
            lmbda = dataset.dataInfo.expectedBG + rvs[i]
            if tpe == "SLv2":
                lmbda = data.A[i] + rvs[i] + data.C[i] * rvs[i]**2 / data.B[i]**2
                # lmbda += dataset.dataInfo.thirdMoments / diag[i]**2 * rvs[i]**2
            lmbda = max ( 0., lmbda )
            obs = scipy.stats.poisson.rvs ( lmbda )
            D = self.createEMStatsDict ( dataset )
            if self.fixedbackgrounds:
                D["newObs"]=dataset.dataInfo.expectedBG
            else:
                D["newObs"]=obs
                D["lmbda"]=lmbda
            if self.compute_ps:
                p = computePForDataSet ( dataset, obs )
                self.comments["new_p"]="p-value (Gaussian nuisance) of newObs"
                D["new_p"]=p
                newZ = computeZFromP ( p )
                D["new_Z"]=newZ
            D["type"]=tpe
            self.comments["type"]="result type (None, SLv1, SLv2, pyhf)"
            expRes.datasets[i].dataInfo.observedN = obs
            label = dataset.globalInfo.id + ":" + dataset.dataInfo.dataId
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
                                # print ( f"@@9 old hi_data {hi,lo}" ) 
                                center = (hi+lo)/2.
                                delta = hi - center
                                data["hi_data"][idd] = center + delta * self.fudge
                                data["lo_data"][idd] = center - delta * self.fudge
                                # print ( f"@@9 new hi {computer.likelihoodComputer.workspaces[iws]['channels'][ich]['samples'][ism]['modifiers'][im]['data']['hi_data'][idd]}" )

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
        from smodels.base.runtime import _deltas_rel_default
        from smodels.statistics.statsTools import StatsComputer
        from smodels.experiment.datasetObj import CombinedDataSet
        import pyhf
        cdataset = CombinedDataSet ( expRes )
        computer = StatsComputer.forPyhf( cdataset, srNsigDict,
                _deltas_rel_default )
        if abs ( self.fudge - 1. ) > 1e-5:
            self.fudgePyhfModel ( expRes, computer )
        srs_in_workspaces = list(expRes.globalInfo.jsonFiles.values())
        anaId = expRes.globalInfo.id
        #ic ( anaId )
        # import sys, IPython; IPython.embed( colors = "neutral" ); sys.exit()
        for ws_i, (ws, srs) in enumerate(zip(
                    computer.likelihoodComputer.workspaces, srs_in_workspaces) ):
            ## srs are the names of the signal regions
            model = ws.model()
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
            #if anaId == "ATLAS-SUSY-2018-31":
            #    ic ( srs )
            #    ic ( sampleDictPyhf )
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
                D["type"]="pyhf"
                if self.compute_ps:
                    p = computePForDataSet ( dataset, newObs )
                    self.comments["new_p"]="p-value (Gaussian nuisance) of newObs"
                    D["new_p"]=p
                    newZ = computeZFromP ( p )
                    D["new_Z"]=newZ
                label = anaId + ":" + dataset.dataInfo.dataId
                self.addToStats ( label, D, dataset.globalInfo )
                ## as the very last measure, we replace the observation with
                ## the fake observation
                dataset.dataInfo.observedN = newObs
        #if anaId == "ATLAS-SUSY-2018-31":
        #    import sys, IPython; IPython.embed( colors = "neutral" ); sys.exit()


    def fakeBackgrounds ( self, listOfExpRes ):
        """ thats the method that samples the backgrounds """
        ret = []
        self.log ( "now fake backgrounds" )
        for expRes in listOfExpRes:
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
        self.log ( "done faking the backgrounds" )
        return ret

    def cleanTxNameData ( self, txnd ):
        txnd.y_values=numpy.array ( txnd.y_values, dtype=numpy.float32 )

        if txnd.dimensionality == 1:
            return txnd
        txnd.tri._points = numpy.array ( txnd.tri._points, dtype=numpy.float32 )
        return txnd

    def filter ( self ):
        """ filter the list fo experimental results.
        :param outfile: store result in outfile (a pickle file)
        :param nofastlim: remove fastlim results
        :param onlyvalidated: remove non-validated results
        :param nosuperseded: remove superseded results
        :param remove_orig: remove original values
        :param remove_nonagg: remove non-aggregated results
        """
        if not ( self.nofastlim or self.onlyvalidated or self.nosuperseded or self.remove_orig or self.remove_nonagg or self.noupperlimits ):
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
        self.dbversion = self.dbversion + ".playedback"
        self.db.txt_meta.databaseVersion = self.db.databaseVersion + ".playedback"
        self.db.pcl_meta.databaseVersion = self.db.databaseVersion + ".playedback"
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
            if not filecmp.cmp ( f, self.rundir+"/"+os.path.basename ( f ) ):
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
            print ( "[expResModifier] WARNING fixing backgrounds but not signals. Sounds weird" )
        if self.fixedbackgrounds and self.fudge > 1e-2:
            print ( "[expResModifier] WARNING fixing backgrounds but fudge factor is not zero. Sounds weird" )
        if self.build:
            from smodels.experiment.txnameObj import TxNameData
            TxNameData._keep_values = True
            from smodels.experiment.databaseObj import Database
            print ( f"[expResModifier] starting to build database at {self.database}." )
            combinationsmatrix, status = getYamlMatrix()
            if not combinationsmatrix or status != 0:
                logger.error("Combination matrix not loaded correctly.")
            print ( f"[expResModifier] loading database {self.database}" )
            db = Database ( self.database, combinationsmatrix=combinationsmatrix )
            print ( f"[expResModifier] built database at {self.database}. Exiting." )
            sys.exit()
        if self.rundir == None:
            print ( f"[expResModifier] setting rundir to {os.getcwd()}" )
            self.rundir = os.getcwd()
        if type(self.rundir)==str and not "/" in self.rundir and \
                not self.rundir.startswith("."):
            self.rundir = f"{os.environ['HOME']}/{self.rundir}"
        statsname = self.suffix + ".dict"
        if self.outfile is not None:
            if self.outfile == "":
                self.outfile = self.suffix+".pcl"
            if self.playback not in [ None, "" ]:
                self.playback ( self.playback, self.outfile )
                statsname = "playback.dict"

            if not self.outfile.endswith(".pcl") and self.outfile != None:
                print ( f"[expResModifier] warning, shouldnt the name of your outputfile ``{self.outfile}'' end with .pcl?" )
        #else: # outfile is None
        #    statsname = None
        self.filter ( )
        if self.dontsample:
            print ( "[expResModifier] we were asked to not sample, so we exit now." )
            sys.exit()

        if self.extract_stats:
            er = self.extractStats()
        else:
            if not self.playback:
                er = self.modifyDatabase ( )
        if statsname is not None:
            self.saveStats( statsname )

        if self.check:
            self.check ( )

        if self.interactive:
            self.interact ( er )

        if self.upload:
            self.upload()

        if self.symlink:
            self.symlink ( )

        self.finalize()

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
    argparser.add_argument ( '--maxmassdist',
            help='maximum euclidean distance in mass space to add the signal in the UL maps [400.]',
            type=float, default=400. )
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
    argparser.add_argument ( '-C', '--compute_ps',
            help='compute p-values for all SRs', action='store_true' )
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
    modifier = ExpResModifier( args.__dict__ )
