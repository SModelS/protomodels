#!/usr/bin/env python3

""" a simple class that performs global modifications on a list of results.
Used to ``take out potential signals'' i.e. put all observations to values
expected from background, by sampling the background model. """

# https://link.springer.com/content/pdf/10.1007/JHEP02(2015)004.pdf

import copy, os, sys, time, subprocess, math, numpy, shutil
import scipy.spatial
sys.path.insert( 0, "../" )
sys.path.append('../smodels')
from csetup import setup
setup()
#sys.path.insert(0,"/scratch-cbe/users/wolfgan.waltenberger/git/protomodels/")
from scipy import stats
from builder.protomodel import ProtoModel
from builder.manipulator import Manipulator
from helpers import computeP
from smodels.theory.model import Model
from smodels.share.models.SMparticles import SMList
from smodels.particlesLoader import BSMList
from smodels.theory.theoryPrediction import theoryPredictionsFor
from smodels.tools.simplifiedLikelihoods import Data, UpperLimitComputer
from smodels.tools.physicsUnits import fb
from smodels.theory import decomposer
from smodels.tools.smodelsLogging import logger
from smodels.experiment.databaseObj import Database

logger.setLevel("ERROR")

class ExpResModifier:
    def __init__ ( self, dbpath, Zmax, rundir, keep, nproc, fudge,
                   suffix: str, lognormal = False, fixedsignals = False,
                   fixedbackgrounds = False, seed = None,
                   maxmassdist = 400., compute_ps = False ):
        """
        :param dbpath: path to database
        :param Zmax: upper limit on an individual excess
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
        self.db = None
        self.comments = {} ## comments on entries in dict
        self.lognormal = lognormal
        self.dbpath = dbpath
        self.maxmassdist = maxmassdist
        self.compute_ps = compute_ps
        self.hasFiltered = False
        self.fixedsignals = fixedsignals
        self.fixedbackgrounds = fixedbackgrounds
        self.protomodel = None
        self.rundir = setup( rundir )
        self.keep = keep
        self.nproc = nproc
        self.fudge = fudge
        self.logfile = "modifier.log"
        self.suffix = suffix
        if Zmax == None:
            Zmax = 100
        self.Zmax = Zmax
        self.startLogger()
        self.stats = {}
        self.logCall()
        self.setSeed ( seed )

    def setSeed( self, seed ):
        if seed is None:
            return
        from ptools import helpers
        helpers.seedRandomNumbers( seed )
        self.pprint ( f"setting random seed to {args.seed}" )

    def interact ( self, listOfExpRes ):
        import IPython
        IPython.embed( using=False )

    def extractStats ( self ):
        """ dont produce a new fake database, extract a stats dict
            from an existing database. """
        self.info ( "extracting stats" )
        picklefile = self.rundir + "/" + self.dbpath
        if self.rundir in self.dbpath:
            picklefile = self.dbpath
        self.pprint ( f"Extracting stats from {picklefile}" )
        if self.db == None:
            self.db = Database ( picklefile )
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
                        D[txname]=list ( txname.txnameData.y_values )

                for i in [ "observedN", "origN", "expectedBG", "lmbda", "bgError",
                           "origUpperLimit", "origExpectedUpperLimit", "upperLimit",
                           "expectedUpperLimit" ]:
                    if hasattr ( info, i ):
                        D[i] = getattr ( info, i )
                self.addToStats ( label, D )

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
        while not allpositive:
            ret = copy.deepcopy ( expected )
            ctr += 1
            x = float("inf")
            if x_ != None:
                x = x_
            D = {}
            while x > self.Zmax:
                x = 0.
                if not self.fixedbackgrounds:
                    x = self.drawNuisance() * self.fudge # draw but once from standard-normal
                # x = stats.norm.rvs() * self.fudge # draw but once from standard-normal
            D["x"] = x
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
                self.log ( "WARNING seems like I am having a hard time getting all "\
                        "values of %s positive." % globalInfo.id )

        label = globalInfo.id + ":ul:" + txname.txName
        D["fudge"]=self.fudge
        self.addToStats ( label, D )
        self.log ( "computed new UL result %s:%s, x=%.2f" % \
                   ( globalInfo.id, txname.txName, x ) )
        if x > 3.5:
            self.log ( "WARNING high UL x=%.2f!!!" % x )
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
                i = '"%s"' % i
            args += i + " "
        f.write ( "[expResModifier.py-%s] %s\n" % \
                  ( time.strftime("%H:%M:%S"), args.strip() ) )
        # f.write ("[slurm.py] %s\n" % " ".join ( sys.argv ) )
        f.close()

    def pprint ( self, *args ):
        """ logging """
        print ( "[expResModifier] %s" % ( " ".join(map(str,args))) )
        with open( self.logfile, "a" ) as f:
            f.write ( "[modifier] %s\n" % ( " ".join(map(str,args)) ) )

    def startLogger ( self ):
        subprocess.getoutput ( "mv %s modifier.old" % self.logfile )
        self.log ( "starting at %s with zmax of %s" % \
                   ( time.asctime(), self.Zmax ) )
        self.log ( "arguments were %s" % ( " ".join ( sys.argv ) ) )

    def log ( self, *args ):
        """ logging to file """
        # logfile = "walker%d.log" % self.walkerid
        with open( self.logfile, "a" ) as f:
            f.write ( "[modifier] %s\n" % ( " ".join(map(str,args)) ) )

    def info ( self, *args ):
        """ logging to file, but also write to screen """
        self.log ( *args )
        print ( "[modifier] %s" % ( " ".join(map(str,args)) ) )

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
            self.pprint ( "When trying to construct protomodel, %s does not exist" % filename )
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
        self.log ( " `- does currentslha exist? %s" % \
                   os.path.exists ( ma.M.currentSLHA ) )
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

    def modifyDatabase ( self, outfile="", pmodel="" ):
        """ modify the database, possibly write out to a pickle file
        :param outfile: if not empty, write the database into file
        :param suffix: suffix to append to database version
        :param pmodel: if not empty, then this is the file name of the signal
                       model. in this case fake a signal
        :returns: the database
        """
        spmodel = f"protomodel is '{pmodel}'"
        if pmodel == "":
            spmodel = "no protomodel given"
        self.info ( f"starting to create {outfile} from {self.dbpath}. suffix is '{self.suffix}', {spmodel}." )
        if self.db == None:
            self.db = Database ( self.dbpath )
        self.dbversion = self.db.databaseVersion
        listOfExpRes = self.removeEmpty ( self.db.expResultList ) ## seems to be the safest bet?
        self.produceProtoModel ( pmodel, self.db.databaseVersion )
        # print ( "pm produced", os.path.exists ( self.protomodel.currentSLHA ) )
        self.log ( "%d results before faking bgs" % len(listOfExpRes) )
        updatedListOfExpRes = self.fakeBackgrounds ( listOfExpRes )
        # print ( "fb produced", os.path.exists ( self.protomodel.currentSLHA ) )
        self.log ( "%d results after faking bgs" % len(updatedListOfExpRes) )
        updatedListOfExpRes = self.addSignals ( updatedListOfExpRes )
        self.log ( "%d results after adding signals" % len(updatedListOfExpRes) )
        if hasattr ( self.db, "subs" ): ## for smodels 2.1
            self.db.subs[0].expResultList = updatedListOfExpRes
            self.db.subs = [ self.db.subs[0] ]
        else:
            self.db.expResultList = updatedListOfExpRes
        newver = self.db.databaseVersion + self.suffix
        self.db.txt_meta.databaseVersion = newver
        self.db.pcl_meta.databaseVersion = newver
        self.pprint ( "Constructed fake database with %d (of %d) results" % \
                ( len(updatedListOfExpRes), len(listOfExpRes) ) )
        if outfile != "":
            self.db.createBinaryFile( outfile )
        return self.db

    def sampleEfficiencyMap ( self, dataset ):
        """ for the given dataset,
        sample from background and put the value as observed """
        orig = dataset.dataInfo.observedN
        exp = dataset.dataInfo.expectedBG
        err = 0.
        if not self.fixedbackgrounds:
            err = dataset.dataInfo.bgError * self.fudge
        D = { "origN": orig, "expectedBG": exp, "bgError": err, "fudge": self.fudge,
              "lumi": float(dataset.globalInfo.lumi * fb) }
        if self.compute_ps:
            p = computeP ( orig, exp, err )
            self.comments["orig_p"]="p-value (Gaussian nuisance) of original observation"
            D["orig_p"]=p
        S, origS = float("inf"), float("nan")
        while S > self.Zmax:
            # lmbda = stats.norm.rvs ( exp, err )
            lmbda = exp
            if not self.fixedbackgrounds:
                lmbda = self.drawNuisance ( exp, err )
            dataset.dataInfo.lmbda = lmbda
            if lmbda < 0.:
                lmbda = 0.
            obs = lmbda
            toterr = 0.
            if not self.fixedbackgrounds:
                obs = stats.poisson.rvs ( lmbda )
                toterr = math.sqrt ( err**2 + exp )
            S, origS = 0., 0.
            if toterr > 0.:
                S = ( obs - exp ) / toterr
                origS = ( orig - exp ) / toterr
            if S < self.Zmax:
                self.log ( "effmap replacing old nobs=%d (bg=%.2f+/-%.2f, lmbda=%.2f, S=%.2f) with nobs=%d for %s:%s" % \
                    ( orig, exp, err, lmbda, S, obs, dataset.globalInfo.id, dataset.dataInfo.dataId  ) )
                dataset.dataInfo.observedN = obs
        if S > 3.5:
            self.log ( "WARNING!!! high em S=%.2f!!!!" % S )
        D["Sbg"]=S
        self.comments["Sbg"]="the significance of the observation, bg only"
        D["S"]=S
        self.comments["S"]="the significance of the observation, taking into account the signal"
        D["origS"]=origS
        self.comments["origS"]="the significance of the original observation"
        self.comments["lmbda"]="Poissonian lambda of the fake background"
        D["lmbda"]=lmbda
        D["newObs"]=obs
        self.comments["newObs"]="the new fake observation"
        if self.compute_ps:
            p = computeP ( obs, exp, err )
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
            print ( f"[expResModifier] warning, no txnames for {label}." )
        D["txns"]=",".join(txnames )
        self.comments["txns"]="list of txnames that populate this signal region / analysis"
        self.addToStats ( label, D )
        return dataset

    def addSignalForEfficiencyMap ( self, dataset, tpred, lumi ):
        """ add a signal to this efficiency map. background sampling is
            already taken care of """
        txns = list ( map ( str, tpred.txnames ) )
        txns.sort()
        self.log ( "add EM matching tpred %s/%s %s: %s" % \
                ( tpred.analysisId(), tpred.dataId(), ",".join(txns), \
                  tpred.xsection.value ) )
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
        D["sigN%s" % txnsc ] = sigN
        D["obsBg"]=self.stats[label]["newObs"]
        err = dataset.dataInfo.bgError * self.fudge
        dataset.dataInfo.sigN = sigN ## keep track of signal
        if sigN == 0:
                self.log ( " `- signal sigN=%d re obsN=%d too small. skip." % \
                           ( sigN, orig ) )
                dataset.dataInfo.origUpperLimit = dataset.dataInfo.upperLimit
                dataset.dataInfo.origExpectedUpperLimit = dataset.dataInfo.expectedUpperLimit
                D["newObs"]=orig
                self.addToStats ( label, D )
                return dataset
        ## the signal is less than permille of bg?
        if orig > 0. and sigN / orig < 1e-3:
                self.log ( " `- signal sigN=%d re obsN=%d too small. skip." % \
                           ( sigN, orig ) )
                dataset.dataInfo.origUpperLimit = dataset.dataInfo.upperLimit
                dataset.dataInfo.origExpectedUpperLimit = dataset.dataInfo.expectedUpperLimit
                D["newObs"]=orig
                self.addToStats ( label, D )
                return dataset
        self.log ( " `- effmap adding sigN=%d to obsN=%d -> newObs=%d" % \
                   ( sigN, orig, orig + sigN ) )
        dataset.dataInfo.trueBG = orig ## keep track of true bg
        dataset.dataInfo.observedN = orig + sigN
        D["newObs"]=dataset.dataInfo.observedN
        exp = dataset.dataInfo.expectedBG
        err = dataset.dataInfo.bgError * self.fudge
        toterr = math.sqrt ( err**2 + exp )
        S = 0.
        if toterr > 0.:
            S = ( dataset.dataInfo.observedN - exp ) / toterr
        D["S"]=S
        self.comments["S"]="the significance of the observation, taking into account the signal"
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
        self.addToStats ( label, D )
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

    def addToStats ( self, label, Dict ):
        """ add the content of dictionary Dict to the stats,
            under the label "label" """
        if not label in self.stats:
            # we dont yet have an entry, so lets start
            self.stats[label]=Dict
            return
        # we have an entry, so we add
        for k,v in Dict.items():
            self.stats[label][k]=v

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
        from smodels.tools.physicsUnits import fb
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
        # print ( "  `-- adding %s to %s" % ( sigmaN, txname ), type(txname), type(txname.txnameData) )
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
            # print ( "    `--- adding %s %s" % ( oldv, sigmaN ) )
            txnd.y_values[yi]=oldv + sigmaN
            hasAdded += 1
            if hasAdded == 0:
                self.pprint ( "warning: signal was not added in {tpred.analysisId()}:{txname.txName}" )
            D[f"signalpoints{txname.txName}"]=hasAdded
            D[f"totalpoints{txname.txName}"]=len(txnd.y_values)
            self.comments["signalpointsTx"]="number of grid points that got the signal injected"
            self.comments["totalpointsTx"]="total number of grid points in that map"
            self.addToStats ( label, D )
        return txnd

    def addSignalForULMap ( self, dataset, tpred, lumi ):
        """ add a signal to this UL result. background sampling is
            already taken care of """
        from smodels.tools.physicsUnits import fb
        from ptools import helpers
        txns = list ( map ( str, tpred.txnames ) )
        txns.sort()
        self.log ( "add UL matching tpred %s: <%s> %s {%s}" % \
                ( tpred.analysisId(), tpred.xsection.value, \
                  tpred.PIDs, ",".join(txns) ) )
        #print ( " `- add UL matching tpred %s: %s[%s] ds:%s" % \
        #        ( tpred.analysisId(), tpred.xsection.value, \
        #          tpred.PIDs, dataset ) )
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
            #print ( "  `-- adding %s to %s" % ( sigmaN, txname ) )
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
                # print ( "    `--- adding %s %s" % ( oldv, sigmaN ) )
                txnd.y_values[yi]=oldv + sigmaN
                hasAdded += 1
            if hasAdded == 0:
                self.pprint ( f"warning: no signal was added in {tpred.analysisId()}:{txname.txName}, closest was point {minPt} at d={minDist:.2f}" )
            D[f"signalpoints{txname.txName}"]=hasAdded
            D[f"totalpoints{txname.txName}"]=len(txnd.y_values)
            self.comments["signalpointsTx"]="number of grid points that got the signal injected"
            self.comments["totalpointsTx"]="total number of grid points in that map"
            self.addToStats ( label, D )
            dataset.txnameList[i].txnameData = txnd
            dataset.txnameList[i].sigmaN = sigmaN
        return dataset

    def saveStats ( self, statsname = None ):
        """ write out the collected stats, so we can discuss experimentalists'
            conservativeness """
        filename = "%s/db%s.dict" % ( self.rundir, self.suffix )
        if statsname != None:
            filename = statsname
        self.log ( f"saving stats to {filename}" )
        meta = { "dbpath": self.dbpath, "Zmax": self.Zmax,
                 "database": self.dbversion, "fudge": self.fudge,
                 "protomodel": '"%s"' % self.protomodel, "timestamp": time.asctime(),
                 "lognormal": self.lognormal, "maxmassdist": self.maxmassdist,
                 "fixedsignals": self.fixedsignals,
                 "fixedbackgrounds": self.fixedbackgrounds }
        from smodels.tools import runtime
        if hasattr ( runtime, "_cap_likelihoods" ):
            meta["_cap_likelihoods"]=runtime._cap_likelihoods
        if hasattr ( runtime, "_drmax" ):
            meta["_drmax"]=runtime._drmax
        if hasattr ( runtime, "_experimental" ):
            meta["_experimental"]=runtime._experimental
        self.pprint ( f"saving stats to {filename}" )
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
        from smodels.tools.physicsUnits import fb, GeV
        model = Model ( BSMList, SMList )
        model.updateParticles ( inputFile=self.protomodel.currentSLHA )
        mingap=10*GeV
        sigmacut = 0.02*fb
        self.topos = decomposer.decompose ( model, sigmacut, minmassgap=mingap )

    def addSignalsSingleProc ( self, listOfExpRes ):
        """ thats the method that adds a typical signal """
        # print ( "adding signals", os.path.exists ( self.protomodel.currentSLHA ) )
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
        self.log ( "now add the signals from %s, %d topologies: %s" % \
                   ( self.getPModelName(), ctr, els ) )
        addedUL, addedEM = 0, 0
        print ( f"{len(listOfExpRes)} results: ", end="" )
        for l,expRes in enumerate(listOfExpRes):
            print ( ".", flush=True, end="" )
            tpreds = theoryPredictionsFor ( expRes, self.topos, useBestDataset=False,
                                            combinedResults=False )
            if tpreds == None:
                ret.append ( expRes )
                continue
            lumi = expRes.globalInfo.lumi
            #self.pprint ( "adding a signal for %s (lumi %s)" % \
            #              ( expRes.id(), lumi ) )
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
            #self.pprint ( "adding a signal for %s (lumi %s)" % \
            #              ( expRes.id(), lumi ) )
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
        self.log ( "now add the signals from %s, K=%s, %d topos, %d procs" % \
                   ( self.getPModelName(), len(self.topos), self.nproc ) )
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

    def fakeBackgrounds ( self, listOfExpRes ):
        """ thats the method that samples the backgrounds """
        ret = []
        self.log ( "now fake backgrounds" )
        for expRes in listOfExpRes:
            for i,dataset in enumerate(expRes.datasets):
                dt = dataset.dataInfo.dataType
                if dt == "upperLimit":
                    expRes.datasets[i] = self.bgUpperLimit ( dataset )
                elif dt == "efficiencyMap":
                    expRes.datasets[i] = self.sampleEfficiencyMap ( dataset )
                else:
                    print ( "[expResModifier] dataset type %s unknown" % dt )
            ret.append ( expRes )
        self.log ( "done faking the backgrounds" )
        return ret

    def cleanTxNameData ( self, txnd ):
        txnd.y_values=numpy.array ( txnd.y_values, dtype=numpy.float32 )

        if txnd.dimensionality == 1:
            return txnd
        txnd.tri._points = numpy.array ( txnd.tri._points, dtype=numpy.float32 )
        return txnd

    def filter ( self, outfile, nofastlim, onlyvalidated, nosuperseded,
                       remove_orig, remove_nonagg ):
        """ filter the list fo experimental results.
        :param outfile: store result in outfile (a pickle file)
        :param nofastlim: remove fastlim results
        :param onlyvalidated: remove non-validated results
        :param nosuperseded: remove superseded results
        :param remove_orig: remove original values
        :param remove_nonagg: remove non-aggregated results
        """
        self.log ( "starting to filter %s. suffix is %s." % \
                   ( outfile, self.suffix ) )
        if self.db == None:
            self.db = Database ( self.dbpath )
        listOfExpRes = self.db.expResultList ## seems to be the safest bet?
        if remove_nonagg:
            from smodels_utils.helper.databaseManipulations import filterNonAggregatedFromList
            n = len(listOfExpRes )
            listOfExpRes = filterNonAggregatedFromList ( listOfExpRes, verbose=True )
            print ( f"[modifier] nonaggregated filter: from {n} to {len(listOfExpRes)}" )
            if len(listOfExpRes ) <  n:
                self.hasFiltered = True
        newList = []
        for er in listOfExpRes:
            addThisOne = True
            if nofastlim:
                if hasattr ( er.globalInfo, "contact" ) and "fastlim" in er.globalInfo.contact:
                    print ( " `- skipping fastlim %s" % er.globalInfo.id )
                    addThisOne = False
                    self.hasFiltered = True
            if nosuperseded:
                if hasattr ( er.globalInfo, "supersededBy" ):
                    print ( " `- skipping superseded %s" % er.globalInfo.id )
                    addThisOne = False
                    self.hasFiltered = True
            if hasattr ( er.globalInfo, "private" ) and er.globalInfo.private in [ "True", True ]:
                    print ( " `- skipping private %s" % er.globalInfo.id )
                    addThisOne = False
                    self.hasFiltered = True
            if not addThisOne:
                self.hasFiltered = True
                continue
            if onlyvalidated:
                newDs = []
                hasIssued = 0
                for ds in er.datasets:
                    txnew = []
                    for txn in ds.txnameList:
                        if txn.validated == False:
                            if hasIssued == 0:
                                print ( " `- skipping non-validated %s/%s/%s" % \
                                    ( txn.txName, ds.dataInfo.dataId, er.globalInfo.id ) )
                            if hasIssued == 1:
                                print ( " `- (suppressed more, similar messages)" )
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
            if remove_orig:
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
        if outfile != "":
            self.db.createBinaryFile( outfile )

    def playback ( self, playbackdict, outfile ):
        """ playback the mods described in playbackdict """
        self.pprint ( "WARNING: playback functionality has not yet been validated!!" )
        with open ( playbackdict, "rt" ) as h:
            lines = h.readlines()
            h.close()
        ## first line goes directly into database
        line = lines.pop(0)
        D = eval ( line )
        if self.db == None:
            self.db = Database ( self.dbpath )
        for k,v in D.items():
            if k == "dbpath":
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
        self.db.dbpath = outfile
        self.dbversion = self.dbversion + ".playedback"
        self.db.txt_meta.databaseVersion = self.db.databaseVersion + ".playedback"
        self.db.pcl_meta.databaseVersion = self.db.databaseVersion + ".playedback"
        self.pprint ( f"writing to {outfile}" )
        self.db.createBinaryFile ( outfile )

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
        # cmd = "cp %s ./modifier.log %s" % ( args.outfile, self.rundir )
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
        cmd = "ln -s %s/%s %s/default.pcl" % ( self.rundir, args.outfile, self.rundir )
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

    def check ( self, picklefile ):
        """ check the picklefile """
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
        print ( "were good", self.db.databaseVersion )


epilog="""
Examples:
=========

Fake SM-only database:
----------------------
./expResModifier.py -R $RUNDIR -d original.pcl -s fake1

Database with a fake signal:
----------------------------
./expResModifier.py -R $RUNDIR -d original.pcl -s signal1 -P pmodel9.py

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

if __name__ == "__main__":
    import argparse
    from argparse import RawTextHelpFormatter
    argparser = argparse.ArgumentParser(
                        description='experimental results modifier. used to take out potential signals from the database by setting all observations to values sampled from the background expectations. can insert signals, too.', formatter_class = RawTextHelpFormatter,
                        epilog=epilog )
    argparser.add_argument ( '-d', '--database',
            help='database to use [../../smodels-database]',
            type=str, default="../../smodels-database" )
    argparser.add_argument ( '-o', '--outfile',
            help='file to write out database pickle. If left empty, then outfile is <suffix>.pcl [""]',
            type=str, default="" )
    argparser.add_argument ( '-s', '--suffix',
            help='suffix for database version ["fake1"]',
            type=str, default="fake1" )
    argparser.add_argument ( '-R', '--rundir',
            help='override rundir [None]',
            type=str, default=None )
    argparser.add_argument ( '-f', '--fudge',
            help='fudge factor [1.0]',
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
            help='use lognormal, not Gaussian for nuisances',
            action='store_true' )
    argparser.add_argument ( '--fixedsignals',
            help='fix the contributions from the signals, dont draw from Poissonian',
            action='store_true' )
    argparser.add_argument ( '--fixedbackgrounds',
            help='fix the contributions from the backgrounds, use central values for all.',
            action='store_true' )
    argparser.add_argument ( '-M', '--max',
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
            help='build the original pickle file with all relevant info, then exit (use --database to specify path)', action='store_true' )
    argparser.add_argument ( '-c', '--check',
            help='check the pickle file <outfile>', action='store_true' )
    argparser.add_argument ( '-C', '--compute_ps',
            help='compute p-values for all SRs', action='store_true' )
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
    if args.fixedbackgrounds and not args.fixedsignals:
        print ( "[expResModifier] WARNING fixing backgrounds but not signals. Sounds weird" )
    if args.fixedbackgrounds and args.fudge > 1e-2:
        print ( "[expResModifier] WARNING fixing backgrounds but fudge factor is not zero. Sounds weird" )
    if args.build:
        from smodels.experiment.txnameObj import TxNameData
        TxNameData._keep_values = True
        from smodels.experiment.databaseObj import Database
        print ( f"[expResModifier] starting to build database at {args.database}." )
        db = Database ( args.database )
        print ( f"[expResModifier] built database at {args.database}. Exiting." )
        sys.exit()
    if args.rundir == None:
        print ( f"[expResModifier] setting rundir to {os.getcwd()}" )
        args.rundir = os.getcwd()
    if type(args.rundir)==str and not "/" in args.rundir and \
            not args.rundir.startswith("."):
        args.rundir = "/scratch-cbe/users/wolfgan.waltenberger/" + args.rundir
    if args.outfile == "":
        args.outfile = args.suffix+".pcl"

    modifier = ExpResModifier( args.database, args.max, args.rundir, args.keep, \
                               args.nproc, args.fudge, args.suffix, args.lognormal,
                               args.fixedsignals, args.fixedbackgrounds, args.seed,
                               args.maxmassdist, args.compute_ps )

    statsname = None
    if args.playback not in [ None, "" ]:
        modifier.playback ( args.playback, args.outfile )
        statsname = "playback.dict"

    if not args.outfile.endswith(".pcl"):
        print ( "[expResModifier] warning, shouldnt the name of your outputfile ``%s'' end with .pcl?" % args.outfile )
    if args.nofastlim or args.onlyvalidated or args.nosuperseded or args.remove_orig or args.remove_nonagg:
        modifier.filter ( args.outfile, args.nofastlim, args.onlyvalidated,
                          args.nosuperseded, args.remove_orig, args.remove_nonagg )
    if args.dontsample:
        print ( "[expResModifier] we were asked to not sample, so we exit now." )
        sys.exit()

    if args.extract_stats:
        er = modifier.extractStats()
    else:
        if not args.playback:
            er = modifier.modifyDatabase ( args.outfile, args.pmodel )

    modifier.saveStats( statsname )

    if args.check:
        modifier.check ( args.outfile )

    if args.interactive:
        modifier.interact ( er )

    if args.upload:
        modifier.upload()

    if args.symlink:
        modifier.symlink ( args.outfile )

    modifier.finalize()
