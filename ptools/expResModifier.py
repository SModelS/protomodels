#!/usr/bin/env python3

""" a simple class that performs global modifications on a list of results.
Used to ``take out potential signals'' i.e. put all observations to values
expected from background, by sampling the background model. """

# https://link.springer.com/content/pdf/10.1007/JHEP02(2015)004.pdf

import copy, os, sys, time, subprocess, math, numpy, shutil
import scipy.spatial
sys.path.insert( 0, "../" )
from csetup import setup
setup()
#sys.path.insert(0,"/scratch-cbe/users/wolfgan.waltenberger/git/protomodels/")
from scipy import stats
from builder.protomodel import ProtoModel
from builder.manipulator import Manipulator
from smodels.theory.model import Model
from smodels.share.models.SMparticles import SMList
from smodels.particlesLoader import BSMList
from smodels.theory.theoryPrediction import theoryPredictionsFor
from smodels.tools.simplifiedLikelihoods import Data, UpperLimitComputer
from smodels.theory import decomposer

class ExpResModifier:
    def __init__ ( self, dbpath, Zmax, rundir, keep, nproc, fudge,
                   suffix: str, lognormal = False, fixedsignals = False ):
        """
        :param dbpath: path to database
        :param Zmax: upper limit on an individual excess
        :param suffix: suffix to use, e.g. fake, signal, etc
        :param lognormal: if True, use lognormal for nuisances, else Gaussian
        :param fixedsignals: if True, then use the central value of theory prediction
                             as the signal yield, dont draw from Poissonian
        """
        self.comments = {} ## comments on entries in dict
        self.lognormal = lognormal
        self.dbpath = dbpath
        self.fixedsignals = fixedsignals
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

    def interact ( self, listOfExpRes ):
        import IPython
        IPython.embed( using=False )

    def extractStats ( self ):
        """ dont produce a new fake database, extract a stats dict
            from an existing database. """
        picklefile = self.rundir + "/" + self.dbpath
        if self.rundir in self.dbpath:
            picklefile = self.dbpath
        self.pprint ( f"Extracting stats from {picklefile}" )
        db = Database ( picklefile )
        self.dbversion = db.databaseVersion
        listOfExpRes = db.expResultList
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
                        D[txname]=txname.txnameData.y_values

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

    def computeNewObserved ( self, txname, globalInfo ):
        """ given expected upper limit, compute a fake observed limit
            by sampling the non-truncated Gaussian likelihood """
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
            D = {}
            while x > self.Zmax:
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

    def finalize ( self ):
        """ finalize, for the moment its just deleting slha files """
        print ( "[expResModifier] finalize" )
        if self.keep:
            return
        if hasattr ( self, "protomodel" ) and self.protomodel is not None:
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
        shutil.copyfile ( filename, self.rundir+"/signal.py" )
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

    def modifyDatabase ( self, outfile="", pmodel="" ):
        """ modify the database, possibly write out to a pickle file
        :param outfile: if not empty, write the database into file
        :param suffix: suffix to append to database version
        :param pmodel: if not empty, then this is the file name of the signal
                       model. in this case fake a signal
        :returns: the database
        """
        self.log ( "starting to create %s. suffix is %s protomodel is %s." % \
                   ( outfile, self.suffix, pmodel ) )
        db = Database ( self.dbpath )
        self.dbversion = db.databaseVersion
        # listOfExpRes = db.getExpResults( useSuperseded=True, useNonValidated=True )
        listOfExpRes = db.expResultList ## seems to be the safest bet?
        self.produceProtoModel ( pmodel, db.databaseVersion )
        # print ( "pm produced", os.path.exists ( self.protomodel.currentSLHA ) )
        self.log ( "%d results before faking bgs" % len(listOfExpRes) )
        updatedListOfExpRes = self.fakeBackgrounds ( listOfExpRes )
        # print ( "fb produced", os.path.exists ( self.protomodel.currentSLHA ) )
        self.log ( "%d results after faking bgs" % len(updatedListOfExpRes) )
        updatedListOfExpRes = self.addSignals ( updatedListOfExpRes )
        self.log ( "%d results after adding signals" % len(updatedListOfExpRes) )
        db.expResultList = updatedListOfExpRes
        newver = db.databaseVersion + self.suffix
        db.txt_meta.databaseVersion = newver
        db.pcl_meta.databaseVersion = newver
        self.pprint ( "Constructed fake database with %d (of %d) results" % \
                ( len(updatedListOfExpRes), len(listOfExpRes) ) )
        if outfile != "":
            db.createBinaryFile( outfile )
        return db

    def sampleEfficiencyMap ( self, dataset ):
        """ for the given dataset,
        sample from background and put the value as observed """
        orig = dataset.dataInfo.observedN
        exp = dataset.dataInfo.expectedBG
        err = dataset.dataInfo.bgError * self.fudge
        D = { "origN": orig, "expectedBG": exp, "bgError": err, "fudge": self.fudge }
        S, origS = float("inf"), float("nan")
        while S > self.Zmax:
            # lmbda = stats.norm.rvs ( exp, err )
            lmbda = self.drawNuisance ( exp, err )
            dataset.dataInfo.lmbda = lmbda
            if lmbda < 0.:
                lmbda = 0.
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
        D["S"]=S
        self.comments["S"]="the significance of the observation"
        D["origS"]=origS
        self.comments["origS"]="the significance of the original observation"
        self.comments["lmbda"]="Poissonian lambda of the fake background"
        D["lmbda"]=lmbda
        D["newObs"]=obs
        self.comments["newObs"]="the new fake observation"
        D["obsBg"]=obs
        self.comments["obsBg"]="the new fake observation, background component"
        D["toterr"]=toterr
        ## origN stores the n_observed of the original database
        dataset.dataInfo.origN = orig
        label = dataset.globalInfo.id + ":" + dataset.dataInfo.dataId
        self.addToStats ( label, D )
        return dataset

    def addSignalForEfficiencyMap ( self, dataset, tpred, lumi ):
        """ add a signal to this efficiency map. background sampling is
            already taken care of """
        txns = list ( map ( str, tpred.txnames ) )
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

    def addSignalForULMap ( self, dataset, tpred, lumi ):
        """ add a signal to this UL result. background sampling is
            already taken care of """
        from smodels.tools.physicsUnits import fb
        txns = list ( map ( str, tpred.txnames ) )
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
        D["txns"]=",".join(txns)
        self.comments["sigmaN"]="the added theory prediction (in fb), for UL maps"
        ## sigmaN is the predicted production cross section of the signal,
        ## in fb
        def distance ( v1, v2 ):
            """ compute distance between v1 and v2 """
            ret = 0.
            v1,v2 = list(v1),list(v2)
            if len(v1)*2 == len(v2):
                v1 = v1*2
            for _1,_2 in zip ( v1, v2 ):
                ret+= ( _1 - _2 )**2
            ret = math.sqrt (ret )
            return ret

        for i,txname in enumerate(dataset.txnameList):
            if not self.txNameIsIn ( txname, tpred ):
                continue
            #print ( "  `-- adding %s to %s" % ( sigmaN, txname ) )
            txnd = txname.txnameData
            etxnd = txname.txnameDataExp
            coordsTpred = txnd.dataToCoordinates ( tpred.mass, txnd._V, txnd.delta_x ) ## coordinates of tpred
            minDist = float("inf") ## for the closest point we store the numbers
            for yi,y in enumerate(txnd.y_values):
                pt = txnd.tri.points[yi] ## the point in the rotated coords
                dist = distance ( pt, coordsTpred )
                if dist > 400.: ## change y_values only in vicinity of protomodel
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
                    D={}## remove
                    D["yold"]=oldo
                    D["dist"]=dist
                    self.comments["dist"]="distance of closest point to protomodel"
                    #D["ptFXME"]=pt
                    #D["coordstpredFXME"]=coordsTpred
                    D["mass"]=str(tpred.mass) ## store as string
                    if hasExpected:
                        D["yexp"]=oldv
                        self.comments["yexp"]="expected y value (fb) closest to signal protomodel for UL map"
                    self.comments["yold"]="old y value (fb) closest to signal protomodel for UL map"
                    self.comments["ynew"]="new y value (fb) closest to signal protomodel for UL map"
                    D["ynew"]=oldv+sigmaN
                # print ( "    `--- adding %s %s" % ( oldv, sigmaN ) )
                txnd.y_values[yi]=oldv + sigmaN
            self.addToStats ( label, D )
            dataset.txnameList[i].txnameData = txnd
            dataset.txnameList[i].sigmaN = sigmaN
        return dataset

    def saveStats ( self ):
        """ write out the collected stats, so we can discuss experimentalists'
            conservativeness """
        filename = "%s/db%s.dict" % ( self.rundir, self.suffix )
        self.log ( f"saving stats to {filename}" )
        meta = { "dbpath": self.dbpath, "Zmax": self.Zmax,
                 "database": self.dbversion, "fudge": self.fudge,
                 "protomodel": '"%s"' % self.protomodel, "timestamp": time.asctime(),
                 "lognormal": self.lognormal }
        with open ( filename,"wt" ) as f:
            f.write ( str(meta)+"\n" )
            if len(self.comments)>0:
                f.write ( "# explanations on the used variables:\n" )
                f.write ( "# =====================================\n" )
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
                       remove_orig ):
        """ filter the list fo experimental results.
        :param outfile: store result in outfile (a pickle file)
        :param nofastlim: remove fastlim results
        :param onlyvalidated: remove non-validated results
        :param nosuperseded: remove superseded results
        :param remove_orig: remove original values
        """
        self.log ( "starting to filter %s. suffix is %s." % \
                   ( outfile, self.suffix ) )
        db = Database ( self.dbpath )
        listOfExpRes = db.expResultList ## seems to be the safest bet?
        newList = []
        for er in listOfExpRes:
            addThisOne = True
            if nofastlim:
                if hasattr ( er.globalInfo, "contact" ) and "fastlim" in er.globalInfo.contact:
                    print ( " `- skipping fastlim %s" % er.globalInfo.id )
                    addThisOne = False
            if nosuperseded:
                if hasattr ( er.globalInfo, "supersededBy" ):
                    print ( " `- skipping superseded %s" % er.globalInfo.id )
                    addThisOne = False
            if hasattr ( er.globalInfo, "private" ) and er.globalInfo.private in [ "True", True ]:
                    print ( " `- skipping private %s" % er.globalInfo.id )
                    addThisOne = False
            if not addThisOne:
                continue
            if onlyvalidated:
                newDs = []
                for ds in er.datasets:
                    txnew = []
                    for txn in ds.txnameList:
                        if txn.validated == False:
                            print ( " `- skipping non-validated %s/%s/%s" % \
                                    ( txn.txName, ds.dataInfo.dataId, er.globalInfo.id ) )
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
                for iD,ds in enumerate(er.datasets):
                    for it,txn in enumerate(ds.txnameList):
                        #txn.txnameData = self.cleanTxNameData ( txn.txnameData )
                        for label in [ "figureUrl", "dataUrl" ]:
                            if hasattr ( txn, label ):
                                delattr ( txn, label )
                        if hasattr ( txn.txnameData, "origdata" ):
                            del er.datasets[iD].txnameList[it].txnameData.origdata
                        if txn.txnameDataExp != None:
                            #txn.txnameDataExp = self.cleanTxNameData ( txn.txnameDataExp )
                            if hasattr ( txn.txnameDataExp, "origdata" ):
                                del er.datasets[iD].txnameList[it].txnameDataExp.origdata
            if not addThisOne:
                continue
            newList.append ( er )
        db.expResultList = newList
        if outfile != "":
            db.createBinaryFile( outfile )

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
        db = Database ( picklefile )
        listOfExpRes = db.getExpResults()
        for er in listOfExpRes:
            datasets = er.datasets
            for ds in datasets:
                txnl = ds.txnameList
                for txn in txnl:
                    x = txn.txnameData.dataType
        print ( "were good", db.databaseVersion )


epilog="""
Examples:
=========

Fake SM-only database:
----------------------
./expResModifier.py -R $RUNDIR -d original.pcl -s fake1

Database with a fake signal:
----------------------------
./expResModifier.py -R $RUNDIR -d original.pcl -s signal1 -P pmodel9.py

Build a database:
-----------------
./expResModifier.py -B -d ../../../smodels-database

Just filter the database:
-------------------------
./expResModifier.py -d ./original.pcl --remove_orig --nofastlim --onlyvalidated --nosuperseded --dontsample -o test.pcl

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
    argparser.add_argument ( '--dontsample',
            help='do not sample at all, only filter',
            action='store_true' )
    argparser.add_argument ( '-l', '--lognormal',
            help='use lognormal, not Gaussian for nuisances',
            action='store_true' )
    argparser.add_argument ( '--fixedsignals',
            help='fix the contributions from the signals, dont draw from Poissonian',
            action='store_true' )
    argparser.add_argument ( '-M', '--max',
            help='upper limit on significance of individual excess [None]',
            type=float, default=None )
    argparser.add_argument ( '--seed',
            help='set a random number seed [None]',
            type=int, default=None )
    argparser.add_argument ( '-N', '--nproc',
            help='number of parallel processes, for signal adding [1]',
            type=int, default=1 )
    argparser.add_argument ( '-P', '--pmodel',
            help='supply filename of a pmodel, in which case create a signal-infused database [""]',
            type=str, default="" )
    argparser.add_argument ( '-v', '--verbose',
            help='print results to stdout', action='store_true' )
    argparser.add_argument ( '-I', '--interactive',
            help='interactive mode', action='store_true' )
    argparser.add_argument ( '-B', '--build',
            help='build the original pickle file with all relevant info, then exit (use --database to specify path)', action='store_true' )
    argparser.add_argument ( '-c', '--check',
            help='check the pickle file <outfile>', action='store_true' )
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
    if args.seed is not None:
        from ptools import helpers
        helpers.seedRandomNumbers( args.seed )
        print ( f"[expResModifier] setting random seed to {args.seed}" )
    if args.build:
        from smodels.experiment.txnameObj import TxNameData
        TxNameData._keep_values = True
        from smodels.experiment.databaseObj import Database
        print ( f"[expResModifier] starting to build database at {args.database}." )
        db = Database ( args.database )
        print ( f"[expResModifier] built database at {args.database}. Exiting." )
        sys.exit()
    if type(args.rundir)==str and not "/" in args.rundir:
        args.rundir = "/scratch-cbe/users/wolfgan.waltenberger/" + args.rundir
    if args.outfile == "":
        args.outfile = args.suffix+".pcl"
    from smodels.experiment.databaseObj import Database
    modifier = ExpResModifier( args.database, args.max, args.rundir, args.keep, \
                               args.nproc, args.fudge, args.suffix, args.lognormal,
                               args.fixedsignals )

    if not args.outfile.endswith(".pcl"):
        print ( "[expResModifier] warning, shouldnt the name of your outputfile ``%s'' end with .pcl?" % args.outfile )
    if args.nofastlim or args.onlyvalidated or args.nosuperseded or args.remove_orig:
        modifier.filter ( args.outfile, args.nofastlim, args.onlyvalidated,
                          args.nosuperseded, args.remove_orig )
        modifier.symlink ( args.outfile )

    if args.dontsample:
        print ( "[expResModifier] we were asked to not sample, so we exit now." )
        sys.exit()

    if args.extract_stats:
        er = modifier.extractStats()
    else:
        er = modifier.modifyDatabase ( args.outfile, args.pmodel )

    modifier.saveStats()

    if args.check:
        modifier.check ( args.outfile )

    if args.interactive:
        modifier.interact ( er )

    if args.upload:
        modifier.upload()

    modifier.symlink ( args.outfile )
    modifier.finalize()
