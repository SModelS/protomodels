#!/usr/bin/env python3

""" a first start at the random walk idea """

__all__ = [ "RandomWalker" ]

import pickle, sys, time, math, socket, os
import subprocess
import numpy, colorama

sys.path.insert(0,f"{os.environ['HOME']}/git/protomodels/")
try:
    sys.path.insert(0,f"{os.environ['HOME']}/git/smodels/")
    import smodels
except:
    from ptools import setPath
sys.path.insert(0,f"/scratch-cbe/users/{os.environ['USER']}/git/smodels-utils/protomodels/")
sys.path.insert(0,"../")
sys.path.insert(0,"../../")
sys.path.insert(0,"../smodels/")
from walker.hiscores import Hiscores
from builder.protomodel import ProtoModel
from builder.manipulator import Manipulator
from tester.predictor import Predictor
from tester.critic import Critic
from ptools.sparticleNames import SParticleNames
from pympler.asizeof import asizeof
from smodels.base.smodelsLogging import logger
from typing import Callable, Dict, Union
from os import PathLike
from base.loggerbase import LoggerBase
from ptools.helpers import prettyPrint

logger.setLevel("ERROR")

def __cleanDirectory ():
    ## is this used?
    subprocess.getoutput ( "mkdir -p tmp" )
    subprocess.getoutput ( "mv .cur*slha tmp/" )
    subprocess.getoutput ( "mv walker*.log tmp/" )
    subprocess.getoutput ( "mv exceptions.log tmp/" )

class RandomWalker ( LoggerBase ):
    def __init__ ( self, walkerid : int = 0, nsteps : int = 10000,
            strategy : str = "aggressive",
            cheatcode : int = 0, dbpath : PathLike = "./database.pcl",
            expected : bool = False, select : str = "all",
            catch_exceptions : bool = True, rundir : Union[PathLike,None] = None,
            do_srcombine : bool = False,
            record_history : bool = False, seed : Union[int,None] = None,
            stopTeleportationAfter : int = -1 ):
        """ initialise the walker
        :param nsteps: maximum number of steps to perform, negative is infinity
        :param cheatcode: cheat mode. 0 is no cheating, 1 is with ranges, 2
                      is the Z323 model.
        :param expected: remove possible signals from database
        :param select: select only subset of results (all for all, em for
                efficiency maps only, ul for upper limits only, alternatively
                select for txnames via e.g. "txnames:T1,T2"
        :param catch_exceptions: should we catch exceptions
        :param do_srcombine: if true, then also perform combinations, either via
                           simplified likelihoods or via pyhf
        :param record_history: if true, attach a history recorder class
        :param seed: random seed, int or None
        :param stopTeleportationAfter: int or None. we stop teleportation after
                this step nr.  If negative or None, we dont teleport at all
        """

        #call the super class of the random walker i.e Loggerbase
        super ( RandomWalker, self ).__init__ ( walkerid )
        dbpath = os.path.expanduser ( dbpath )
        if type(walkerid) != int or type(nsteps) != int or type(strategy)!= str:
            self.pprint ( f"Wrong call of constructor: {walkerid}, {nsteps}, {strategy}" )
            sys.exit(-2)
        self.walkerid = walkerid ## walker id, for parallel runs
        self.rundir = rundir
        if rundir == None:
            self.rundir = "./"

        if seed is not None:
            from ptools import helpers
            helpers.seedRandomNumbers(seed + walkerid )
            self.pprint ( f"setting random seed to {seed}" )

        #Initialize Predictor
        self.predictor =  Predictor( self.walkerid, dbpath=dbpath,
                              expected=expected, select=select, do_srcombine=do_srcombine )
        self.critic =  Critic( self.walkerid, dbpath=dbpath,
                              expected=expected, select=select, do_srcombine=do_srcombine )

        #Initialize Hiscore (with access to the predictor)
        picklefile = f"{self.rundir}/H{walkerid}.cache"
        save_hiscores = True
        self.hiscoreList = Hiscores ( walkerid, save_hiscores=save_hiscores,
                picklefile=picklefile, backup=False, predictor=self.predictor )
        self.hiscoreList.nkeep = 1

        #Initialize ProtoModel and Manipulator:
        protomodel = ProtoModel( self.walkerid, keep_meta = True,
                dbversion = self.predictor.database.databaseVersion )

        self.manipulator = Manipulator ( protomodel, strategy,
                        do_record = record_history, seed = seed )
        self.catch_exceptions = catch_exceptions
        self.maxsteps = nsteps
        if stopTeleportationAfter == None:
            stopTeleportationAfter = -1
        # stopTeleportationAfter = self.maxsteps/3.
        self.stopTeleportationAfter = stopTeleportationAfter
        if record_history:
            from ptools.history import History
            self.recorder = History ( f"{self.rundir}/history{walkerid}.list" )
            self.manipulator.do_record = True
        jobid = "unknown"
        if "SLURM_JOBID" in os.environ:
            jobid = os.environ["SLURM_JOBID"]
        self.pprint ( f"Ramping up with slurm jobid {jobid}" )

        if cheatcode <= 0:
            self.takeStep() # the first step should be considered as "taken"
            #Set current TL and K values to threshold values
            self.currentTL = -0.1
            self.currentK = -20.0
        else:
            self.manipulator.cheat ( cheatcode )
            #self.predictor.predict(self.protomodel)
            self.predict(self.manipulator)
            if type(self.manipulator.M.TL) != type(None):
                self.pprint ( f"Cheat model gets TL={self.manipulator.M.TL:.2f}, "\
                              f"K={self.manipulator.M.K:.2f}" )
            # self.printStats ( substep=4 )
            self.manipulator.backupModel()
            self.hiscoreList.newResult ( self.manipulator )
            self.printStats ( substep=5 )
            self.currentK = self.manipulator.M.K
            self.currentTL = self.manipulator.M.TL

    def setWalkerId ( self, Id ):
        self.walkerid = Id
        self.manipulator.setWalkerId ( Id )

    @classmethod
    def fromProtoModel( cls, protomodel : ProtoModel, **args : Dict ):
        """ create a RandomWalker from a ProtoModel. Continue walking
            from that model """
        ret = cls( **args )
        ret.manipulator.M = protomodel
        if "walkerid" in args:
            ret.manipulator.setWalkerId ( args["walkerid"] )
        ret.manipulator.backupModel()
        return ret

    def extractArguments ( func : Callable, **args : Dict ) -> Dict:
        """ from args, extract all the entries that are parameters of func """
        pm = {}
        import inspect
        sig = inspect.signature ( func )
        for i in sig.parameters.keys():
            if i != "self" and i in args:
                pm[i]=args[i]
        return pm

    @classmethod
    def fromDictionary( cls, dictionary : Union[PathLike,Dict], **args : Dict ):
        """ create a RandomWalker from a hiscore dictionary. Continue walking
            from the model in that dictionary
        :param dictionary: either a dictionary, or a string containing a dictionary,
        or the path to a dictionary
        """
        if type(dictionary) == str and dictionary.endswith ( ".dict" ):
            if not os.path.exists ( dictionary ):
                logger.error ( f"argument {dictionary} is a string, but doesnt work as pathname" )
                sys.exit()
            try:
                logger.info ( f"trying to interpret {dictionary} as a path... " )
                f = open ( dictionary, "rt" )
                tmp = eval ( f.read() )
                f.close()
                if type(tmp) == dict:
                    dictionary = tmp
                if type(tmp) == list and type(tmp[0])==dict:
                    dictionary = tmp[0]
                    logger.info ( f" ... seemed to work!" )
            except Exception as e:
                logger.error  ( f"could not interpret the content of {dictionary}: {e}" )

        ret = cls( **args ) ## simply pass on all the arguments

        pm = RandomWalker.extractArguments ( ProtoModel.__init__, **args )
        ret.manipulator.M = ProtoModel( **pm )
        ret.manipulator.initFromDict ( dictionary )
        if "walkerid" in args:
            ret.manipulator.setWalkerId ( args["walkerid"] )
        ret.manipulator.M.createNewSLHAFileName()
        # ret.printStats ( substep=3 )
        ret.manipulator.backupModel()
        return ret

    @property
    def protomodel(self):
        return self.manipulator.M

    @protomodel.setter
    def protomodel(self, protomodel):
        self.manipulator.M = protomodel

    def printStats ( self, substep=None ):
        """ print the stats, i.e. number of unfrozen particles.
            for debugging. """
        #print(f"best Combo {self.manipulator.M.bestCombo}, step {self.manipulator.M.step}")
        nUnfrozen = len( self.protomodel.unFrozenParticles() )
        nTotal = len ( self.protomodel.particles )
        pidsp = self.protomodel.unFrozenParticles()
        pidsp.sort()
        namer = SParticleNames ( False )

        prtcles = ", ".join ( map ( namer.asciiName, pidsp ) )
        if self.manipulator.M.bestCombo:
            pidsbc = list ( self.manipulator.getAllPidsOfBestCombo() )
            pidsbc.sort()
            prtclesbc = ", ".join ( map ( namer.asciiName, pidsbc ) )
            self.pprint ( f"Step {self.protomodel.step} has {nUnfrozen}/{nTotal} unfrozen particles: {prtcles} [in best combo: {prtclesbc}]" )
            if len(pidsbc)>0 and not set(pidsbc).issubset ( set(pidsp) ):
                self.pprint ( f"  `-- error! best combo pids ({pidsbc}) arent subset of masses pids ({pidsp})!" )
                self.manipulator.M.bestCombo = None

    def predict ( self, manipulator : Manipulator ):
        """ Calls predictor.predict to get the theory predictions for model. Loops for 5 times till model.muhat is close to 1.0 """
        #print(f"Adress of manip : {id(manipulator)}")
        model = manipulator.M
        muhat_converge = False
        for i in range(5):
            predict = self.predictor.predict(model)
            if predict: #returns False if no preds are found or TL is None (i.e no comb found)
                #print(f"i {i}, muhat {model.muhat}, convergence {abs(model.muhat - 1.0)}")
                if abs(model.muhat - 1.0) < 1e-02:
                    self.pprint(f"Step {model.step} converged at loop {i} with muhat {model.muhat}!")
                    muhat_converge = True
                    break
                manipulator.rescaleSignalBy(model.muhat) #?
            else:
                break # Rescale signal by a significant number?

        if not muhat_converge:  #reverting step
            if predict:
                self.pprint ( f"Step {model.step} did not converge to muhat 1.0, model muhat is {model.muhat}. Going back to previous step." )
            else:
                self.pprint ( f"Step {model.step} did not converge to muhat 1.0. Model did not find any prediction." )
            return False

        return True

    def onestep ( self ):
        #Add one step
        self.protomodel.step+=1
        self.pprint ( "Step %d begins." % ( self.protomodel.step ) )
        self.printStats( )
        #Remove data about best combo
        self.protomodel.cleanBestCombo()
        # self.printStats( substep=11 )
        printMemUsage = False
        if printMemUsage:
            self.pprint ( f"memory footprint (kb): walker {asizeof(self)/1024}, model {asizeof(self.protomodel)/1024}" )

        #Trim the model, so we start only with the relevant particles for the
        #best combination in the previous step:
        if self.protomodel.bestCombo:
            self.log ( "freeze pids that arent in best combo, we dont need them:" )
            nfrozen = self.manipulator.freezePidsNotInBestCombo()
            self.log ( " `- froze %d particles not in best combo" % nfrozen )
        # self.printStats( substep=12 )

        #Take a step in the model space:
        self.manipulator.randomlyChangeModel()
        # self.printStats( substep=13 )

        nUnfrozen = len( self.protomodel.unFrozenParticles() )
        ## number of pids in best combo, as a check

        #Try to create a simpler model
        #(merge pre-defined particles if their mass difference is below dm)
        protomodelSimp = self.manipulator.simplifyModel(dm=200.0)
        manipulatorSimp = None
        if protomodelSimp: manipulatorSimp = Manipulator ( protomodelSimp, strategy="aggressive",do_record = False, seed = 0 )
        boolProtoSimp = False

        # self.printStats( substep=14 )

        if self.catch_exceptions:
            try:
                if not self.predict(self.manipulator):
                    #print("return")
                    self.protomodel.K = None
                    return #??
                if protomodelSimp:
                    #print(f"Address of manip before call {id(manipulatorSimp)}")
                    boolProtoSimp = self.predict(manipulatorSimp) #!rewrite
            except Exception as e:
                self.pprint ( f"@@@ caught exception @@@" )
                self.pprint ( f"{type(e)} ``{str(e)}'' encountered when trying to predict. lets revert and not count it as a step." )
                import traceback
                self.pprint ( f"traceback says:: {traceback.format_exc()}" )
                self.pprint ( f"model is:" )
                d = self.manipulator.writeDictFile(None)
                self.pprint ( f"{str(d)}" )
                self.pprint ( f"@@@ end exception @@@" )
                if False:
                    import tempfile
                    f = tempfile.mktemp ( suffix=".slha", prefix="failed", dir="./" )
                    self.manipulator.M.writeSLHAFile ( f )
                self.manipulator.restoreModel()
                self.manipulator.M.step -= 1 # we dont count that step.
                import traceback ## FIXME print to file!!
                traceback.print_exc()
                return
        else:
            if not self.predict(self.manipulator):
                self.protomodel.K = None
                return      #??
            if protomodelSimp:
                boolProtoSimp = self.predict(manipulatorSimp)

        #Now keep the model with highest score:
        if protomodelSimp and boolProtoSimp:
            if self.manipulator.M.K is None or (protomodelSimp.K is not None
                        and (protomodelSimp.K > self.manipulator.M.K)):
                self.manipulator.M = protomodelSimp


        #If no combination could be found, return
        if self.manipulator.M.TL is None or self.manipulator.M.K is None:
            return

        if len(self.manipulator.M.rvalues) > 1:
            self.log ( "Top r values are: %.2f, %.2f" % \
                       ( self.manipulator.M.rvalues[0], self.manipulator.M.rvalues[1] ) )

        self.log ( "Step %d: found highest TL: %.2f" % \
                   ( self.protomodel.step, self.protomodel.TL ) )

        nUnfrozen = len ( self.protomodel.unFrozenParticles() )
        self.pprint ( "best combo for strategy ``%s'' is %s: %s: [K=%.2f, TL=%.2f, %d unfrozen]" % \
            ( self.manipulator.strategy, self.protomodel.letters, self.protomodel.description, self.protomodel.K, self.protomodel.TL, nUnfrozen ) )
        smaxstp = "%s" % self.maxsteps
        if self.maxsteps < 0:
            smaxstp = "inf"

        #For low scoring models, teleport to a high score model:
        if self.checkIfToTeleport( pmax=0.5, norm = 10.0 ):
            # if we teleport the rest becomes irrelevant
            return
        self.printStats( )

        self.log ( "Step %d check if result goes into hiscore list" % \
                   ( self.protomodel.step ) )
        srs = "%s" % ", ".join ( [ "%.2f" % x for x in self.protomodel.rvalues[:3] ] )
        self.log ( "r values before calling .newResult are at %s" % srs )
        self.hiscoreList.newResult ( self.manipulator ) ## add to high score list
        srs = "%s" % ", ".join ( [ "%.2f" % x for x in self.protomodel.rvalues[:3] ] )
        self.log ( "r values after calling .newResult are at %s" % srs )
        self.log ( "done check for result to go into hiscore list" )
        self.log ( "Step %d/%s finished." % ( self.protomodel.step, smaxstp ) )

    def checkIfToTeleport ( self, pmax=0.1, norm = 10.0 ):
        """ check if we should teleport to a high score model. If yes, then we
            should then also perform the teleportation. The teleportation is
            done only if the model has a score smaller then the best score in
            hiscoreList.  The teleportation probability is given by
            pmax*(1-exp^(K-bestK)/norm), so pmax is the maximum probability
            (when K -> -infinity).

        :param pmax: Maximum probability for teleportation.
        :param norm: Normalization for K distance.
        """
        if self.protomodel.step > self.stopTeleportationAfter:
            self.log ( "teleportation is turned off after step #%d" % \
                       self.stopTeleportationAfter )
            return False
        #self.log ( "teleportation turned off" )
        #return False
        import random
        bestK = self.hiscoreList.globalMaxK()
        if bestK < 1.:
            self.log ( "bestK is smaller than one. no teleporting." )
            return False
        ourK = -2.
        if hasattr ( self.manipulator.M, "K" ) and self.manipulator.M.K > -2:
            ourK = self.manipulator.M.K
        #The current model already is the best, do nothing.
        if ourK >= bestK:
            return False
        #Otherwise compute the teleportation probability:
        dK = ( ourK - bestK ) / norm
        prob = pmax*(1. - math.exp( dK ))
        a = random.uniform ( 0., 1. )
        doTP = ( a < prob ) ## do teleport, yes or no
        sDoTP = "a>p: dont teleport."
        if doTP:
            sDoTP = "a<p: do teleport."
        self.log ( "check if to teleport, Kmax=%.2f, ours is=%.2f, p=%.2f, a=%.2f, %s" % \
                   ( bestK, ourK, prob, a, sDoTP ) )
        if doTP:
            self.manipulator.teleportToHiscore()
        return doTP

    def takeStep ( self ):
        """ take the step, save it as last step """
        ## Backup model
        self.manipulator.backupModel()
        # Update current K and TL values
        self.currentK = self.protomodel.K
        self.currentTL = self.protomodel.TL
        self.manipulator.record( "take step" )

    def decideOnTakingStep ( self ):
        """ depending on the ratio of K values, decide on whether to take the step or not.
            If ratio > 1., take the step, if < 1, let chance decide. """
        K = self.currentK
        if K == None: # if the old is none, we do everything
            self.takeStep()
            return

        newK = self.protomodel.K
        if newK == None:
            # if the new is none, but the old isnt, we go back
            self.manipulator.restoreModel( reportReversion=True )
            return

        if newK > K:
            self.highlight ( "info", f"K: {prettyPrint(K)} -> {prettyPrint(newK)}: check critics." )

            if self.critic.predict_critic(self.protomodel, keep_predictions=True):
                self.pprint ( "Passed both critics, taking the step." )
                self.takeStep()
            else:
                self.pprint ( "Failed at least one critic, the step is reverted." )
                self.manipulator.restoreModel( reportReversion=True )

        else:
            import random

            u = random.uniform(0.,1.)
            ratio = numpy.exp(.5*( newK - K))
            if u > ratio:
                self.pprint ( f"u={u:.2f} > {ratio:.2f}; K: {prettyPrint(K)} -> {prettyPrint(newK)}: revert." )
                self.manipulator.restoreModel( reportReversion=True )
            else:
                self.highlight ( "info", f"K: {prettyPrint(K)} -> {prettyPrint(newK)}; u={u:.2f} <= {ratio:.2f}: check critics." )   #SN: <+ and not > right?

                if self.critic.predict_critic(self.protomodel, keep_predictions=True):
                    self.pprint ( "Passed both critics, taking the step." )
                    self.takeStep()
                else:
                    self.pprint ( "Failed at least one critic, the step is reverted." )
                    self.manipulator.restoreModel( reportReversion=True )

    def record ( self ):
        """ if recorder is defined, then record. """
        ## do we have a history recorder?
        if not hasattr ( self, "recorder" ):
            return
        self.recorder.add ( self.manipulator )

    def walk ( self, catchem=False ):
        """ Now perform the random walk """
        # self.printStats ( substep = 2 )
        self.manipulator.backupModel()
        if len ( self.manipulator.M.unFrozenParticles( withLSP=False ) ) < 1:
            ## start with unfreezing a random particle
            self.manipulator.randomlyUnfreezeParticle(force = True)
            self.manipulator.backupModel()

        while self.maxsteps < 0 or self.protomodel.step<self.maxsteps:

            if not catchem:
                self.onestep()
            else:
                try:
                    self.onestep()
                except Exception as e:
                    # https://bioinfoexpert.com/2016/01/18/tracing-exceptions-in-multiprocessing-in-python/
                    self.pprint ( "taking a step resulted in exception: %s, %s" % \
                                  (type(e), e ) )
                    import traceback
                    traceback.print_stack( limit=None )
                    except_type, except_class, tb = sys.exc_info()
                    extracted = traceback.extract_tb(tb)
                    for point in extracted:
                        self.pprint ( f"extracted: {point}" )
                    with open( f"{self.rundir}/exceptions.log","a") as f:
                        f.write ( f"{time.asctime()}: taking a step resulted in exception: {type(e)}, {e}\n" )
                        f.write ( f"   `- exception occured in walker #{self.protomodel.walkerid}\n" )
                        import traceback
                        f.write ( f"traceback: {str(traceback.format_exc())}\n" )
                    sys.exit(-1)

            #If no combination was found, go back
            if self.protomodel.K is None:
                print("returned")
                self.manipulator.restoreModel(reportReversion=True)
                continue

            # obtain the ratio of posteriors
            self.decideOnTakingStep ()
            self.record()
        self.manipulator.M.delCurrentSLHA()
        self.pprint ( "Was asked to stop after %d steps" % self.maxsteps )

if __name__ == "__main__":
    masses = { 1000022: 47.1, 1000023: 449.8, 1000024: 135.47, 1000037: 398.7,
               1000006: 543.7 }
    ssms = {(1000022, 1000022): 0.188322, (1000023, 1000023): 0.208443, (1000022, 1000023): 0.213852, (1000024, 1000024): 1.402063, (-1000024, 1000024): 1.402063, (-1000024, -1000024): 1.402063, (1000022, 1000024): 1.173816, (-1000024, 1000022): 1.173816, (1000023, 1000024): 1.363297, (-1000024, 1000023): 1.363297, (1000037, 1000037): 1.414218, (-1000037, 1000037): 1.414218, (-1000037, -1000037): 1.009287, (1000022, 1000037): 1.257227, (-1000037, 1000022): 1.257227, (1000023, 1000037): 2.62731, (-1000037, 1000023): 2.62731, (1000024, 1000037): 0.0, (-1000037, 1000024): 0.0, (-1000024, 1000037): 0.0, (-1000037, -1000024): 0.0, (1000006, 1000006): 0.050047, (-1000006, 1000006): 0.050047, (-1000006, -1000006): 0.050047, (1000006, 1000022): 0.055586, (-1000006, 1000022): 0.055586}
    decays = {1000022: {}, 1000023: {(1000022, 25): 1.0}, 1000024: {(1000022, 24): 1.0}, 1000037: {(1000022, 24): 0.353533, (1000024, 23): 0.646467}, 1000006: {(1000022, 6): 1.0}}
    D = {'masses': masses, 'ssmultipliers': ssms, 'decays': decays }
    dbpath = "official"
    select = "txnames:electroweakinos,stops"
    select = "all"
    walker = RandomWalker.fromDictionary ( D, walkerid = 0, dbpath = dbpath,
            do_srcombine = True, select = select )
    walker.walk()
