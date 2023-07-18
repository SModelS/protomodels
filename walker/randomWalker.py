#!/usr/bin/env python3

""" a first start at the random walk idea """

__all__ = [ "RandomWalker" ]

import pickle, sys, time, math, socket, os
import subprocess
import numpy, colorama

try:
    import smodels
except:
    from ptools import setPath
sys.path.insert(0,f"/scratch-cbe/users/{os.environ['USER']}/git/smodels-utils/protomodels/")
sys.path.insert(0,"../")
from walker.hiscore import Hiscore
from builder.protomodel import ProtoModel
from builder.manipulator import Manipulator
from tester.predictor import Predictor
from ptools.sparticleNames import SParticleNames
from pympler.asizeof import asizeof
from smodels.base.smodelsLogging import logger
from typing import Callable, Dict, Union
from os import PathLike

logger.setLevel("ERROR")

def __cleanDirectory ():
    ## is this used?
    subprocess.getoutput ( "mkdir -p tmp" )
    subprocess.getoutput ( "mv .cur*slha tmp/" )
    subprocess.getoutput ( "mv walker*.log tmp/" )
    subprocess.getoutput ( "mv exceptions.log tmp/" )

class RandomWalker:
    def __init__ ( self, walkerid : int = 0, nsteps : int = 10000, 
            strategy : str = "aggressive", dump_training : bool = False, 
            cheatcode : int = 0, dbpath : PathLike = "./database.pcl", 
            expected : bool = False, select : str = "all", 
            catch_exceptions : bool = True, rundir : Union[PathLike,None] = None, 
            nevents : int = 100000, do_combine : bool = False, 
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
        :param nevents: number of MC events when computing cross-sections
        :param do_combine: if true, then also perform combinations, either via
                           simplified likelihoods or via pyhf
        :param record_history: if true, attach a history recorder class
        :param seed: random seed, int or None
        :param stopTeleportationAfter: int or None. we stop teleportation after
                this step nr.  If negative or None, we dont teleport at all
        """
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
                              expected=expected, select=select, do_combine=do_combine )

        #Initialize Hiscore (with access to the predictor)
        picklefile = f"{self.rundir}/H{walkerid}.hi"
        save_hiscores = True
        self.hiscoreList = Hiscore ( walkerid, save_hiscores=save_hiscores, 
                picklefile=picklefile, backup=False, predictor=self.predictor )
        self.hiscoreList.nkeep = 1

        #Initialize ProtoModel and Manipulator:
        protomodel = ProtoModel( self.walkerid, keep_meta = True, nevents = nevents, 
                dbversion = self.predictor.database.databaseVersion )

        self.manipulator = Manipulator ( protomodel, strategy, 
                        do_record = record_history, seed = seed )
        self.catch_exceptions = catch_exceptions
        self.maxsteps = nsteps
        if stopTeleportationAfter == None:
            stopTeleportationAfter = -1
        # stopTeleportationAfter = self.maxsteps/3.
        self.stopTeleportationAfter = stopTeleportationAfter
        self.accelerator = None
        if record_history:
            from ptools.history import History
            self.recorder = History ( f"{self.rundir}/history{walkerid}.list" )
            self.manipulator.do_record = True
        jobid = "unknown"
        if "SLURM_JOBID" in os.environ:
            jobid = os.environ["SLURM_JOBID"]
        self.pprint ( "Ramping up with slurm jobid %s" % jobid )

        if cheatcode <= 0:
            self.takeStep() # the first step should be considered as "taken"
            #Set current Z and K values to threshold values
            self.currentZ = -0.1
            self.currentK = -20.0
        else:
            self.manipulator.cheat ( cheatcode )
            self.predictor.predict(self.protomodel)
            self.pprint ( "Cheat model gets Z=%.2f, K=%.2f" % \
                          ( self.manipulator.M.Z, self.manipulator.M.K ) )
            # self.printStats ( substep=4 )
            self.manipulator.backupModel()
            self.hiscoreList.newResult ( self.manipulator )
            self.printStats ( substep=5 )
            self.currentK = self.manipulator.M.K
            self.currentZ = self.manipulator.M.Z
        if dump_training:
            from accelerator import Accelerator
            ## we use the accelerator only to dump the training data
            self.accelerator = Accelerator ( walkerid= walkerid,
                                dump_training= True,
                                is_trained = False  )

    def hostname ( self ):
        return socket.gethostname()

    def setWalkerId ( self, Id ):
        self.walkerid = Id
        self.manipulator.setWalkerId ( Id )
        if self.accelerator != None:
            self.accelerator.walkerid = Id

    @classmethod
    def fromProtoModel( cls, protomodel : ProtoModel, **args : Dict ):
        """ create a RandomWalker from a ProtoModel. Continue walking
            from that model """
        ret = cls( **args )
        ret.manipulator.M = protomodel
        if "walkerid" in args:
            ret.manipulator.setWalkerId ( args["walkerid"] )
        ret.manipulator.backupModel()
        if dump_training:
            ## we use the accelerator only to dump the training data
            from accelerator import Accelerator
            ret.accelerator = Accelerator ( walkerid= walkerid, dump_training=True,
                                        is_trained = False )
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

    def pprint ( self, *args ):
        """ logging """
        if not hasattr ( self, "walkerid" ):
            self.walkerid=-1
        print ( "[walk:%d:%s-%s] %s" % ( self.walkerid, self.hostname(), time.strftime("%H:%M:%S"), " ".join(map(str,args))) )
        self.log ( *args )

    @property
    def protomodel(self):
        return self.manipulator.M

    @protomodel.setter
    def protomodel(self, protomodel):
        self.manipulator.M = protomodel

    def printStats ( self, substep=None ):
        """ print the stats, i.e. number of unfrozen particles.
            for debugging. """
        nUnfrozen = len( self.protomodel.unFrozenParticles() )
        nTotal = len ( self.protomodel.particles )
        pidsp = self.protomodel.unFrozenParticles()
        pidsp.sort()
        namer = SParticleNames ( False )

        prtcles = ", ".join ( map ( namer.asciiName, pidsp ) )
        pidsbc = list ( self.manipulator.getAllPidsOfBestCombo() )
        pidsbc.sort()
        prtclesbc = ", ".join ( map ( namer.asciiName, pidsbc ) )
        self.pprint ( "Step %d has %d/%d unfrozen particles: %s [in best combo: %s]" % \
              ( self.protomodel.step, nUnfrozen, nTotal, \
                prtcles, prtclesbc ) )
        if len(pidsbc)>0 and not set(pidsbc).issubset ( set(pidsp) ):
            self.pprint ( "  `-- error! best combo pids arent subset of masses pids!!!" )
            self.manipulator.M.bestCombo = None

    def predict ( self, model : Union [ ProtoModel, None ] = None ):
        """ convenience function """
        if model == None:
            self.predictor.predict(self.manipulator.M)
        else:
            self.predictor.predict(model)

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
            self.pprint ( "memory footprint (kb): walker %d, model %d, accelerator %d" %\
                    ( asizeof(self)/1024,asizeof(self.protomodel)/1024,asizeof(self.accelerator)/1024 ))

        #Trim the model, so we start only with the relevant particles for the
        #best combination in the previous step:
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
        #(merge pre-defined particles of their mass difference is below dm)
        protomodelSimp = self.manipulator.simplifyModel(dm=200.0)

        # self.printStats( substep=14 )

        if self.catch_exceptions:
            try:
                self.predict()
                if protomodelSimp:
                    self.predict(protomodelSimp)
            except Exception as e:
                self.pprint ( f"{type(e)} ``{str(e)}'' encountered when trying to predict. lets revert and not count it as a step." )
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
            self.predict()
            if protomodelSimp:
                self.predict(protomodelSimp)

        #Now keep the model with highest score:
        if protomodelSimp:
            if self.manipulator.M.K is None or (protomodelSimp.K is not None
                        and (protomodelSimp.K > self.manipulator.M.K)):
                self.manipulator.M = protomodelSimp


        #If no combination could be found, return
        if self.manipulator.M.Z is None:
            return

        #the muhat multiplier gets multiplied into the signal strengths
        self.manipulator.rescaleSignalBy(self.protomodel.muhat)

        self.log ( "Top r values after rescaling are: %.2f, %.2f" % \
                   ( self.manipulator.M.rvalues[0], self.manipulator.M.rvalues[1] ) )

        self.log ( "Step %d: found highest Z: %.2f" % \
                   ( self.protomodel.step, self.protomodel.Z ) )

        nUnfrozen = len ( self.protomodel.unFrozenParticles() )
        self.pprint ( "best combo for strategy ``%s'' is %s: %s: [K=%.2f, Z=%.2f, %d unfrozen]" % \
            ( self.manipulator.strategy, self.protomodel.letters, self.protomodel.description, self.protomodel.K, self.protomodel.Z, nUnfrozen ) )
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
        # Update current K and Z values
        self.currentK = self.protomodel.K
        self.currentZ = self.protomodel.Z
        self.manipulator.record( "take step" )

    def highlight ( self, msgType = "info", *args ):
        """ logging, hilit """
        col = colorama.Fore.GREEN
        print ( "%s[walk:%d] %s%s" % ( col, self.walkerid, " ".join(map(str,args)), colorama.Fore.RESET ) )

    def decideOnTakingStep ( self ):
        """ depending on the ratio of K values, decide on whether to take the step or not.
            If ratio > 1., take the step, if < 1, let chance decide. """
        import random
        ratio = 1.
        K = self.currentK
        newK = self.protomodel.K
        if K > -20. and newK < K:
            ratio = numpy.exp(.5*( newK - K))

        if ratio >= 1.:
            self.highlight ( "info", "K: %.3f -> %.3f: r=%.4f, take the step" % \
                             ( self.currentK, self.protomodel.K, ratio ) )
            if self.protomodel.K > 0. and self.protomodel.K < 0.7 * self.currentK:
                self.pprint ( " `- weirdly, though, K decreases. Please check." )
                sys.exit(-2)
            self.takeStep()
        else:
            u=random.uniform(0.,1.)
            if u > ratio:
                self.pprint ( "u=%.2f > %.2f; K: %.2f -> %.2f: revert." % (u,ratio,self.currentK,
                                self.protomodel.K) )
                self.manipulator.restoreModel( reportReversion=True )
                if hasattr ( self, "oldgrad" ) and self.accelerator != None:
                    self.accelerator.grad = self.oldgrad
            else:
                self.pprint ( "u=%.2f <= %.2f ; %.2f -> %.2f: take the step, even though old is better." % (u, ratio,self.currentK,self.protomodel.Z) )
                self.takeStep()

    def log ( self, *args ):
        """ logging to file """
        with open( "%s/walker%d.log" % ( self.rundir, self.walkerid ), "a" ) as f:
            f.write ( "[randomWalker-%s] %s\n" % ( time.strftime("%H:%M:%S"), " ".join(map(str,args)) ) )

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
                        self.pprint ( "extracted: %s" % point )
                    with open("%s/exceptions.log" % self.rundir,"a") as f:
                        f.write ( "%s: taking a step resulted in exception: %s, %s\n" % \
                                  (time.asctime(), type(e), e ) )
                        f.write ( "   `- exception occured in walker #%s\n" % \
                                  self.protomodel.walkerid )
                    sys.exit(-1)

            #If no combination was found, go back
            if self.protomodel.K is None:
                self.manipulator.restoreModel()
                continue

            # obtain the ratio of posteriors
            self.decideOnTakingStep ()
            self.record()
        self.manipulator.M.delCurrentSLHA()
        self.pprint ( "Was asked to stop after %d steps" % self.maxsteps )

if __name__ == "__main__":
    D = {'masses': {1000022: 299., 1000006: 1159.7, 1000001: 875.37}, 
         'ssmultipliers': {(1000022, 1000022): 0.06, (1000006, 1000006): 0.482, (-1000006, 1000006): 0.482, (-1000006, -1000006): 0.482, (1000006, 1000022): 0.482, (-1000006, 1000022): 0.482, (1000001, 1000001): 0.83, (-1000001, 1000001): 0.83, (-1000001, -1000001): 0.83, (1000001, 1000022): 0.83, (-1000001, 1000022): 0.83, (1000001, 1000006): 0.83, (-1000001, 1000006): 0.83, (-1000006, 1000001): 0.83, (-1000006, -1000001): 0.83}, 
         'decays': {1000022: {}, 1000006: {(1000022, 6): 1.0}, 1000001: {(1000022, 1): 1.0}}, 
         'K': 6.34, 'Z': 3.14, 'step': 16, 'timestamp': 'Wed May 24 23:18:37 2023', 'walkerid': 2}
    dbpath = "/users/wolfgan.waltenberger/git/smodels-database"
    dbpath = "official"
    walker = RandomWalker.fromDictionary ( D, walkerid = 0, dbpath = dbpath, 
                nevents = 5000, do_combine = True ) 
    walker.walk()
