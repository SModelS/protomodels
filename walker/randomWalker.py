#!/usr/bin/env python3

""" a first start at the random walk idea """

__all__ = [ "RandomWalker" ]

import signal
import pickle, sys, time, math, socket, os
import subprocess
import colorama
import numpy as np

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
from ptools import helpers
from ptools.helpers import prettyPrint

try:
    from smodels.statistics.pyhfInterface import setBackend
    o = setBackend ( "pytorch" )
    if not o:
        logger.warning ( "could not set backend to pytorch, falling back to numpy" )
except ImportError as e:
    logger.warning ( "could not set backend to pytorch: {e} (are you using smodels >=3.0.1)?" )

logger.setLevel("ERROR")

class RandomWalker ( LoggerBase ):
    def __init__ ( self, rvars : dict ):
#        walkerid : Union[str,int] = 0, nsteps : int = 10000,
#            strategy : str = "aggressive",
#            cheatcode : Union[str,int] = "no_cheat", dbpath : PathLike = "./database.pcl",
#            expected : bool = False, select : str = "all", cap_ssm = 100,
#            catch_exceptions : bool = True, rundir : Union[PathLike,None] = None,
#            do_srcombine : bool = False, test_param_space = False, run_mcmc = False,
#            record_history : bool = False, seed : Union[int,None] = None,
#            stopTeleportationAfter : int = -1,
#            templateSLHA : os.PathLike = "template_default.slha",
#            allowN1N1Prod : bool = False, susy_mode : bool = False,
#            use_initialiser : Union[str,bool] = False ):
        """ initialise the walker

        rvars ( dict ):
          - nsteps: maximum number of steps to perform, negative is infinity
          - cheatcode: cheat mode. 0 or "no_cheat" is no cheating, else
        cheatcode is path to model.
          - expected: remove possible signals from database
          - select: select only subset of results (all for all, em for
                efficiency maps only, ul for upper limits only, alternatively
                select for txnames via e.g. "txnames:T1,T2"
          - cap_ssm: set the maximum value for all signal strength multipliers (default=100)
          - catch_exceptions: should we catch exceptions
          - do_srcombine: if true, then also perform combinations, either via
                           simplified likelihoods or via pyhf
          - test_param_space: if true, run with constant K and TL (=1.0)
          - run_mcmc: if true, run mcmc walk without changing dimensions
          - record_history: if true, attach a history recorder class
          - seed: random seed, int or None
          - stopTeleportationAfter: int or None. we stop teleportation after
            this step nr.  If negative or None, we dont teleport at all
          - templateSLHA: the template file that is used
          - allowN1N1Prod: allow N1 N1 production mode
          - susy_mode: susy mode (penalty for ssms away from unity)
          - use_initialiser: if string, then interpret it as path to database
            dictionary file that we use for the initialiser. only works if
            not using cheatcode. if false, then dont use initialiser.
        """
        rvars = self.defaults ( rvars )
        globals().update ( rvars ) # doesnt work for all
        dbpath = rvars["dbpath"]
        use_initialiser = rvars["use_initialiser"]
        jmax = rvars["jmax"]
        continueFrom = rvars["continueFrom"]
        cheatcode = rvars["cheatcode"]
        do_srcombine = rvars["do_srcombine"]
        stopTeleportationAfter = rvars["stopTeleportationAfter"]
        run_mcmc = rvars["run_mcmc"]

        #call the super class of the random walker i.e Loggerbase
        super ( RandomWalker, self ).__init__ ( walkerid )
        dbpath = os.path.expanduser ( dbpath )
        if type(walkerid) != int or type(nsteps) != int or type(strategy)!= str:
            self.pprint ( f"Wrong call of constructor: {walkerid}, {nsteps}, {strategy}" )
            sys.exit(-2)
        self.walkerid = walkerid ## walker id, for parallel runs
        self.templateSLHA = templateSLHA
        self.rundir = rundir
        self.cap_ssm = cap_ssm
        if rundir == None:
            self.rundir = "./"
        self.random_seed = np.random.seed()
        self.test_param_space = test_param_space
        if seed is not None:
            self.random_seed = seed
            helpers.seedRandomNumbers(self.random_seed + walkerid )
            self.pprint ( f"setting random seed to {self.random_seed}" )
        helpers.mkdir ( "dictfiles" )
        self.dictfile = f"dictfiles/pmodel_{walkerid}.dict"
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
                dbversion = self.predictor.database.databaseVersion,
                templateSLHA = templateSLHA, allowN1N1Prod = allowN1N1Prod,
                susy_mode = susy_mode )

        self.manipulator = Manipulator ( protomodel, strategy,
                        do_record = record_history, seed = self.random_seed )
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
        # self.pprint ( f"Ramping up with slurm jobid {jobid} using template {templateSLHA} allowN1N1 {allowN1N1Prod} susy_mode {susy_mode}" )
        self.pprint ( f"Ramping up with slurm jobid {jobid}" )
        self.pprint ( f"It is {time.asctime()}" )
        self.pprint ( f"template {templateSLHA} allowN1N1 {allowN1N1Prod} susy_mode {susy_mode}" )
        
        #keep track of log llhd ratio
        self.trace_logllhdratio = []
        self.run_mcmc = run_mcmc
        self.use_initialiser = use_initialiser
        self.initialiser = None
        if self.use_initialiser not in [ False, None ]:
            if cheatcode not in [ "no_cheat", "", "none", None, 0 ]:
                logger.error ( f"use_initialiser {use_initialiser} specified, but also cheatcode {cheatcode} defined" )
                logger.error ( f"cheatcode takes precedence" )
            else:
                from walker.initialiser import Initialiser
                self.initialiser  = Initialiser ( self.walkerid, self.use_initialiser,
                       allowN1N1Prod = allowN1N1Prod, verbose = False,
                       dbpath = dbpath, templateName = templateSLHA )
                #init_model = self.initialiser.propose()
                init_model = self.initialiser.bestOfN(5)
                if init_model != None:
                    self.manipulator.initFromDict ( init_model )
        if self.run_mcmc: self.highlight("info", "Running MCMC walk")
        if cheatcode in [ "no_cheat", "", "none", None, 0 ]:
            self.takeStep() # the first step should be considered as "taken"
            #Set current TL and K values to threshold values
            self.currentTL = -0.1
            self.currentK = -20.0
        else:
            self.manipulator.cheat ( cheatcode )
            if self.test_param_space:
                self.manipulator.M.K, self.manipulator.M.TL = 1.0,1.0
                self.currentK = self.manipulator.M.K
                self.currentTL = self.manipulator.M.TL
            
            else:
                self.predict(self.manipulator)
                if type(self.manipulator.M.TL) != type(None) and type(self.manipulator.M.K) != type(None):
                    self.log ( f"Cheat model gets TL={self.manipulator.M.TL:.2f}, "\
                                  f"K={self.manipulator.M.K:.2f}" )
                    # self.printStats ( substep=4 )
                    self.manipulator.backupModel()
                    self.hiscoreList.newResult ( self.manipulator )
                self.printStats ( substep=5 )
                #self.manipulator.M.K = 1.0
                #self.manipulator.M.TL = 1.0
                self.currentK = self.manipulator.M.K
                self.currentTL = self.manipulator.M.TL
                if self.run_mcmc: self.currentBestCombo = set(self.manipulator.M.description.split(','))


    def setWalkerId ( self, Id ):
        self.walkerid = Id
        self.manipulator.setWalkerId ( Id )

    @classmethod
    def fromProtoModel( cls, protomodel : ProtoModel, args : Dict ):
        """ create a RandomWalker from a ProtoModel. Continue walking
            from that model """
        ret = cls( args )
        ret.manipulator.M = protomodel
        if "walkerid" in args:
            ret.manipulator.setWalkerId ( args["walkerid"] )
        ret.manipulator.backupModel()
        return ret

    def extractArguments ( func : Callable, args : Dict ) -> Dict:
        """ from args, extract all the entries that are parameters of func """
        pm = {}
        import inspect
        sig = inspect.signature ( func )
        for i in sig.parameters.keys():
            if i != "self" and i in args:
                pm[i]=args[i]
        return pm

    @classmethod
    def defaults ( cls, rvars ):
        """ define the defaults """
        defs = { "walkerid": "default", "use_initialiser": False, "jmin": 0, "continueFrom": "",
                 "cheatcode": "no_cheat", "nsteps": 1000, "dbpath": "./database.pcl",
                 "strategy": "aggressive", "rundir": "./", "cap_ssm": 100.,
                 "test_param_space": False, "seed": None, "expected": False,
                 "record_history": False, "catch_exceptions": True,
                 "stopTeleportationAfter": -1, "templateSLHA": "template_default.slha",
                 "allowN1N1Prod": False, "susy_mode": False, "use_initialiser": False,
                 "do_srcombine": True, "run_mcmc": False }
        for k,v in defs.items():
            if not k in rvars:
                rvars[k]=v
        if not "jmax" in rvars:
            rvars["jmax"]=rvars["jmin"]+1
        return rvars

    @classmethod
    def fromDictionary( cls, dictionary : Union[PathLike,Dict], rvars : Dict ):
        """ create a RandomWalker from a hiscore dictionary. Continue walking
            from the model in that dictionary
        :param dictionary: either a dictionary, or a string containing a dictionary,
        :param rvars: the rvars dictionary see constructor
        or the path to a dictionary
        """
        rvars = cls.defaults ( rvars )
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

        ret = cls( rvars ) ## simply pass on all the arguments

        pm = RandomWalker.extractArguments ( ProtoModel.__init__, rvars )
        ret.manipulator.M = ProtoModel( **pm )
        ret.manipulator.initFromDict ( dictionary )
        if "walkerid" in rvars:
            ret.manipulator.setWalkerId ( rvars["walkerid"] )
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
            self.log ( f"Step {self.protomodel.step} has {nUnfrozen}/{nTotal} unfrozen particles: {prtcles} [in best combo: {prtclesbc}]" )
            if len(pidsbc)>0 and not set(pidsbc).issubset ( set(pidsp) ):
                self.pprint ( f"  `-- error! best combo pids ({pidsbc}) arent subset of masses pids ({pidsp})!" )
                self.manipulator.M.bestCombo = None

    def writeToDictFile(self, proto_dict):
            # Load existing list or initialize a new one
            if os.path.exists(self.dictfile):
                with open (self.dictfile, "rt" ) as f:
                    txt = f.read()
                    try:
                        dicts = eval( txt )
                    except (SyntaxError,ValueError) as e:
                        logger.error ( f"when trying to read {self.dictfile}: {e}" )
                        dicts = []
                    
            else: dicts = []
            dicts.append(proto_dict)
            # helpers.mkdir ( os.path.dirname ( self.dictfile ) )
            with open (self.dictfile, "wt" ) as f:
                f.write (f"{dicts}")


    def predict ( self, manipulator : Manipulator) -> bool:
        """ Calls predictor.predict to get the theory predictions for model. 
        Loops for 5 times till model.muhat is close to 1.0 
        
        :returns: true if worked
        """
        #print(f"Adress of manip : {id(manipulator)}")
        model = manipulator.M
        if self.test_param_space:
            model.K = 1.0
            model.TL = 1.0
            proto_dict = manipulator.getPmodelDict()
            self.log(f"Protomodel: {proto_dict}")
            self.writeToDictFile(proto_dict)
            return True

        muhat_converge = False
        previousMuhat = None
        for i in range(5):
            predict = self.predictor.predict(manipulator, run_mcmc=self.run_mcmc)
            if predict: #returns False if no preds are found or TL is None (i.e no comb found)
                #print(f"i {i}, muhat {model.muhat}, convergence {abs(model.muhat - 1.0)}")
                if abs(model.muhat - 1.0) < 1e-02:
                    self.log(f"Step {model.step} converged at loop {i} with muhat {model.muhat}!")
                    muhat_converge = True
                    #proto_dict = manipulator.getPmodelDict()
                    #self.log(f"Protomodel: {proto_dict}")
                    break
                previousMuhat = model.muhat
                if model.muhat == 0.0: break
                manipulator.rescaleSignalBy(model.muhat, cap_ssm = self.cap_ssm) #?
            else:
                break # Rescale signal by a significant number?

        if not muhat_converge:  #reverting step
            model.K = None
            model.TL = None
            proto_dict = manipulator.getPmodelDict()
            if predict:
                self.log ( f"Step {model.step} did not converge to muhat 1.0, model muhat is {previousMuhat}. Going back to previous step." )
                self.log(f"Protomodel: {proto_dict}")   
                self.writeToDictFile(proto_dict)
            else:
                self.log ( f"Step {model.step} did not converge to muhat 1.0. Model did not find any prediction." )
                self.log(f"Protomodel: {proto_dict}") 
                self.writeToDictFile(proto_dict)
            return False

        return True

    def onestep ( self ):
        #Add one step
        self.protomodel.step+=1
        self.pprint (f"Step {self.protomodel.step} begins.")
        self.printStats( )
        #Remove data about best combo
        self.log("Clean best combo")
        self.protomodel.cleanBestCombo()
        # self.printStats( substep=11 )
        printMemUsage = False
        if printMemUsage:
            self.pprint ( f"memory footprint (kb): walker {asizeof(self)/1024}, model {asizeof(self.protomodel)/1024}" )

        #Trim the model, so we start only with the relevant particles for the
        #best combination in the previous step -> doing in predictor now
        # self.printStats( substep=12 )

        #Take a step in the model space:
        self.log("Randomly change model")
        self.manipulator.randomlyChangeModel(run_mcmc = self.run_mcmc, cap_ssm = self.cap_ssm)
        self.manipulator.reassignPID()
        # self.printStats( substep=13 )

        nUnfrozen = len( self.protomodel.unFrozenParticles() )
        ## number of pids in best combo, as a check

        #Try to create a simpler model
        #(merge pre-defined particles if their mass difference is below dm)
        if not self.run_mcmc:
            self.log("Try to simplify model")
            protomodelSimp = self.manipulator.simplifyModel(dm=200.0)
        else: protomodelSimp = None
        manipulatorSimp = None
        if protomodelSimp:
            manipulatorSimp = Manipulator ( protomodelSimp, strategy="aggressive",do_record = False, seed = self.random_seed )
            manipulatorSimp.reassignPID()
        boolProtoSimp = False

        # self.printStats( substep=14 )

        if self.catch_exceptions:
            try:
                if not self.predict(self.manipulator):
                    self.protomodel.K = None
                    return 
                if protomodelSimp:
                    boolProtoSimp = self.predict(manipulatorSimp) 
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
                    self.manipulator.M._writeSLHAFile ( f )
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
                        and (protomodelSimp.K >= self.manipulator.M.K)):
                self.log("Accepting the simplified model")
                self.manipulator.proposal_ratio['q_total'] *= self.manipulator.proposal_ratio['merge']['q']
                self.manipulator.M = protomodelSimp

        #If no combination could be found, return
        if self.manipulator.M.TL is None or self.manipulator.M.K is None:
            return

        if len(self.manipulator.M.rvalues) > 1:
            self.log ( f"Top r values are: {self.manipulator.M.rvalues[0]:.2f}, {self.manipulator.M.rvalues[1]:.2f}" )

        self.log ( f"Step {int(self.protomodel.step)}: found highest TL: {self.protomodel.TL:.2f}" )

        nUnfrozen = len ( self.protomodel.unFrozenParticles() )
        self.log ( f"Best combo is {self.protomodel.letters}: {self.protomodel.description}: [K={self.protomodel.K:.2f}, TL={self.protomodel.TL:.2f}, {int(nUnfrozen)} unfrozen]" )

        #For low scoring models, teleport to a high score model:
        if self.checkIfToTeleport( pmax=0.5, norm = 10.0 ):
            # if we teleport the rest becomes irrelevant
            return
        self.printStats( )

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
            self.log ( f"teleportation is turned off after step #{int(self.stopTeleportationAfter)}" )
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
        a = np.random.uniform ( 0., 1. )
        doTP = ( a < prob ) ## do teleport, yes or no
        sDoTP = "a>p: dont teleport."
        if doTP:
            sDoTP = "a<p: do teleport."
        self.log ( f"check if to teleport, Kmax={bestK:.2f}, ours is={ourK:.2f}, p={prob:.2f}, a={a:.2f}, {sDoTP}" )
        if doTP:
            self.manipulator.teleportToHiscore()
        return doTP

    def takeStep ( self ):
        """ take the step, save it as last step """
        if not self.test_param_space:
            ## possibly add to hiscore list
            self.log ( f"Step {self.protomodel.step} check if result goes into hiscore list" )
            #srs = ", ".join ( [ f"{x:.2f}" for x in self.protomodel.rvalues[:3] ] )    #protomodel.rvalues were used before to find the max allowed mu
            #self.log ( f"r values before calling .newResult are at {srs}" )
            self.hiscoreList.newResult ( self.manipulator ) ## add to high score list
            #srs = ", ".join ( [ f"{x:.2f}" for x in self.protomodel.rvalues[:3] ] )
            #self.log ( f"r values after calling .newResult are at {srs}" )
            self.log ( "done check for result to go into hiscore list" )
        ## Backup model
        self.manipulator.backupModel()
        # Update current K and TL values
        self.currentK = self.protomodel.K
        self.currentTL = self.protomodel.TL
        if self.run_mcmc: self.currentBestCombo = set(self.protomodel.description.split(','))
        self.manipulator.record( "take step" )

    def decideOnTakingStep ( self):
        """ depending on the ratio of K values, decide on whether to take the step or not.
            If ratio > 1., take the step, if < 1, let chance decide. """
        K = self.currentK
        log_llhdRatio_current = self.currentTL
        if K == None: # if the old is none, we do everything
            self.takeStep()
            return
        
        if self.run_mcmc:       #check if combination of results changes during mcmc walk
            newcombo = set(self.protomodel.description.split(','))
            if newcombo != self.currentBestCombo:
                self.log("Best Combination of results changed. Go back to previous model")
                proto_dict = self.manipulator.getPmodelDict(acc=False, critic_acc=False)
                self.log(f"Protomodel: {proto_dict}")
                self.manipulator.restoreModel( reportReversion=True )
                return
        newK = self.protomodel.K
        log_llhdRatio_new = self.protomodel.TL
        
        if self.test_param_space:
            self.takeStep()
            self.log("Testing parameter space, K,TL = 1.0. Take step")
            return
        
        if newK == None:
            # if the new is none, but the old isnt, we go back
            self.manipulator.restoreModel( reportReversion=True )
            return
        if log_llhdRatio_current >= 1.0 and log_llhdRatio_new < 1.0:
            #If current log likeihood ratio >= 1.0, but the new step has a log llhd ratio < 1.0, return to previous protomodel
            self.log(f"Previous TL: {log_llhdRatio_current} >= 1.0, New TL: {log_llhdRatio_new} < 1.0. Going back to previous model")
            proto_dict = self.manipulator.getPmodelDict()
            self.log(f"Protomodel: {proto_dict}")
            self.writeToDictFile(proto_dict)
            self.manipulator.restoreModel( reportReversion=True )
            return
        #K = 2 log (L1/L0) + 2 log(prior)
        #acceptance ratio = (new_L1/new_l0)/(current_L1/current_L0) * (new_prior)/(old_prior) * proposal_ratio
        #acceptance ratio = exp( (log(new_L1/new_l0) + log(new_prior)) - (log(current_L1/current_L0) - log(old_prior)) + log(proposal_ratio))
        # acceptance ratio = exp(1/2 (newK - K) + log (proposal_ratio))
        log_prop_ratio = np.log(self.manipulator.proposal_ratio['q_total'])
        acceptance_ratio = np.exp(0.5*( newK - K) + log_prop_ratio)
        self.log(f"Step {self.protomodel.step}: Acceptance ratio: {acceptance_ratio}, K ratio: {0.5*(newK-K)}, Log proposal ratio: {log_prop_ratio}")
        #print(f"Step {self.protomodel.step}: Acceptance ratio {acceptance_ratio}, K {K}, newK {newK}")
        if acceptance_ratio > 1:
            self.highlight ( "info", f"Acceptance ratio > 1.0. K: {prettyPrint(K)} -> {prettyPrint(newK)}; Check Critics." )

            cr, cr_answer = self.critic.predict_critic(self.protomodel, keep_predictions=True)
            if cr:
                self.highlight ( "info", "Passed both critics, taking the step." )
                proto_dict = self.manipulator.getPmodelDict(acc=True, critic_acc=True)
                self.log(f"Protomodel: {proto_dict}")
                self.writeToDictFile(proto_dict)
                self.takeStep()
            else:
                self.highlight ( "info", f"Failed critic: {cr_answer}; the step is reverted." )
                proto_dict = self.manipulator.getPmodelDict(acc=True, critic_acc=False)
                self.log(f"Protomodel: {proto_dict}")
                self.writeToDictFile(proto_dict)
                self.manipulator.restoreModel( reportReversion=True )

        else:
            #Draw random number u
            from scipy.stats import uniform
            u = uniform.rvs(loc=0.,scale=1.,size=1)[0]
            #print(f"Step {self.protomodel.step}: u {u}, Acceptance ratio {acceptance_ratio}, K {K}, newK {newK}")
            if u > acceptance_ratio:
                self.highlight ("info", f"u={u:.2f} > {acceptance_ratio:.2f}; K: {prettyPrint(K)} -> {prettyPrint(newK)}: revert." )
                proto_dict = self.manipulator.getPmodelDict(acc=False, critic_acc=False)
                self.log(f"Protomodel: {proto_dict}")
                self.writeToDictFile(proto_dict)
                self.manipulator.restoreModel( reportReversion=True )
            else:
                self.highlight ( "info", f"u={u:.2f} <= {acceptance_ratio:.2f};K: {prettyPrint(K)} -> {prettyPrint(newK)}; Check Critics." )   #SN: <+ and not > right?

                cr, _ = self.critic.predict_critic(self.protomodel, keep_predictions=True)
                if cr:
                    self.log ( "Passed both critics, taking the step." )
                    proto_dict = self.manipulator.getPmodelDict(acc=True, critic_acc=True)
                    self.log(f"Protomodel: {proto_dict}")
                    self.writeToDictFile(proto_dict)
                    self.trace_logllhdratio.append(log_llhdRatio_new - log_llhdRatio_current)
                    self.takeStep()
                else:
                    self.log ( "Failed at least one critic, the step is reverted." )
                    proto_dict = self.manipulator.getPmodelDict(acc=True, critic_acc=False)
                    self.log(f"Protomodel: {proto_dict}")
                    self.writeToDictFile(proto_dict)
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
            self.log("Only LSP present. Forcing to unfreeze random particle")
            self.manipulator.propose_model = self.manipulator.M.copy()
            self.manipulator.proposal_ratio = {'add_par':{'q':1.0}, 'rem_par':{'q':1.0}, 'br':{'q':1.0}, 'ssm':{'q':1.0}, 'q_total':1.0}
            
            unfrozenParticle = self.manipulator.randomlyUnfreezeParticle()
            self.manipulator.run_mcmc = False
            self.manipulator.proposal_density( move='add_par', force_move=True)
            self.manipulator.backupModel()
        
        old_handler = signal.getsignal(signal.SIGTERM)
        
        def handle_termination(signum=None, frame=None):
            """Handles both SLURM termination signals and manual interruptions."""
            if signum in [ signal.SIGTERM, signal.SIGINT ]:
                self.highlight("info", f"Saving current protomodel to pmodel{self.walkerid}.dict")
                self.manipulator.restoreModel( reportReversion=True )
                helpers.mkdir ( "Pmodels" )
                self.manipulator.writeDictFile(outfile=f"Pmodels/pmodel{self.walkerid}.dict", step=self.manipulator.M.step - 1)
            if old_handler is signal.SIG_DFL:
                # Default behavior for SIGINT is to raise KeyboardInterrupt,
                # which usually exits with code 130.
                print("Exiting gracefully...")
                sys.exit(130)
            elif old_handler is signal.SIG_IGN:
                print("Old handler ignored SIGINT, continuing.")
            else:
                # Call previous custom handler
                old_handler(signum, frame)
            # os.kill ( os.getpid(), signal.SIGKILL )
            # sys.exit(0)
        # Register signal handlers for graceful shutdown
        signal.signal(signal.SIGTERM, handle_termination)  # SLURM termination signal
        # signal.signal(signal.SIGKILL, handle_termination)  # OOM memory error signal
        signal.signal(signal.SIGINT, handle_termination)   # Manual interruption (Ctrl+C)
        signal.signal(signal.SIGUSR1, handle_termination)  # SLURM preemption signal

        while self.maxsteps < 0 or self.protomodel.step<self.maxsteps:
            if not catchem:
                try:
                    self.onestep()
                    if self.manipulator.M.step % 1000 == 0: # log every 1000th step
                        self.manipulator.writeDictFile(outfile=f"Pmodels/pmodel{self.walkerid}.dict", step=self.manipulator.M.step )
                except KeyboardInterrupt:
                    self.highlight(f"error", "Interrupted by user/ Killed by slurm.")
                    handle_termination()
                    
            else:
                try:
                    self.onestep()
                except Exception as e:
                    # https://bioinfoexpert.com/2016/01/18/tracing-exceptions-in-multiprocessing-in-python/
                    self.pprint ( f"taking a step resulted in exception: {type(e)}, {e}" )
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
                    self.highlight("info", f"Saving last protomodel to pmodel{self.walkerid}.dict")
                    self.manipulator.restoreModel( reportReversion=True )
                    helpers.mkdir ( "Pmodels" )
                    self.manipulator.writeDictFile(outfile=f"Pmodels/pmodel{self.walkerid}.dict", step=self.manipulator.M.step - 1)
                    sys.exit(-1)

            #If no combination was found, go back
            if self.protomodel.K is None:
                self.log("K is none, return to previous model")
                self.manipulator.restoreModel(reportReversion=True)
                continue

            # obtain the ratio of posteriors
            self.decideOnTakingStep ()
            self.record()
            smaxstp = f"{self.maxsteps}"
            if self.maxsteps < 0:
                smaxstp = "inf"
            self.pprint ( f"Step {self.protomodel.step}/{smaxstp} finished." )
        
        self.manipulator.M.delCurrentSLHA()
        self.pprint ( f"Was asked to stop after {self.maxsteps} steps" )
        self.pprint(f"Writing last protomodel to pmodel{self.walkerid}.dict")
        helpers.mkdir ( "Pmodels" )
        self.manipulator.writeDictFile(outfile=f"Pmodels/pmodel{self.walkerid}.dict")

def model1():
    masses = {1000022: 46.514732, 1000023: 350, 1000024: 350 }
    
    ssms = {(1000022, 1000022): 0.191043, (1000023, 1000023): 0.204701, 
        (1000022, 1000023): 0.0, (1000024, 1000024): 1.829998, 
        (-1000024, 1000024): 1.829998, (-1000024, -1000024): 1.829998, 
        (1000022, 1000024): 1.329093, (-1000024, 1000022): 1.329093, 
        (1000023, 1000024): 0.0, (-1000024, 1000023): 0.0 }
    
    decays = {1000022: {}, 1000023: {(1000022, 25): 1.0}, 1000024: {(1000022, 24): 1.0}, 1000037: {(1000022, 24): 1.0}}
    
    D = {'masses': masses, 'ssmultipliers': ssms, 'decays': decays }
    return D

def model2():
    D={
        'masses': {
            1000022: 35.00333930498704,
            1000023: 67.46989452437126,
            1000024: 86.39134410603818,
            1000006: 466.77198888781453
        },
        'ssmultipliers': {
            (1000022, 1000022): 3.076810727833559,
            (1000022, 1000023): 0.45615298790558134,
            (1000022, 1000024): 0.2420703888833409,
            (1000023, 1000024): 13.70830463164829,
            (1000023, 1000023): 0.2420703888833409,
            (-1000006, 1000006): 0.7140157447902225
        },
        'decays': {
            1000022: {},
            1000023: {},
            1000006: {
                (1000022, 6): 0.17524601220007244,
                (1000023, 6): 0.35964987901179013,
                (1000024, 5): 0.4651041087881375
            }
        },
        'templateSLHA': 'naturalEwkino_nodegeneracy.slha',
        'allowN1N1Prod': True,
        'susy_mode': True,
    }
    return D

if __name__ == "__main__":
    dbpath = "../../smodels-database/"
    dbpath = "official"
    select = "txnames:electroweakinos,electroweakinos_offshell"
    select = "all"
    D = model2()
    #walker = RandomWalker( walkerid=0, nsteps = 1000,
    #                dbpath=dbpath, cheatcode=1, select=select, do_srcombine = True )
    rvars = { "walkerid": 0, "dbpath": dbpath, "do_srcombine": True,
              "select": select, "templateSLHA": "templateNaturalEwkino.slha",
              "allowN1N1Prod": True, "susy_mode": False, "use_initialiser": "./ini.dict" }
               
    walker = RandomWalker.fromDictionary ( D, rvars )
    walker.walk()
