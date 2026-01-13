#!/usr/bin/env python3

""" Class that encapsulates the manipulations we perform on the protomodels,
    so that the protomodel class is a data-centric class, and this one
    an algorithm-centric class.
"""

__all__ = [ "Manipulator" ]

import copy, os, sys, itertools, colorama, random
import numpy as np
from typing import Union, Dict, List, Tuple, Set
from unum import Unum
from os import PathLike
import tempfile
from scipy.stats import norm, lognorm, uniform

from smodels.base.physicsUnits import fb, TeV, GeV
from smodels.base.crossSection import LO
from smodels.matching.theoryPrediction import TheoryPrediction
import smodels

from smodels_utils.helper.terminalcolors import *

from base.loggerbase import LoggerBase
from base.constants import smMasses, smWidths
from base.runEnviron import RunEnviron
from builder.protomodel import ProtoModel
from ptools.sparticleNames import SParticleNames
from ptools.helpers import nround, getAllPidsOfTheoryPred, py_dumps, mkdir, \
         formatObject

LSP = ProtoModel.LSP

class Manipulator ( LoggerBase ):
    """ contains the protomodel manipulation algorithms. """

    # walledpids are particle ids that have a minimum mass requirement,
    # given as the values in the dictionary
    walledpids = { 1000001 : 310, 1000002 : 310, 1000003 : 310, 1000004 : 310,
                   1000021 : 310, 1000023 : 20, 1000024 :  100, 1000037 : 100,
                   1000025 : 20 }
    ## forbiddenparticles are particle ids that we do not touch in this run
    # forbiddenparticles = []

    # decayless particles: list of particles that are allowed to not
    # have any branching ratio parameters defined. LSP, because it is the DMC,
    # N2 and C1 because they may have natural offshell boson decays defined,
    # ( so they dont have any free parameters)
    # This happens eg when using templateNaturalEwkino.slha
    # FIXME might have to make this smarter later on

    def __init__ ( self, protomodel : Union[ProtoModel,Dict,PathLike],
            environ : RunEnviron, verbose : bool = False,
            do_record : bool = False, seed : Union[bool,int] = None,
            nth : int = 0, walkerid : Union[None,int] = None,
            initTestStats : bool = False ):
        """
        :param protomodel: is either a protomodel, or a hiscore dictionary,
        or a path to a protomodel
        :param strategy: what combination strategy?
        currently we only have "aggressive"
        :param do_record: do record actions taken
        :param seed: random seed
        :param nth: if initialisation from hiscores dict file, initialise
        from nth entry in that dict file
        :param walkerid: usually taken from protomodel (None), but can
        also specify
        :param initTestStats: if True, set also test statistics K and TL

        Example usage:

        .. code-block:: python3

            >>> # instantiate a model plus manipulator from hiscore file
            >>> m = Manipulator( "hiscores.dict" )
            >>> # get the predictions
            >>> predictor = Predictor(0, do_srcombine=True )
            >>> predictor.predict ( m.M, keep_predictions=True )
            >>> # print test statistics,
            >>> # best combination, most constraining analyses, etc
            >>> m.describe()
        """
        super(Manipulator, self).__init__ ( walkerid if walkerid is not None else 0 )
        if type ( protomodel ) == ProtoModel:
            ## make sure we log correctly asap
            self.walkerid = protomodel.walkerid
        self.namer = SParticleNames ( False )
        self.run_mcmc = False
        self.M = protomodel
        self.propose_model = None

        if type(protomodel) == dict:
            self.M = ProtoModel ( walkerid = walkerid, environ = environ )
            self.initFromDict ( protomodel, initTestStats = initTestStats )
        if type(protomodel) == str:
            self.M = ProtoModel ( walkerid = walkerid, environ = environ )
            if protomodel.endswith ( ".dict" ) or protomodel.endswith ( ".truth" ):
                self.initFromDictFile ( protomodel, nth = nth )
        self.walkerid = self.M.walkerid
        self.seed = seed
        self.verbose = verbose

        #need the below attribute to be defined here, in case of using hiscoreCLI to manipulate the protomodel
        self.proposal_ratio = {'add_par':{'q':1.0}, 'rem_par':{'q':1.0}, 'br':{'q':1.0}, 'ssm':{'q':1.0}, 'q_total':1.0}

        #Store a canonical order for the masses. So the ordering in each tuple is enforced
        #consider Xu and Xs as an extension of Xd, Xc not considered since charm-tagging results exist
        self.canonicalOrder =  [( 1000001, 1000002 ),(1000001, 1000003),( 1000005, 2000005 ),( 1000006, 2000006 ),
                          ( 1000023, 1000025 ), ( 1000024, 1000037 ),( 1000025, 1000035 )]

        #Define groups of particles to be merged if their mass difference is below a given mass gap
        self.mergerCandidates =   [(1000001, 1000002, 1000003), ( 1000005, 2000005 ), (1000006, 2000006),
                                    ( 1000024, 1000037 ), ( 1000023, 1000025 )]
        self.do_record = do_record ## if True, then record changes
        self.recording = [] ## just in case

    def forcedMassDegeneratePids ( self, pid : int,
           include_self : bool = True ) -> List:
        """ determine which pids are forced mass degenerate to ours

        :returns: list of all mass degenerate pids
        """
        ret = [pid]
        for pids in self.M.forced_degeneracies:
            if pid in pids:
                ret = list(pids)
        if not include_self:
            ret.remove(pid)
        return ret

    def getClosestPair ( self, pids ):
        """ of <n> PIDs, identify the two that are closest in mass """
        if len(pids)<2:
            return None
        dmin = float("inf")
        pair = (0,0)
        for pid1 in pids:
            for pid2 in pids:
                if pid1 == pid2:
                    continue
                dm = abs ( self.M.masses[pid2] - self.M.masses[pid1] )
                if dm < dmin:
                    dmin = dm
                    pair = ( pid1, pid2 )
        return pair,dmin

    def checkIfOffshell(self, pid, protomodel=None):
        if protomodel is None:
            protomodel = self.M
        offshell = False
        if 1000023 in protomodel.unFrozenParticles() or 1000024 in protomodel.unFrozenParticles():
            if pid == 1000023 and (protomodel.masses[pid] - protomodel.masses[LSP]) < (smMasses["Z"] + smWidths["Z"]): offshell = True
            elif pid == 1000024 and (protomodel.masses[pid] - protomodel.masses[LSP]) < (smMasses["W"] + smWidths["W"]): offshell = True
            else: offshell = False
        else: offshell = False

        return offshell

    def teleportToHiscore ( self ):
        """ without further ado, discard your current model and start
            fresh with the hiscore model. """
        ## FIXME this is currently not used.
        fname = "hiscores.dict"
        if not os.path.exists ( fname ):
            self.pprint ( "no other walkers found (could not find hiscores.dict)" )
            return
        try:
            with open ( fname, "rt" ) as f:
                dicts = eval ( f.read() )
        except (EOFError,SyntaxError,NameError) as e:
            self.pprint ( f"when trying to teleport, found {e}. cancel teleportation." )
            # can happen if it is just being written. in this case dont teleport
            return
        ith = 0
        choices = []
        f = 1
        for i in range(len(dicts)-1,-1,-1):
            choices += [i]*f
            f=f*2
        ith = int(np.random.choice ( choices ))
        self.log ( f"teleporting, we have {len(dicts)} dicts" )
        self.log ( f"choosing the {int(ith)}th entry, it has a K of {dicts[ith]['K']:.2f}" )
        step = self.M.step
        nth = f"{int(ith)}th"
        longforms = { 1: "first", 2: "second", 3: "third" }
        if ith in longforms:
            nth = longforms[ith]
        self.record ( f"teleporting to {nth} hiscore" )
        self.initFromDict ( dicts[ith], initTestStats=True )
        self.M.step = step ## continue counting!
        self.M.bestCombo = None

    def shiftAllMassesBy ( self, dm : float, lsp : bool = False ) -> list:
        """ shift the masses of all unfrozen particles by dm [GeV]
        (simple covenience function)

        :param dm: the shift of all masses, in GeV
        :param lsp: if true, then shift also the lsp

        :returns: list of all pids that got shifted
        """
        self.log(f"Shifting all other masses by {dm:.2f}")
        shifted = set()
        for pid,m in self.M.masses.items():
            if pid is LSP and not lsp:
                continue
            self.M.masses[pid]=m+dm
            shifted.add ( pid )
            if self.M.masses[pid] > self.M.maxMass:
                self.M.masses[pid] = self.M.maxMass
        return shifted

    def changeMassAccToSSM(self, pid_pair: Tuple, new_ssm: Union[int, float], sqrts:Union[int,float,Unum] = 13):
        """
        Change mass of a protomodel pid_pair according to an input signal strength multiplier keeping protomodel xsec
        at input sqrts the same. Update the protomodel pid_pair's mass and ssm
        :param pid_pair: Tuple of pids whose mass and ssm needs to be changed. Ex:(-1000006,1000006)
        :param new_ssm: The input signal strength multipier
        :param sqrts: The sqrts at which the protomodel xsec needs to remain the same. Default is 13(TeV)
        """

        model = self.M
        model_xsecs = model.getXsecs()[0]   #get all the protomodel xsec
        model_xs = 0.
        model_mass = 0.
        change = False

        if type(sqrts) == Unum: sqrts = sqrts.asNumber(TeV)
        #print("here")
        for xsec in model_xsecs:
            if (xsec.pid == pid_pair) and abs(xsec.info.sqrts.asNumber(TeV) - sqrts) <= 1e-07:  #get the protomodel xsec at pid_pair and sqrts
                model_xs = xsec.value.asNumber(fb)
                change = True

        if not change: return


        #print("Model Xsec: ", model_xs)
        new_susy_xs = (model_xs/new_ssm)*fb         #the new susy xsec according to the new_ssm
        #print("Susy Xsec: ", new_susy_xs)

        from ptools.xsecFit import XSecFitter
        func = XSecFitter(pid_pair, sqrts)

        model_mass = func.getValueFromFit(new_susy_xs, inverse=True)        #get the mass according to the new_susy_xsec from the xsecFitter
        if model_mass is None:
            print(f"No ref xsec available for {pid_pair}")
            return
        #print("new mass: ", model_mass)

        #update the protomodel's ssm and mass for the pid_pair
        model.ssmultipliers[pid_pair] = new_ssm
        if pid_pair[0] in model.masses: model.masses[pid_pair[0]] = model_mass.asNumber(GeV)
        else: model.masses[pid_pair[1]] = model_mass.asNumber(GeV)

    def changeSSMAccToMass(self, pid_pair: Tuple, new_mass: Union[int, float, Unum], sqrts:Union[int,float,Unum] = 13):
        """
        Change signal strength multiplier of a protomodel pid_pair according to an input mass keeping the protomodel xsec
        at input sqrts the same. Update the protomodel pid_pair's mass and ssm
        :param pid_pair: Tuple of pids whose mass and ssm needs to be changed. Ex:(-1000006,1000006)
        :param new_mass: The input mass
        :param sqrts: The sqrts at which the protomodel xsec needs to remain the same
        """

        model = self.M
        model_mass = new_mass

        if type(new_mass) is not Unum: model_mass = new_mass*GeV
        if type(sqrts) is Unum: sqrts = sqrts.asNumber(TeV)

        from ptools.xsecFit import xsecFitter
        func = xsecFitter(pid_pair, sqrts)

        susy_xs = func.getValueFromFit(model_mass, inverse=False)         #get the susy_xsec according to the new_mass from the xsecFitter
        if susy_xs is None:
            #print(f"No ref xsec for {pid_pair} at mass {model_mass}")
            return

        #print("new susy xsec: ", susy_xs)

        susy_xs = susy_xs.asNumber(fb)
        model_xsecs = model.getXsecs()[0]                                 #get all the protomodel xsec
        for xsec in model_xsecs:
            if (xsec.pid == pid_pair) and abs(xsec.info.sqrts.asNumber(TeV) - sqrts) <= 1e-07:  #get the protomodel xsec at pid_pair and sqrts
                model_xs = xsec.value.asNumber(fb)
            else: return

        new_ssm = (model_xs/susy_xs)        #the new_ssm according to the new_susy_xsec
        #print("New SSM: ", new_ssm)

        #update the protomodel's ssm and mass for the pid_pair
        model.ssmultipliers[pid_pair] = new_ssm
        if pid_pair[0] in model.masses: model.masses[pid_pair[0]] = model_mass.asNumber(GeV)
        else: model.masses[pid_pair[1]] = model_mass.asNumber(GeV)


    def getPmodelDict (self, get_xsecs : bool = False, acc : bool = False,
            critic_acc : bool = False) -> Dict:
        """ get the Pmodel dictionary
        :param critic_all: if true, passed the critic
        :param acc: if true, passed the acceptance ratio check
        """
        if type(self.M) == type(None):
            ## there is nothing to write
            self.log("No protomodel")
            return

        proto_dict = self.M.dict(sort_dict=True)
        if not get_xsecs and 'xsecs[fb]' in proto_dict.keys():
            del proto_dict['xsecs[fb]']

        addTexts = False

        if acc and critic_acc:
            proto_dict['Accepted'] = 0
            if addTexts:
                proto_dict['Accepted_Text'] = "passed all"
        elif acc and not critic_acc:
            proto_dict['Accepted'] = 1
            if addTexts:
                proto_dict['Accepted_Text'] = "passed acceptance ratio not critic"
        else:
            proto_dict['Accepted'] = 2
            if addTexts:
                proto_dict['Accepted_Text'] = "failed acceptance ratio"

        proto_dict['K'] = self.M.K
        proto_dict['TL'] = self.M.TL
        proto_dict['step'] = self.M.step

        return proto_dict

    @classmethod
    def writeDictionaryToFile ( cls, filename : os.PathLike,
           obj : dict, appendMode : bool = False ) -> dict:
        """ class method, write the dictionary obj in a formatted manner to file
        filename

        :param appendMode: if true, assume it is one of a list of dictionaries.
        intent.
        :returns: dictionary that was written out
        """


        mode,comma = "wt",""
        if appendMode:
            mode,comma = "at",","
        level = 1 if appendMode else 0
        # mkdir ( os.path.dirname ( filename ) )
        with open ( filename, mode ) as f:
            d = py_dumps ( obj, level = level )
            if appendMode:
                d = " "*4 + d
            f.write ( f"{d}{comma}\n" )
            #if not appendMode:
            #    f.write ( "\n" )
            f.close()
        return d

    def writeDictFile ( self, outfile : Union[str,None] = "pmodel.dict", step=None,
            cleanOut : bool = True, comment : str = "", appendMode : bool = False,
            ndecimals : int = 6 ) -> Dict:
        """ write out the dict file to outfile

        :param outfile: output file, but replacing %t with int(time.time()). If None,
                        then dont write file, just create dictionary object
        :param step: The step at which the protomodel is present. If None, step=self.M.step
        :param cleanOut: clean the dictionary from defaults, remove meta info
        :param comment: add a comment field
        :param ndecimals: number of digits after decimal when rounding
        :param appendMode: if true, append to file, and add comma after dictionary.
                           if false, overwrite, and no comma at the end.
        :returns: the dictionary with the protomodel
        """
        if type(self.M) == type(None):
            ## there is nothing to write
            return {}
        D = copy.deepcopy ( self.M.dict() )
        frozen = self.M.frozenParticles()
        if cleanOut:
            origMasses = self.M.dict()["masses"]
            ## but with a bit of cleaning!
            for k,v in origMasses.items():
                if v > 5e5:
                    D["masses"].pop(k)
                    if k in D["decays"]:
                        D["decays"].pop(k)
                else:
                    D["masses"][k]=round(v,ndecimals)
            for k,decays in self.M.dict()["decays"].items():
                for i,v in decays.items():
                    if not k in D["decays"]:
                        continue
                    if not i  in D["decays"][k]:
                        continue
                    if v < 1e-7:
                        D["decays"][k].pop(i)
                    else:
                        D["decays"][k][i]=round(v,ndecimals)
            for k,v in self.M.dict()["ssmultipliers"].items():
                ## if any of the pids is frozen, we dont write out
                hasFrozenPid = False
                for pid in k:
                    if abs(pid) in frozen:
                        hasFrozenPid = True
                if hasFrozenPid: #  or abs ( v - 1.) < 1e-5:
                    D["ssmultipliers"].pop(k)
                else:
                    D["ssmultipliers"][k]=round(v,ndecimals)
        if hasattr ( self, "seed" ) and self.seed != None:
            D["seed"]=self.seed
        D["TL"]=nround(self.M.TL,ndecimals)
        D["K"]=nround(self.M.K,ndecimals)
        D["walkerid"]=self.M.walkerid
        if step: D["step"]= step
        else: D["step"] =  self.M.step
        if not cleanOut:
            import time
            D["timestamp"]=time.asctime()
            D["walkerid"]=self.M.walkerid
            D["step"]=self.M.step
            D["protomodels_version"]=self.M.protomodels_version
            D["smodels_version"]=smodels.installation.version()
            D["database_version"]=self.M.dbversion
            D["templateSLHA"]=self.M.environ.templateName
            D["allowN1N1Prod"]=self.M.environ.allowN1N1Prod
            D["susy_mode"]=self.M.environ.susy_mode
        D["description"]=self.M.description
        if hasattr ( self.M, "ul_critic" ):
            D["ul_critic"]=self.M.ul_critic
        if hasattr ( self.M, "llhd_critic" ):
            D["llhd_critic"]=self.M.llhd_critic
        if hasattr ( self.M, "description" ) and D["description"]=="":
            D.pop ( "description" )
        if len(comment)>0:
            D["comment"]=comment
        if outfile == None:
            return D
        import time
        fname = outfile.replace("%t", str(int(time.time())) )
        self.writeDictionaryToFile ( fname, D, appendMode )
        return D

    def pidInList ( self, pid, lst, signed ):
        """ is pid in lst """
        if signed:
            return pid in lst
        return pid in lst or -pid in lst

    def initFromDictFile ( self, filename : PathLike, initTestStats : bool = False,
           nth : int = 0 ) -> bool:
        """ setup the protomodel from dictionary in file <filename>.
            If it is a list of dictionaries, take the 1st entry.
        :param filename: name of file
        :param initTestStats: if True, set also test statistics K and TL
        :param nth: if we find a list of models, pick the nth. 0 = 1st. If nth
                    does not exist, return False
        :returns: true, if successful
        """
        if not os.path.exists ( filename ):
            self.pprint ( f"filename {filename} does not exist!" )
            return False
        with open ( filename, "rt" ) as f:
            txt=f.read()
            txt = txt.replace( "null", "float('nan')" ) # json has 'null'
            D = eval ( txt )
        if type(D) == list:
            if len(D)<nth+1:
                self.pprint ( f"asking for {nth}th entry, but we only have {len(D)}" )
                return False
            if type(D[nth]) != dict:
                self.pprint ( f"{nth}th entry in list is not a dictionary" )
                return False
            self.initFromDict ( D[nth], filename, initTestStats )
            return
        if type(D) == dict:
            self.initFromDict ( D, filename )
            return True
        self.pprint ( f"dont understand content of file {filename}" )
        return False

    def diff ( self, other ):
        """ diff between our protomodel and <other> protomodel
        :returns: dictionary of differences
        """
        ret = {}
        keys = set ( self.M.__dict__.keys() )
        keys = keys.union ( set ( other.__dict__.keys() ) )
        for key in keys:
            if not key in self.M.__dict__:
                ret[key]="missing in self"
                continue
            if not key in other.__dict__:
                ret[key]="missing in other"
                continue
            selfV, otherV = getattr ( self.M, key ), getattr ( other, key )
            if type(selfV) != type(otherV):
                ret[key]=f"different types {type(selfV)} != {type(otherV)}"
                continue
            if type(selfV) in [ int, str ]:
                if selfV != otherV:
                    ret[key]=f"self is {selfV} other is {otherV}"
                continue
            if type(selfV) in [ float ]:
                if selfV != otherV:
                    ret[key]=f"self is {selfV} other is {otherV}"
                continue
            if type(selfV) in [ dict ]:
                if selfV != otherV:
                    ret[key]=f"dictionaries differ"
                continue
            if type(selfV) in [ list, set ]:
                if selfV != otherV:
                    ret[key]=f"containers differ"
                continue

        return ret

    #def initModel(self):

    #    self.propose_model = self.manipulator.M.copy()
    #    self.manipulator.proposal_ratio = {'add_par':{'q':1.0}, 'rem_par':{'q':1.0}, 'br':{'q':1.0}, 'ssm':{'q':1.0}, 'q_total':1.0}

    #    unfrozenParticle = self.manipulator.randomlyUnfreezeParticle()

    #    self.manipulator.proposal_density( move='add_par', force_unfreeze=True)
    #    self.manipulator.backupModel()

    def initFromDict ( self, D : Dict, filename : str = "",
            initTestStats : bool = False ):
        """ setup the protomodel from dictionary D.
        :param D: dictionary, as defined in pmodel*.dict files.
        :param filename: name of origin. not necessary, only for logging.
        :param initTestStats: if True, set also test statistics K and TL
        """
        scom = ""
        if "comment" in D:
                scom = f": {D['comment']}"
        if filename == "":
            line = "initializing from dictionary: "
            for k,v in D["masses"].items():
                line += f"{self.namer.asciiName(k)}, "
            line = line[:-2]
            self.pprint ( line )
        else:
            sK=""
            if "K" in D:
                sK = f" K={formatObject ( D['K'], 1 )}, "
            self.highlight ( "info", f"starting with {sK}{filename}{scom}" )
        if self.walkerid != None:
            self.M.walkerid = self.walkerid
        #Reset all model attributes:
        self.M.initializeModel()
        #Set attributes to dictionary values:
        for k,v in D["masses"].items():
            self.M.masses[k]=v
        for k,v in D["ssmultipliers"].items():
            self.M.ssmultipliers[k]=v
        for mpid,decays in D["decays"].items():
            if not mpid in self.M.decays:
                self.M.decays[mpid]={}
            for dpid,v in decays.items():
                self.M.decays[mpid][dpid]=v
        if "step" in D: ## keep track of number of steps
            self.M.step = D["step"]
        #if "walkerid" in D:
        #    self.M.walkerid = D["walkerid"]
        if initTestStats:
            if "TL" in D:
                self.M.TL = D["TL"]
            if "K" in D:
                self.M.K = D["K"]
            if "database_version" in D:
                self.M.dbver = D["database_version"]
        if "xsecs[fb]" in D:
            tmp = D["xsecs[fb]"]
            xsecs = []
            from smodels.base.crossSection import XSection
            for ss,value in tmp.items():
                xsec = XSection()
                xsec.value = value*fb
                xsec.info.sqrts = ss[1]*TeV
                xsec.info.label = f"{ss[1]} TeV (LO) [from_dict_file]"
                # xsec.info.label = "from_dict_file"
                xsec.info.order = 0
                xsec._pid = ss[0]
                xsecs.append ( xsec )
            self.M._xsecMasses = copy.deepcopy ( self.M.masses )
            self.M._xsecSSMs = copy.deepcopy ( self.M.ssmultipliers )
            self.M._stored_xsecs = ( xsecs, "loaded from dict file" )
            # self.M.computeXSecs()

    def cheat ( self, mode : Union[int,str] = "no_cheat" ):
        """ cheating, i.e. starting with models that are known to work well
        :param mode: if string, then this is path to cheat model file.
        if integer, then cheat model is f'Pmodels/pmodel{mode}.dict'
        """

        if mode in [ "no_cheat", "", None, 0 ]: ## no cheating
            return
        filename = mode
        if type(mode)==int:
            filename = f"Pmodels/pmodel{mode}.dict"
        if not os.path.exists ( filename ):
            fname = f"{os.getcwd()}/{filename}"
            self.highlight ( "red",
                    f"cheat mode started with {mode}, but no {fname} found" )
            return
            # sys.exit(-1)
        # scom = ""
        with open ( filename, "rt" ) as f:
            m = eval ( f.read() )
            if type(m)==list:
                self.pprint ( f"With {filename} a path to list of models was supplied. Will pick the first." )
                m = m[0]
        self.initFromDict ( m, filename )

    def checkForNans ( self ):
        """ check protomodel for NaNs, for debugging only """
        for pid,m in self.M.masses.items():
            if np.isnan ( m ):
                self.pprint ( f"Checking for nans: mass of {int(pid)} is nan" )

    def get ( self ):
        """ since the shallowcopy business does not work as expected,
        here is a trivial way to overwrite the original protomodel.
        use as: protomodel = manipulator.get()
        """
        return self.M

    def setWalkerId ( self, Id : int ):
        """ set the walker id of protomodel """
        self.M.walkerid = Id

    def printCombo ( self, combo : Union[None,List[TheoryPrediction]] = None,
            detailed : bool = False ):
        """ pretty print prediction combos.
            If None, print best combo
        :param combo: None, to print the best combo, else print that combo
        :param detailed: if true, print more detailed report
        """
        print ( "best combo:" )
        if combo == None:
            combo = self.M.bestCombo
        for i in combo:
            txns = ",".join ( set ( map ( str, i.txnames ) ) )
            dId = i.dataId() if i.dataId() != None else "UL"
            print ( f" `- {GREEN}{i.analysisId()}:{dId}: {txns}{RESET}" )
            line = "        "
            if detailed:
                import math
                dtype = i.dataType()
                if dtype == "efficiencyMap":
                    dI = i.dataset.dataInfo
                    eBG = dI.expectedBG
                    if eBG == int(eBG):
                        eBG=int(eBG)
                    bgErr = dI.bgError
                    if bgErr == int(bgErr):
                        bgErr=int(bgErr)
                    toterr = math.sqrt ( bgErr**2 + eBG )
                    line += f"obs={dI.observedN} exp={eBG:.2f}+-{bgErr}"
                    if toterr > 0.:
                        line += f" Z={RED}{(dI.observedN - eBG ) / toterr :.1f}*sigma{RESET}"
                    print ( line )
                if dtype in [ "upperLimit", "combined" ]:
                    try:
                        eUL = i.getUpperLimit ( expected = True ).asNumber(fb)
                    except Exception as e:
                        eUL = i.getUpperLimit ( evaluationType = True ).asNumber(fb)
                    oUL = i.getUpperLimit ( ).asNumber(fb)
                    sigma_exp = eUL / 1.96 # the expected scale, sigma
                    Z = ( oUL - eUL ) / sigma_exp
                    line += f"obs={oUL:.1f}*fb exp={eUL:.1f}*fb Z={RED}{Z:.1f}*sigma{RESET}"
                    print ( line )

            allpids = list( getAllPidsOfTheoryPred ( i ) )
            pidline = f"        pids:"
            for pid in allpids[:2]:
                pidline += f" {self.namer.asciiName(pid)}"
            if len(pidline) > 80:
                pidline=f"{pidline[:76]} ..."
            if len(allpids)>3:
                          pidline += f" ..."
            print ( pidline )

    def printAllTheoryPredictions ( self, detailed : bool = False ):
        """ pretty print all theory predictions for the model
        :param detailed: if true, give more details
        """
        print ( "theory predictions:" )
        combo = self.M.ul_critic_tpList
        for c in combo:
            i = c["tp"]
            dId = i.dataId() if i.dataId() != None else "UL"
            txns = ",".join ( set ( map ( str, i.txnames ) ) )
            print ( f" - {i.analysisId()}:{dId}:{txns}" )
            if detailed:
                robs, rexp = "n/a", "n/a"
                if c['robs'] is not None:
                    robs = f"{c['robs']:.1f}"
                if c['rexp'] is not None:
                    rexp = f"{c['rexp']:.1f}"
                print ( f"   - robs={robs} rexp={rexp}" )
            allpids = list ( getAllPidsOfTheoryPred ( i ) )
            pidline  = f"        pids:"
            for pid in allpids[:2]:
                pidline += f" {self.namer.asciiName(pid)}"
                if len(pidline) > 80:
                    pidline=f"{pidline[:76]} ..."
            if len(allpids)>3:
                pidline += ( " ..." )
            print ( pidline )

    def removeIllegalBRs(self, rescaleSSMs=False, protomodel = None):
        """ remove all illegal decays and decays of frozen particles. Renormalize all branchings """

        if protomodel is None:
            protomodel = self.M

        for pid in protomodel.frozenParticles():
            if pid in protomodel.decays:
                protomodel.decays.pop(pid)

        olddecays = list ( protomodel.decays.keys() )
        #Loop over all decays:
        for pid in olddecays:
            if not pid in protomodel.decays:
                ## self.normalizeBranchings is allowed to take out
                ## pids, thus we need to check
                continue
            #Get allowed decay channels:
            openChannels = protomodel.getOpenChannels(pid)
            #Check if any of the existing decays are forbidden:
            delDecays = [dpid for dpid in protomodel.decays[pid] if not  dpid in openChannels]
            #if delDecays != []: print(f"removing decays of {pid}: {delDecays}")
            for dpid in delDecays:
                protomodel.decays[pid].pop(dpid)
            #if delDecays != []: print(f"remaining decays of {pid}: {protomodel.decays[pid]}")
            #Make sure to normalize the branchings
            self.normalizeBranchings(pid, rescaleSSMs=rescaleSSMs, protomodel=protomodel)

    def initBranchings ( self, pid : int ,
            protomodel : Union[ProtoModel,None] = None,
            randomly : bool = True ) -> bool:
        """ Initialize BRs to different open decay channels for pid.
        Either assign 'democratic' BRs or assign random BRs.
        :param pid: particle id
        :param protomodel: if None, then use self.M
        :param randomly: if true, initialize to random values, else democratically

        :returns: true if really initialized, false if could not initialize
        """
        if protomodel is None:
            protomodel = self.M

        #Do not modify the LSP decays
        if pid in self.decaylessParticles:
            return False
        pid_name = self.namer.asciiName(pid)

        #Erase BRs (if any has been stored)
        protomodel.decays[pid] = {}

        #Get the allowed decay channels:
        openChannels = protomodel.getOpenChannels(pid)

        dkeys = set()
        for dpid in openChannels:
            dk = self.M.decay_keys[pid][dpid]
            dkeys.add(dk)

        inv_decay_keys = self.M.inv_decay_keys[pid]

        br_tot = 0.
        # iterate through all decay keys
        for dkey in dkeys:
            br = 1.
            if randomly:
                # if randomly, then draw from exponentials,
                # so that we have a dominating decay
                br = random.expovariate()
            if len(inv_decay_keys[dkey])>1:
                br = br / len(inv_decay_keys[dkey])
            for dpids in inv_decay_keys[dkey]:
                protomodel.decays[pid][dpids]=br
                br_tot += br

        if br_tot == 0.:
            self.error ( f"could not initialize {pid_name}: br_tot={br_tot}." )
            self.error ( f"openChannels were {openChannels}. will freeze." )
            self.freezeParticle ( pid )
            return False
        ## finally normalize
        for dpids,br in protomodel.decays[pid].items():
            protomodel.decays[pid][dpids]=br/br_tot
        return True

    def normalizeBranchings(self, pid : int, rescaleSSMs : bool =False,
            protomodel : Union[ProtoModel,None] = None) -> bool:
        """ normalize branchings of a particle if the total BR is differs
        from 1.0.

        :param pid: Particle to have their branchings normalized.
        If pid = None, normalize all decays.
        :param rescaleSSMs: if True, rescale the corresponding signal strength
        multipliers, so that sigma x br stays the same.

        :returns: true if succesful
        """

        if not protomodel:
            protomodel = self.M

        if not pid in protomodel.decays:
            protomodel.pprint(f"When attempting to normalize: {pid} not in decays")
            return False

        BRtot = sum(protomodel.decays[pid].values())
        if BRtot == 0:
            self.log ( f"the decayless particles are {self.namer.asciiName(self.M.decaylessParticles)} [{self.M.decaylessParticles}]" )
            if pid not in self.M.decaylessParticles:
                #print(f"decay of {pid}: {protomodel.decays[pid]}")
                self.log(f"decay of {pid}: {protomodel.decays[pid]}")
                protomodel.pprint ( f"When attempting to normalize: total BR of ({pid}) is zero, and it is not in decaylessParticles. we need to take out {pid}." )
                ## we need to freeze also <pid> now
                ## (since we have no sensible channels anymore)
                self.freezeParticles ( pid, force=True, protomodel=protomodel )
            return False

        if abs(BRtot-1.0) < 1e-4:
            #BRs are already normalized.
            return True

        self.log ( f"normalized branchings of {self.namer.asciiName(pid)} by {BRtot:.2f}" )

        for dpid in protomodel.decays[pid]:
            protomodel.decays[pid][dpid] *= 1/BRtot

        ## adjust the signal strength multipliers to keep everything else
        ## as it was
        if not rescaleSSMs:
            return True

        #rescaling ssms?
        for pidpair,ssm in protomodel.ssmultipliers.items():
            if pidpair in [ (pid,pid),(-pid,-pid),(-pid,pid),(pid,-pid) ]:
                newssm = min(1e5,ssm*BRtot*BRtot) #Rescale pair production by BRtot^2
            elif (pid in pidpair) or (-pid in pidpair):
                newssm = min(1e5,ssm*BRtot) #Rescale associated production by BRtot
            else:
                continue
            protomodel.ssmultipliers[pidpair]=newssm

        return True

    def initSSMFor(self, pid : int, protomodel : Union[ProtoModel,None] =None,
            ssmSigma : float = 1.0, cap_ssm : float = 100. ):
        """ Initialize SSM multipliers (for pair production of
        particle/anti-particle): new ssm = lognorm.rvs(1.0, ssmSigma)
        """
        p_name = self.namer.asciiName(pid)
        self.log ( f"initSSMFor {p_name}({pid}) with ssmSigma={ssmSigma:.1} cap_ssm={cap_ssm:.1f}" )

        if protomodel is None:
            protomodel = self.M

        unfrozen = protomodel.unFrozenParticles()
        pidpair = protomodel.getAllowedProdModes()

        #num_ssms = np.random.geometric(p=0.5)
        #p = len(self.M.ssmultipliers.keys())/20 ? -> make it difficult to add new prod modes if already lot present?
        #init_pidpair = np.random.choice(pidpair, size=num_ssms, replace=False)

        for ppair in pidpair:
            if abs(ppair[0]) in unfrozen and abs(ppair[1]) in unfrozen:
                if abs(ppair[0]) == pid or abs(ppair[1]) == pid:
                    ssm = float(lognorm.rvs(s = ssmSigma, scale = 1.0))   #center ssm around 1.0, better to have log scale
                    if ssm > cap_ssm: ssm = cap_ssm
                    protomodel.ssmultipliers[ppair] = ssm
                    self.log ( f"setting ssm of {ppair}({self.namer.asciiName(ppair)}) to {ssm:.2f}" )


    def describe ( self, allTheoryPredictions : bool = False ):
        """ lengthy description of protomodel
        :param allTheoryPredictions: if true, list all theory preds
        """
        sK, sTL = formatObject(self.M.K), formatObject(self.M.TL)
        print( f'\nK = {sK}, TL = {sTL}, muhat = {self.M.muhat:1.2f}, mumax={self.M.mumax:1.3g}' )
        print('  * Best Combo:')
        for tp in self.M.bestCombo:
            txns = ",".join ( set ( map ( str, tp.txnames ) ) )
            eUL = "no ULexp"
            anaId = tp.expResult.globalInfo.id
            dt = tp.dataType( short=True )
            fullId = f"{anaId}:{dt}"
            if hasattr ( tp, "expectedUL" ) and type(tp.expectedUL) != type(None):
                eUL = f"UL_exp={tp.expectedUL.asNumber(fb):1.2f}"
            if dt in [ "comb", "combined" ]:
                eUL = f"{tp.getUpperLimit ( expected=True ).asNumber(fb):1.2g}*fb"
            if dt in [ "em", "efficiencyMap" ]:
                dI = tp.dataset.dataInfo
                pred = f"{float ( (tp.xsection.value * tp.dataset.globalInfo.lumi).asNumber() ):1.2g}"
                fullId = f"{anaId}:{tp.dataset.dataInfo.dataId}"
                print(f'     - {fullId} [{txns}] obsN={dI.observedN} expBG={dI.expectedBG}+/-{dI.bgError} pred={pred}' )
            else:
                pred=f"{tp.xsection.value.asNumber(fb):.2g}*fb"
                UL = f"{tp.getUpperLimit().asNumber(fb):.2g}*fb"
                print(f'     - {fullId} [{txns}] pred={pred} UL={UL} eUL={eUL}' )

        print ( )
        print('  * Constraints:')
        for tp in sorted( self.M.ul_critic_tpList, key = lambda x: x['robs'], reverse=True ):
            if not allTheoryPredictions and tp['robs'] < 1.0:
                # if not all theory predictions are asked for, only do r>=1
                continue
            txns = ",".join ( set ( map ( str, tp[2].txnames ) ) )
            eUL = ""
            # eUL = ", no ULexp"
            anaId = tp[2].expResult.globalInfo.id
            anaId+=f":{tp[2].dataType()}"
            dt = tp[2].dataType( short=True )
            if hasattr ( tp[2], "expectedUL" ) and type(tp[2].expectedUL) != type(None):
                eUL = f", UL_exp={tp[2].expectedUL.asNumber(fb):1.2g}*fb"
            if dt in [ "comb", "combined" ]:
                eUL = f", UL_exp={tp[2].getUpperLimit ( expected=True ).asNumber(fb):1.2g}*fb"
            UL=f"{tp[2].upperLimit.asNumber(fb):1.2g}*fb"
            print( f'     - r={tp[2].getRValue():1.2f} {anaId} [{txns}] pred={tp[2].xsection.value.asNumber(fb):1.2g}*fb, UL={UL}{eUL}' )

    def rescaleSignalBy ( self, s : Union[float,None], cap_ssm : float = 100. ):
        """ multiply the signal strength multipliers with s
        :param s: rescale all signal strength multipliers with s.
        if None, rescale with self.M.muhat
        :param cap_ssm: allow only ssms < cap_ssm
        """
        if s == None:
            s = self.M.muhat
        if s == 0.:
            self.log ( "Rescaling by zero? Ignore." )
            return
        if abs ( s - 1.0 ) < 1e-5:
            return
        if s > cap_ssm: s = cap_ssm                                 #Do not rescale ssms by more than 100 (?)
        self.log ( f"rescaling signal by muhat of {s:.2f}" )
        self.M.rvalues = [r*s for r in self.M.rvalues[:]]
        self.M.muhat *= 1./s
        if self.M.mumax: self.M.mumax*= 1./s
        excl_prod_modes = self.getSSMsNotInBestCombo()
        self.M.rescaleXSecsBy(s, excl=excl_prod_modes, cap_ssm=cap_ssm)

        if hasattr(self.M,'ul_critic_tpList') and self.M.ul_critic_tpList is not None:
            for i,tp in enumerate(self.M.ul_critic_tpList[:]):
                rnew = tp['robs']*s
                if tp['rexp']:
                    rexpnew = tp['rexp']*s
                else:
                    rexpnew = tp['rexp']
                tpNew = tp['tp']
                tpNew.xsection *= s #rescale theoryPrediction
                #Remove likelihood and chi2, since they are no longer valid
                #if hasattr(tpNew,'likelihood'):
                #    del tpNew.likelihood
                #if hasattr(tpNew,'chi2'):
                #    del tpNew.chi2
                self.M.ul_critic_tpList[i] = { "robs": rnew,"rexp": rexpnew, "tp": tpNew }
        if hasattr(self.M,'bestCombo') and self.M.bestCombo is not None:
            for tp in self.M.bestCombo:
                tp.xsection *= s #rescale theoryPrediction
                #Remove likelihood and chi2, since they are no longer valid
                #if hasattr(tp,'likelihood'):
                #    del tp.likelihood
                #if hasattr(tp,'chi2'):
                #    del tp.chi2

    def z_model(self, model_current, model_propose, force_move=False):

        a,b,c = 2,4,8
        shift_parameter = 8
        #print(f"Move {move}")
        #total number of particles
        n_par_current = len(model_current.unFrozenParticles() )
        n_par_propose = len(model_propose.unFrozenParticles())

        #total number of non-trivial br
        n_decays_current = sum([len(dc.keys())- 1 for par, dc in model_current.decays.items() if par != LSP ])
        n_decays_propose = sum([len(dc.keys())- 1 for par, dc in model_propose.decays.items() if par != LSP ])

        #total number of production modes
        n_ssms_current = len([ssm for ssm in model_current.ssmultipliers.values() if ssm > 1e-04])
        n_ssms_propose = len([ssm for ssm in model_propose.ssmultipliers.values() if ssm > 1e-04])

        deg_prop = n_par_propose + n_decays_propose + n_ssms_propose
        deg_current = n_par_current + n_decays_current + n_ssms_current
        if self.run_mcmc:
            if deg_prop != deg_current:
                self.highlight("error", f"Error! Dimension Changed during mcmc walk!")
                if n_par_current != n_par_propose:
                    self.log(f"Particle content changed during mcmc walk. Prev num of par {n_par_current}, current num of par {n_par_propose}")
                if n_decays_current != n_decays_propose:
                    self.log(f"Decays changed during mcmc walk. Prev num of decays {n_decays_current}, current num of decays {n_decays_propose}")
                    print(f"current: {model_propose.decays}, {n_decays_propose}")
                    print(f"prev: {model_current.decays}, {n_decays_current}")
                if n_ssms_current != n_ssms_propose:
                    self.log(f"Production modes changed during mcmc walk. Prev num of prod modes {n_ssms_current}, current num of prod modes {n_ssms_propose}")
                return False
            return True
        #print(f"Current deg: {deg_current}")
        #print(f"Proposed deg: {deg_prop}")

        z_current = (n_par_current - shift_parameter)/a + n_decays_current/b + n_ssms_current/c
        z_propose = (n_par_propose - shift_parameter)/a + n_decays_propose/b + n_ssms_propose/c
        prob_12, prob_21 = 1.0, 1.0

        if z_propose > z_current:
            #print(f"z_propose {z_propose} > z_current {z_current}")
            prob_12 = (1 + np.exp(z_current))/(1 + np.exp(z_propose))            #Adding new degrees of freedom, i->i+1
            if force_move: prob_12 = 1.0
            prob_21 = (1 + np.exp(-z_propose))/(1 + np.exp(-z_current))          #i+1 -> i
            #print(f"prob12, prob21 {prob_12}, {prob_21}")
        elif z_propose < z_current:
            #print(f"z_propose {z_propose} < z_current {z_current}")
            prob_12 = (1 + np.exp(-z_current))/(1 + np.exp(-z_propose))          #Removing degrees of freedom, i->i+1
            if force_move: prob_12 = 1.0
            prob_21 = (1 + np.exp(z_propose))/(1 + np.exp(z_current))            #i+1 ->i
            #print(f"prob12, prob21 {prob_12}, {prob_21}")
        else:
            #print(f"z_propose {z_propose} = z_current {z_current}")     #No change in degrees of freedom
            prob_12, prob_21 = 1.0, 1.0

        if prob_12 > 1.0 or prob_21 > 1.0:
            self.highlight("error", f"Proposal probabilities greater than 1.0 while changing the dimensions: {prob_12},{prob_21}.")

        return prob_12, prob_21

    def proposal_density(self, move, force_move=False):         #shift prior from combiner here or vice versa
        """ Define the proposal density for changing the model
        :param move: specify which move we are making
        """
        if force_move:
            self.M = self.propose_model
            self.propose_model = self.M.copy()
            return True
        prob_12, prob_21 = self.z_model(self.M, self.propose_model)

        ## cap TL here or else you freeze in particle content for injected signals
        TL = self.M.TL
        if TL > 49:
            self.highlight("warning", f"TL {formatObject(TL)}>49. we cap at 49 for computing the probability of the proposal!" )
            TL = 49
        ### FIXME add warning if this upper limit gets triggered!

        if TL > 0.0: prob = min(1.0, prob_12/np.sqrt(TL))
        else: prob = min(1.0, prob_12)

        u = np.random.uniform(0,1)
        if u > prob:
            #print(f"Reject move {move}")
            self.log(f"u={u:.2f} > proposal prob {prob:.2f}. Reject Move.")
            self.propose_model = self.M.copy()
            self.proposal_ratio[move] = {'q': 1.0}
            #print(f"q_total after {move}: {self.proposal_ratio['q_total']}")
            return False
        else:
            #print(f"accept move {move}")
            self.log(f"u={u:.2f} <= proposal prob {prob:.2f}. Accept Move.")
            self.M = self.propose_model
            self.propose_model = self.M.copy()
            self.proposal_ratio[move]['q'] *= prob_21/prob_12
            if np.isnan(np.log(self.proposal_ratio[move]['q'])):
                self.log(f"Move {move} has nan value for {self.proposal_ratio[move]['q']}")
                self.proposal_ratio[move]['q'] = 1.0
            for parameter, q_par in self.proposal_ratio[move].items():
                #print(f"Parameter : {parameter}, q_par {q_par}")
                self.proposal_ratio['q_total'] *= q_par
            #print(f"log q_total after {move}: {np.log(self.proposal_ratio['q_total'])}")
            if np.isnan(np.log(self.proposal_ratio['q_total'])):
                self.log(f"Move {move} has nan value for qtotal {self.proposal_ratio['q_total']}")
                self.proposal_ratio['q_total'] = 1.0
            self.log(f"Log of total proposal ratio after {move}: {np.log(self.proposal_ratio['q_total']):.2f}")
            return True

    def randomlyChangeModel(self,force_move : bool = False, probBR : float = 0.2,
            probSS : float = 0.25, probSSingle : float = 0.8, ssmSigma : float = 1.0,
            probMerge : float = 0.05, sigmaFreeze : float = 0.5,
            probMassive : float = 0.3, probMass : float = 0.05, run_mcmc= False, cap_ssm=100):
        """Randomly modify the proto-model following the steps:

        1) A random particle can be unfrozen with a probability
        controlled by sigmaUnFreeze
        2) A random BR can be modified (with probability probBR)
        3) A random signal strenght can be modified
        the probability is controlled by probSS, probSSingle and ssmSigma
        4) Particles can be merged (with probability probMerge)
        5) A random particle can be frozen
        the probability is controlled by sigmaFreeze and probMassive
        6) A random mass can be changed by a maximum value of dx
        with probability of probMass
        """

        self.proposal_ratio = {'add_par':{'q':1.0}, 'rem_par':{'q':1.0}, 'br':{'q':1.0}, 'ssm':{'q':1.0}, 'q_total':1.0}
        self.run_mcmc = run_mcmc
        if self.run_mcmc: old_model = self.M.copy()
        else: self.propose_model = self.M.copy()
        changeDesc = {}
        nChanges = 0

        # If TL < = 0, try to drastically change model, else will be stuck in a model with low TL for many steps
        if self.M.TL == None or self.M.TL <=0 :
            probSS = 1.0
            probBR = 1.0
            probMass = 1.0
            force_move = True
            #do we want to not freeze particles? -> not freezing if less than or equal to 3 particles

        if not self.run_mcmc:
            recentlyUnfrozen = self.randomlyUnfreezeParticle(cap_ssm=cap_ssm)
            if recentlyUnfrozen:
                accept_move = self.proposal_density(move='add_par', force_move=force_move)
                if accept_move:
                    nChanges += 1
                    self.log(f"Accept unfreezing of {self.namer.asciiName(recentlyUnfrozen)} ({recentlyUnfrozen})")
                else:
                    self.log(f"Reject unfreezing {self.namer.asciiName(recentlyUnfrozen)} ({recentlyUnfrozen})")
                    recentlyUnfrozen = None

            frozenParticles = self.randomlyFreezeParticle(recentlyUnfrozen=recentlyUnfrozen)
            if frozenParticles:
                accept_move = self.proposal_density(move='rem_par', force_move=force_move)
                if accept_move:
                    nChanges += 1
                    self.log(f"Accept freezing of {', '.join(map(str,frozenParticles))} ({self.namer.asciiName(frozenParticles)})")
                    #print(f"Protomodel now: {self.M.unFrozenParticles()}")
                else:
                    self.log(f"Reject freezing of {', '.join(map(str,frozenParticles))} ({self.namer.asciiName(frozenParticles)})")

            changes = self.randomlyChangeBranchings(protomodel=self.propose_model, prob=probBR)
            if changes > 0:
                accept_move = self.proposal_density(move='br', force_move=force_move)
                if accept_move :
                    nChanges += 1          #;print(f"Protomodel after: {self.M.decays}")
                    self.log("Accept changes in branchings")
                else: self.log(f"Reject changing branchings")

            changes = self.randomlyChangeSignalStrengths(protomodel=self.propose_model, prob = probSS, probSingle = probSSingle, ssmSigma = ssmSigma, cap_ssm=cap_ssm)
            if changes > 0:
                accept_move = self.proposal_density(move='ssm', force_move=force_move)
                if accept_move:
                    nChanges += 1       #; print(f"Protomodel after: {self.M.ssmultipliers}")
                    self.log("Accept changes in ssm")
                else: self.log("Reject changes in ssm")


        else:
            changes = self.randomlyChangeBranchings(protomodel = self.M, prob=0.4, zeroBRprob = 0., singleBRprob = 0., addBRprob = 0.)
            nChanges += changes
            changes = self.randomlyChangeSignalStrengths(protomodel=self.M, prob=0.4, probSingle=1.0, ssmSigma=ssmSigma)
            nChanges += changes


        if not nChanges: #If nothing has changed, force a random change of masses
            changes = self.randomlyChangeMasses(prob=1.0)
            nChanges += changes
        else: #Change masses with 5% probability
            if self.run_mcmc: self.randomlyChangeMasses(prob = 0.7)
            else: self.randomlyChangeMasses(prob = probMass)

        if self.run_mcmc:
            check_dim_change = self.z_model(old_model, self.M)
            if not check_dim_change:
                self.log("Returning to previous model")
                self.M = old_model
        #print(f"q_total = {self.proposal_ratio['q_total']}")
        #Update cross-sections (if needed)
        self.M.getXsecs()

    def randomlyUnfreezeParticle ( self, cap_ssm : float = 100. ) -> int:
        """ Unfreezes a (random) frozen particle according to gaussian distribution
            with a width of <sigma>.

        :param cap_ssm: allow a maximum change in signal strength multipliers
        by a factor of cap_ssm

        :returns: 1 if a particle got unfrozen, 0 if not.
        """

        #Decide whether to unfreeze according to the number of active particles
        #(always unfreeze if the model only has one particle) -> still keep?
        nUnfrozen = len( self.M.unFrozenParticles() )

        # Randomly select the pid:
        frozen = self.M.frozenParticles()
        if self.M.environ.forbiddenparticles != []:
            frozen = [par for par in self.M.frozenParticles() if par not in self.M.environ.forbiddenparticles]

        if len(frozen)==0:
            return None
        pid = int(np.random.choice ( frozen ))

        if pid in self.M.environ.forbiddenparticles:
            self.log ( f"wanted to unfreeze {self.namer.asciiName(pid)} but its forbidden" )
            return None

        #Check for canonical ordering.
        #If pid matches the heavier state and the lighter state is frozen,
        #unfreeze the lighter state instead
        for pids in self.canonicalOrder:
            if pid == pids[1] and (pids[0] in frozen):
                pid = pids[0] #Unfreeze the lighter state
                break

        self.log ( f"Propose unfreezing {self.namer.asciiName(pid)}({pid})" )
        #print(f"Propose unfreezing {self.namer.asciiName(pid)}" )
        unfrozen = self.unFreezeParticles( pid, protomodel = self.propose_model, cap_ssm=cap_ssm)
        return unfrozen

    def randomlyChangeBranchings ( self, protomodel = None, prob : float =0.2,
            zeroBRprob : float = 0.05, singleBRprob : float = 0.05,
            addBRprob : float = 0.1 ) -> int:
        """ randomly change the branchings of a single particle

        :param prob: Probability for changing a branching ratio
        :param zeroBRprob: With zeroBRprob probability, close decay channel
        :param singleBRprob: With probability singleBRprob,
        keep only one decay channel
        :param addBRprob: With probability addBRprob, add a new decay channel
        for the pid
        :returns: number of changes
        """
        if protomodel is None:
            protomodel = self.M
        uBranch = np.random.uniform(0,1)
        if uBranch < (1-prob):
            self.log("Not changing branchings")
            return 0

        unfrozenparticles = self.M.unFrozenParticles( withLSP=False )
        if len(unfrozenparticles)<2:
            self.log( "Not enough unfrozen particles to change random branching" )
            return 0
        p = int(np.random.choice ( unfrozenparticles ))
        if not p in self.M.decays.keys():
            # this is now allowed to happen
            return 0
            self.highlight ( "error", f"why is {int(p)} not in decays?? {self.M.decays.keys()}" )
            # we dont know about this decay? we initialize with the default!

        return self.randomlyChangeBranchingOfPid ( p, protomodel, zeroBRprob, singleBRprob, addBRprob)

    def record ( self, change : str ):
        """ log the changes that have been performed on the model
        :param change: a string that describes my latest action
        """
        if not self.do_record:
            return
        self.recording.append ( change )
        if len(self.recording)>20: ## make sure this never explodes
            self.recording = self.recording[-20:]


    def randomlyChangeBranchingOfPid ( self, pid : int,
            protomodel : Union[ProtoModel,None] = None,
            zeroBRprob : float = 0.05, singleBRprob : float = 0.05,
            addBRprob : float = 0.1 ) -> int:
        """ randomly change the branching a particle pid
        :returns: number of changes
        """

        if protomodel is None:
            protomodel = self.M

        openChannels = protomodel.getOpenChannels(pid)
        dkeys = set()
        for dpid in openChannels:
            dk = protomodel.decay_keys[pid][dpid]
            dkeys.add(dk)

        dkeys = list(dkeys)
        self.log( f"Trying to change branchings of {pid} ({self.namer.asciiName(pid)})." )
        if len(openChannels) < 2:
            self.log( f"Number of open channels of {pid} is {len(openChannels)}. Cannot change branchings." )
            # not enough channels open to tamper with branchings!
            return 0

        self.proposal_ratio['br']['rem'] = 1.0
        self.proposal_ratio['br']['add'] = 1.0

        dx = 0.1/np.sqrt(len(openChannels)) ## maximum change per channel??

        #Keep only one channel (with probability singleBRprob)
        uSingle = np.random.uniform( 0., 1. )
        if uSingle < singleBRprob:
            p_name = self.namer.asciiName(pid)
            self.log(f"Keeping only one decay channel for {pid} ({p_name}).")
            #Choose random decay key:
            dk = np.random.choice(dkeys)
            #get decay channel assocaiated with key, make sure all channels assocaited with same key get same branchings
            decay_chan = [key for key,value in protomodel.decay_keys[pid].items() if value == dk]
            #Get proposal ratio for removing old br
            #proposal ratio for rem br =  p(i+1 -> i)/ p(i->i+1) = p(add br to i+1 to go to i)/p(rem br to go to i+1)
            #p(add) = p(addBR)
            #p(rem) = p(singleBR)p(not choosing dk) = singleBRprob * (1 - 1/(len(dkeys))
            prob_add = addBRprob
            prob_rem = singleBRprob * (1.0 - 1/len(dkeys))
            self.proposal_ratio['br']['rem'] *= prob_add/prob_rem

            protomodel.decays[pid] = {}
            br = 1.0/len(decay_chan)
            for dpid in decay_chan:
                self.record ( f"change decay of {self.namer.texName(pid,addDollars=True)} -> {self.namer.texName(dpid,addDollars=True)} to {br:.2f}" )
                p_name = self.namer.asciiName(pid)
                dp_name = self.namer.asciiName(dpid)
                self.log ( f"changed decay of {p_name} -> {dp_name} to {br:.2f}" )
                protomodel.decays[pid].update({dpid: br})

            #print(f"Prob to rem {self.proposal_ratio['br']['rem']}")
            return 1

        #Otherwise randomly change each channel(s) (based on the current BR)
        for dk in dkeys:
            oldbr = 0.

            #Check if decay channel already existed:
            decay_chan = [key for key,value in protomodel.decay_keys[pid].items() if value == dk]
            if decay_chan[0] in protomodel.decays[pid]:
                oldbr = self.M.decays[pid][decay_chan[0]]

            if oldbr > 0:
                #Close channel(s) (with zeroBRprob probability)
                uZero = np.random.uniform( 0., 1. )
                if uZero < zeroBRprob:
                    #Get proposal ratio for removing old br
                    #proposal ratio for rem br = p(add)/p(rem) = p(i+1 -> i)/ p(i->i+1)
                    #p(add) = p(addBRprob)
                    #p(rem) = p(zeroBR)
                    self.proposal_ratio['br']['rem'] *= addBRprob/zeroBRprob
                    for dpid in decay_chan:
                        self.record ( f"Removed decay {self.namer.texName(pid,addDollars=True)} -> {self.namer.texName(dpid,addDollars=True)} with br {oldbr:.2f}." )
                        self.log ( f"Removed decay {self.namer.asciiName(pid)} -> {self.namer.asciiName(dpid)} with br {oldbr:.2f}." )
                        if dpid in protomodel.decays[pid]:
                            protomodel.decays[pid].pop(dpid)
                        else:
                            self.debug ( f"{dpid} not in protomodel.decays[{pid}]{protomodel.decays[pid]}, decay_keys are {protomodel.decay_keys[pid]}" )
                    continue

                #Randomly change BR around old value
                #Min,Max = max(0.,oldbr-dx), min(oldbr+dx,1.)
                #br = float(np.random.uniform( Min, Max )/len(decay_chan))
                br = float(norm.rvs ( 1. / len(openChannels), np.sqrt ( .5 / len(openChannels) )  ))
                br = max(0.001, br)
                for dpid in decay_chan:
                    protomodel.decays[pid][dpid] = br
                    self.record ( f"Change branchings of {self.namer.texName(pid,addDollars=True)} -> {self.namer.texName(dpid,addDollars=True)} to {br:.2f}" )
                    self.log ( f"Changed  branchings of {self.namer.asciiName(pid)} -> {self.namer.asciiName(dpid)} to {br:.2f}" )

            else:
                #Add channel(s) (with addBRprob probability)
                uAdd = np.random.uniform( 0., 1. )
                if uAdd < addBRprob:
                    br = float(norm.rvs ( 1. / len(openChannels), np.sqrt ( .5 / len(openChannels) )  ))
                    br = max(0.001, br)
                    #Get proposal ratio for adding new br
                    #proposal ratio for add br = p(rem)/p(add) = p(i+1 -> i)/ p(i->i+1)
                    #p(add) = p(addBR)
                    #p(rem) = p(singleBR*p(rem br) + zeroBR)
                    prob_add = addBRprob
                    prob_rem = singleBRprob * (1.0 - 1/len(dkeys)) + zeroBRprob
                    self.proposal_ratio['br']['add'] *= prob_rem/prob_add
                    for dpid in decay_chan:
                        protomodel.decays[pid][dpid] = br
                        self.log ( f"Added decay of {self.namer.texName(pid,addDollars=True)} -> {self.namer.texName(dpid,addDollars=True)} with br {br:.2f}" )


        #Make sure there is at least one open channel:
        BRtot = sum(protomodel.decays[pid].values())
        if BRtot == 0.0:
            self.log(f"BRtot = 0 for {pid} ({self.namer.asciiName(pid)}). Randomly Choosing one decay channel.")
            dk = np.random.choice(dkeys)
            decay_chan = [key for key,value in protomodel.decay_keys[pid].items() if value == dk]
            br = 1.0/len(decay_chan)
            protomodel.decays[pid] = {}
            for dpid in decay_chan:
                self.record ( f"change decay of {self.namer.texName(pid,addDollars=True)} -> {self.namer.texName(dpid,addDollars=True)} to {br:.2f}" )
                self.log ( f"Changed decay of {self.namer.asciiName(pid)} -> {self.namer.asciiName(dpid)} to {br:.2f}" )
                protomodel.decays[pid].update({dpid: br})

        #Make sure BRsprot add up to 1:
        self.normalizeBranchings(pid, protomodel=protomodel)
        #print(f"Prob to rem {self.proposal_ratio['br']['rem']}")
        #print(f"Prob to add {self.proposal_ratio['br']['add']}")
        return 1

    def randomlyChangeSignalStrengths ( self, protomodel=None, prob : float =0.25,
            probSingle : float =0.8, ssmSigma : float = 1.0,
            cap_ssm : float =100. ) -> int:
        """ randomly change one of the signal strengths according to a gaussian
        distribution centered around the original SSM.

        :param prob: Probability for changing the signal strengths
        :param probSingle: Probability for changing the signal strength of a
        single particle
        :param ssmSigma: Width for the gaussian, put to 1.0?

        :returns: 1 if something got changed, else 0
        """

        uSSM = np.random.uniform(0,1)
        if uSSM < (1-prob):
            self.log("Not changing ssm")
            return 0

        if protomodel is None:
            protmodel = self.M

        if np.random.uniform(0,1) < probSingle:
            return self.randomlyChangeSSOfOneParticle(protomodel = protomodel,
                    ssmSigma=ssmSigma, cap_ssm=cap_ssm)

        unfrozenparticles = self.M.unFrozenParticles( withLSP=self.M.allowN1N1Prod )

        if len(unfrozenparticles)<2:
            self.log ( "Not enough unfrozen particles to change random ssm" )
            return 0

        #first get allowed production modes
        prodModes = protomodel.getAllowedProdModes()
        #filter the production modes which have the unfrozen particles
        pidpair = set()
        for ppair in prodModes:
            if abs(ppair[0]) in unfrozenparticles and abs(ppair[1]) in unfrozenparticles:
                pidpair.add(ppair)

        pidpair = list(pidpair)

        #prob to add a ssm occurs 70% of the time. This is done if none of three random moves below are implemented
        #prob to rem a ssm occurs only 10% of the time, in the first random move implemented below (under a >0.9)
        #Do random moves
        a = np.random.uniform ( 0., 1. )
        if a > .9: ## sometimes, just knock out a random SSM
            prod_list = list(protomodel.ssmultipliers.keys())
            if len(prod_list) == 0: return 0
            random_ind = int(np.random.choice(len(prod_list)))
            randomProd = prod_list[random_ind]
            self.log(f"Remove prod mode {self.namer.texName(randomProd,addDollars=True)}" )
            protomodel.ssmultipliers.pop(randomProd)
            #get proposal ratio
            #proposal ratio for rem ssm = p(i+1 -> i)/ p(i->i+1) = p(add)/p(rem)
            #p(add) = 0.7 * 1/(len(pidpair))  (look code below)
            #p(rem) = 0.1 * p(rem randomProd) = 0.1 * 1/len(prod_list)
            prob_add = 0.7/len(pidpair)
            prob_rem = 0.1/len(prod_list)
            self.proposal_ratio['ssm']['rem'] = prob_add/prob_rem
            #print(f"Prob to rem {self.proposal_ratio['ssm']['rem']}")
            return 1
        if a < .1: ## sometimes, just try to set ssm to 1.
            prod_list = list(protomodel.ssmultipliers.keys())
            if len(prod_list) == 0: return 0
            random_ind = int(np.random.choice(len(prod_list)))
            randomProd = prod_list[random_ind]
            self.log(f"Change ssm of {self.namer.texName(randomProd,addDollars=True)} to 1." )
            protomodel.ssmultipliers[randomProd]=1.
            return 1
        if .1 < a < .2: ## sometimes, just try to set to ssm of different particle
            prod_list = list(protomodel.ssmultipliers.keys())
            if len(prod_list) == 0: return 0
            random_ind = int(np.random.choice(len(prod_list)))
            randomProd = prod_list[random_ind]
            v = np.random.choice ( list ( protomodel.ssmultipliers.values() ) )
            self.log ( f"Change ssm of {self.namer.texName(randomProd,addDollars=True)} to {v:.2f}" )
            protomodel.ssmultipliers[randomProd]= float(v)
            return 1

        #Else either add new ssm or change existing ssm of pair
        #Randomly choose which process pids to change:
        random_ind = int(np.random.choice(len(pidpair)))
        pair = pidpair[random_ind]

        newSSM = 1.0
        if not pair in protomodel.ssmultipliers:
            newSSM = float(lognorm.rvs(s = ssmSigma, scale = 1.0)) #center ssm around 1.0, better to have log scale
            if newSSM > cap_ssm: newSSM = cap_ssm
            protomodel.ssmultipliers[pair] = newSSM
            #get proposal ratio
            #proposal ratio for add ssm = p(i+1 -> i)/ p(i->i+1) = p(rem)/p(add)
            #p(add pair) = p(adding) P(choosing pair) = 0.7*(1/(len(pidpair))
            #p(rem) = p(closing)p(rem pair) = 0.1 * 1/len(ssms)
            prob_add = 0.7/len(pidpair)
            prob_rem = 0.1/len(protomodel.ssmultipliers.keys())
            self.proposal_ratio['ssm']['add'] = prob_rem/prob_add
            #print(f"Adding new pair of ssm {pidpair} with ratio {self.proposal_ratio['ssm']['add']}")
            self.log( f"Add new prod mode {self.namer.texName(pair,addDollars=True)} with ssm {newSSM}" )
        else:
            newSSM = float(lognorm.rvs(s = ssmSigma, scale = 1.0))
            if newSSM > cap_ssm: newSSM = cap_ssm
            protomodel.ssmultipliers[pair] = newSSM
            #self.changeSSM(pair,newSSM)
            self.log ( f"Changing ssm of {self.namer.asciiName(pair[0])},{self.namer.asciiName(pair[1])}: {newSSM:.2f}." )
            self.record ( f"change ssm of {self.namer.texName(pair[0])},{self.namer.texName(pair[1])} to {newSSM:.2f}." )
        return 1

    def randomlyChangeSSOfOneParticle ( self, pid = None, protomodel=None, ssmSigma=1.0, cap_ssm=100. ):
        """ randomly change the SS's consistently for one pid
        :param pid: change for this pid. If None, change of a random pid.
        """
        if protomodel is None:
            protomodel = self.M

        unfrozenparticles = self.M.unFrozenParticles( withLSP=self.M.allowN1N1Prod )

        if len(unfrozenparticles)<2:
            self.log ( "Not enough unfrozen particles to change random ssm" )   #why? we are changing only for 1 particle?
            return 0

        p = int(np.random.choice ( unfrozenparticles ))
        if pid != None: p = pid
        self.log (f"Changing all ssms of {self.namer.asciiName(p)} ({p})" )

        ssms = []
        for dpd,v in protomodel.ssmultipliers.items():
            if p in dpd or -p in dpd:
                newSSM = float(lognorm.rvs(s = ssmSigma, scale = 1.0))
                if newSSM > cap_ssm: newSSM = cap_ssm
                protomodel.ssmultipliers[dpd]= newSSM
                self.log (f"Changing ssm of {self.namer.asciiName(dpd)} ({dpd}) to newSSM" )
                #self.changeSSM ( dpd, newssm )
                ssms.append ( newSSM )

        return 1

    def pidPairIsInSSMs ( self, pids : Tuple ) -> bool:
        """ is a given pid pair in SSMs? """
        if pids[1] < pids[0]:
            pids = ( pids[1], pids[0] )
        return pids in self.M.ssmultipliers.keys()

    def changeSSM ( self, pids : Tuple, newssm, recursive : bool = True,
                    verbose : bool = True, cap_ssm=100. ):
        """ change the signal strength multiplier of pids to newssm,
            if we have stored xsecs, we correct them, also

        :param pids: Tuple of particle ids, e.g. (1000024,1000023)
        :param recursive: if true, then change also for all other signs, e.g.
        (-pids[0],-pids[1]), etc
        :param verbose: if False, then dont mention it. used for recursive.
        """
        if type(pids) != tuple:
            self.highlight ( "error", "when changing SSMs, need to supply PIDs as a tuple!" )
            return
        if len(pids)!= 2:
            self.highlight ( "error", "when changing SSMs, need to supply PIDs as a tuple of two pids!" )
            return
        if pids[1] < pids[0]:
            self.debug ( "warn", "when changing SSMs, pids are wrongly ordered. Reverting them." )
            pids = ( pids[1], pids[0] )

        if not pids in self.M.ssmultipliers:
            self.highlight ( "warn", f"when changing SSMs, cannot find {str(pids)}. not changing anything." )
            return
        oldssm = self.M.ssmultipliers[pids]
        if newssm > cap_ssm:
            newssm = cap_ssm
        if verbose:
            self.record ( f"change ssm of {self.namer.texName(pids,addDollars=True)} to {newssm:.2f}" )
        self.M.ssmultipliers[pids]=newssm
        if (oldssm + newssm) > 0.:
            if 2. * abs ( oldssm - newssm ) / ( oldssm + newssm ) > 1e-4:
                if verbose:
                    self.highlight ( "info", f"changing ssm of {self.namer.asciiName(pids)} from {oldssm:.2f} to {newssm:.2f}" )

        if not recursive:
            return
        ## change for all signs
        if self.pidPairIsInSSMs ( (-pids[0],pids[1]) ):
            self.changeSSM ( (-pids[0],pids[1]), newssm, recursive=False, verbose=False )
        if self.pidPairIsInSSMs ( ( pids[0],- pids[1]) ):
            self.changeSSM ( (pids[0],- pids[1]), newssm, recursive=False, verbose=False )
        if self.pidPairIsInSSMs ( ( - pids[0],- pids[1]) ):
            self.changeSSM ( (-pids[0],- pids[1]), newssm, recursive=False, verbose=False )


    def randomlyFreezeParticle ( self, recentlyUnfrozen = None ) -> list:
        """ freezes a random unfrozen particle according to gaussian distribution with width sigma.
        :param recentlyUnfrozen: do not freeze recently unfrozen particle if not None

        :returns: list of pids that got frozen
        """

        nUnfrozen = len( self.M.unFrozenParticles() )
        #Always keep at least 2 particles
        if nUnfrozen <= 2:
            self.log("Not freezing: Only 2 particles present.")
            return None

        if nUnfrozen <=3 and (self.M.TL == None or self.M.TL <=0 ):
            self.log(f"Not freezing: Only {nUnfrozen} particles and low TL {self.M.TL}.")
            return None

        unfrozen = [pid for pid in self.M.unFrozenParticles( withLSP = False ) if pid != recentlyUnfrozen]

        if len(unfrozen) == 0:
            self.log("Not freezing: Cannot freeze recently unfrozen particle {recentlyUnfrozen}.")
            return None

        pid = int(np.random.choice ( unfrozen ))
        allpids = self.forcedMassDegeneratePids ( pid )
        self.log ( f"randomly freezing particle {pid}->{allpids} {self.M.forced_degeneracies} {self.M.environ.templateName}" )
        frozen = self.freezeParticles ( pid, protomodel=self.propose_model )
        return frozen

    def freezeMostMassiveParticle ( self, protomodel=None) -> list:
        """ freezes the most massive unfrozen particle

        :returns: list of pids of particles that have been frozen
        """
        self.log ( f"freeze most massive particle" )

        if protomodel is None:
            protomodel = self.M

        unfrozen = protomodel.unFrozenParticles( withLSP=False )
        if len(unfrozen)<2:
            return 0                               #freeze only if at least 3 unfrozen particles exist

        pid,minmass=0,0
        for i in unfrozen:
            if self.M.masses[i]>minmass:           #check always with the current model self.M while trying to freeze massive particle?
                minmass = self.M.masses[i]
                pid = i

        protomodel.log ( f"Propose freezing most massive particle {pid}({self.namer.asciiName(pid)}) minmass=({minmass:.1f})" )
        frozen = self.freezeParticles ( pid, protomodel = protomodel)
        return frozen

    def freezeParticles ( self, pid : int, force : bool = False,
            protomodel = None, merge : bool = False,
            reassignPID : Union[None,Tuple[int]] = None ) -> list:
        """ freeze particle pid, and all its forced mass degenerate
        siblings

        :param pid: PID to be frozen
        :param force: If False, will only freeze the particle if it does not
        violate the canonical order (e.g. will not freeze stop1 if stop2 is
        unfrozen) and the model contains at least 3 particles.
        :param reassignPID: if tuple of two pids and not none,
        then reassign pid #1 to #2 in ssmultipliers
        :param merge: something about the proposal density?

        :returns: list of pids that really were frozen out
        """
        frozen = []
        allpids = self.forcedMassDegeneratePids(pid)
        for ipid in allpids:
            frozen += self.freezeParticle ( ipid, force=True,
                    protomodel=protomodel, reassignPID = reassignPID )
        return frozen


    def freezeParticle ( self, pid : int, force : bool = False,
            protomodel = None, merge : bool = False,
            reassignPID : Union[None,Tuple[int]] = None ) -> list:
        """ freeze particle pid, take care of offshell removal, and
            branching normalization

        :param pid: PID to be frozen
        :param force: If False, will only freeze the particle if it does not
        violate the canonical order (e.g. will not freeze stop1 if stop2 is
        unfrozen) and the model contains at least 3 particles.
        :param merge: something about the proposal density?
        :param reassignPID: if tuple of two pids and not none,
        then reassign pid #1 to #2 in ssmultipliers

        :returns: list of pids that really were frozen out
        """
        self.log ( f"freeze {pid}({self.namer.asciiName(pid)})" )

        if protomodel is None:
            protomodel = self.M

        #Check for canonical ordering.
        unfrozen = protomodel.unFrozenParticles( withLSP=False )
        if not force:
            if len(unfrozen) < 2:
                self.log("Not freezing: Only 2 particles present.")
                return []
            #If pid matches the lighter state and the heavier state is unfrozen,
            #do not freeze the particle
            for pids in self.canonicalOrder:
                if pid == pids[0] and pids[1] in unfrozen:
                    self.log(f"Not freezing: Tried to freeze {pids[0]} but {pids[1]} is unfrozen")
                    return []
        #protomodel.log ( f"Freezing {self.namer.asciiName(pid)}" )
        #self.record ( f"freeze {self.namer.texName(pid,addDollars=True)}" )
        #Remove pid from masses, decays and signal multipliers:
        if not force: self.log(f"Propose freezing pid: {pid}({self.namer.asciiName(pid)})")
        else: self.log(f"Freezing pid: {self.namer.asciiName(pid)}({pid})")
        #print(f"Propose freezing pid: {pid}")

        #get total num of frozen and unfrozen par for proposal ratio
        num_unfrozen = len(unfrozen)
        num_frozen = len(protomodel.frozenParticles())
        for pids in self.canonicalOrder:
            if pids[0] in self.M.environ.forbiddenparticles and pids[1] in self.M.environ.forbiddenparticles: continue
            if pids[0] in unfrozen and pids[1] in unfrozen:
                num_unfrozen -= 1       #num of par to freeze is smaller (i.e cannot freeze pids[0] while pids[1] is unfrozen)
            if pids[0] in protomodel.frozenParticles() and pids[1] in protomodel.frozenParticles():
                num_frozen -= 1         #num of par to unfreeze is smaller (i.e cannot unfreeze pids[1] while pids[0] is frozen)

        if self.M.environ.forbiddenparticles != [] :
            num_frozen -= len(self.M.environ.forbiddenparticles)
            if num_frozen < 0: self.log(f"NUm frozen {num_frozen} <0 !"); num_frozen = 1
        #proposal ratio = p(i+1 -> i)/p(i->i+1) = p(add pid from frozen)/p(rem pid from unfrozen) = (1/(n_fr+1))/(1/n_un)
        if merge: self.proposal_ratio['merge']['q'] *= len(unfrozen)/(num_frozen + 1)
        if not force: self.proposal_ratio['rem_par']['q'] *= len(unfrozen)/(num_frozen + 1)
        #print(f"Prob to freeze = {self.proposal_ratio['rem_par']['q']}")

        if  pid in protomodel.masses: protomodel.masses.pop(pid)
        if  pid in protomodel.decays: protomodel.decays.pop(pid)

        newSSMs = {}
        if reassignPID != None:
            for pids, value in protomodel.ssmultipliers.items():
                if reassignPID[0] == abs ( pids[0] ):
                    newpids = ( reassignPID[1], pids[1] )
                    newSSMs[newpids]=value
                if reassignPID[0] == abs ( pids[1] ):
                    newpids = ( pids[0], reassignPID[1] )
                    newSSMs[newpids]=value

        removeSSM = [pids for pids in protomodel.ssmultipliers if (pid in pids or -pid in pids)]
        for pids in removeSSM:
            protomodel.ssmultipliers.pop(pids)
        for pids,value in newSSMs.items():
            protomodel.ssmultipliers[pids]=value

        #Fix branching ratios and rescale signal strengths, so other channels are not affected
        self.removeIllegalBRs(rescaleSSMs=True, protomodel=protomodel)
        return [ pid ]

    def unFreezeParticles (self, pid : int, force : bool = False,
            protomodel = None, cap_ssm : float = 100.,
            mass : Union[None,float] = None ) -> list:
        """ unfreeze particle pid, assign masses, BRs and signal strength
        multipliers. but also unfreeze sibling particles that are
        forced to be mass degenerate

        :param pid: PID to be unfrozen
        :param force: If False, will only unfreeze the particle if it does not
        violate the canonical order (e.g. will not unfreeze stop2 if stop1 is
        frozen).
        :param protomodel: if given, then unFreeze for that protomodel, not for
        your protomodel
        :param mass: if given, unfreeze to that mass

        :returns: list of particles that really got unfrozen
        """
        nmass = self.unFreezeParticle ( pid, force, protomodel, cap_ssm, mass )
        if nmass == None:
            return []
        unfrozen = [ pid ]
        allpids = self.forcedMassDegeneratePids ( pid, include_self = False )
        if allpids != []:
            self.log ( f"unFreezeParticles {pid}->{allpids}" )
        for ipid in allpids:
            x = self.unFreezeParticle ( ipid, force, protomodel, cap_ssm, nmass )
            if x is not None:
                unfrozen += [ ipid ]
        return unfrozen

    def unFreezeParticle (self, pid : int, force : bool = False,
            protomodel = None, cap_ssm : float = 100.,
            mass : Union[None,float] = None ) -> float:
        """ unfreeze particle pid, assign masses, BRs and signal strength
        multipliers.

        :param pid: PID to be unfrozen
        :param force: If False, will only unfreeze the particle if it does not
        violate the canonical order (e.g. will not unfreeze stop2 if stop1 is
        frozen).
        :param protomodel: if given, then unFreeze for that protomodel, not for
        your protomodel
        :param mass: if given, unfreeze to that mass

        :returns: mass of particle that got unfrozen or None
        """

        if protomodel is None:
            protomodel = self.M

        #Check for canonical ordering.
        frozen = protomodel.frozenParticles()
        n_frozen = len(frozen)
        unfrozen = protomodel.unFrozenParticles()
        n_unfrozen = len(protomodel.unFrozenParticles( withLSP=False ))

        if not force:
            #If pid matches the heavier state and the lighter state is frozen,
            #do not unfreeze the particle
            for pids in self.canonicalOrder:
                if pid == pids[1] and pids[0] in frozen:
                    return None

        #get total num of frozen and unfrozen par for proposal ratio
        for pids in self.canonicalOrder:
            if pids[0] in self.M.environ.forbiddenparticles and pids[1] in self.M.environ.forbiddenparticles: continue
            if pids[0] in frozen and pids[1] in frozen:
                n_frozen -= 1                   #num of par to unfreeze is smaller (i.e cannot unfreeze pids[1] while pids[0] is frozen)
            if pids[0] in unfrozen and pids[1] in unfrozen:
                n_unfrozen -= 1               #num of par to freeze is smaller (i.e cannot freeze pids[0] while pids[1] is unfrozen)

        if self.M.environ.forbiddenparticles != []:
            n_frozen -= len(self.M.environ.forbiddenparticles)
            if n_frozen < 0: self.log(f"Num frozen {n_frozen} < 0! "); n_frozen = 1
        #proposal ratio = p(i+1 -> i)/p(i->i+1) = p(rem pid)/p(add pid) = (1/(n_un+1))/(1/n_fr)
        self.proposal_ratio['add_par']['q'] *= n_frozen/(n_unfrozen + 1)
        #print(f"Prob to unfreeze = {self.proposal_ratio['add_par']['q']}")

        #Absolute mass range:
        maxMass = protomodel.maxMass    #2400 GeV
        minMass = protomodel.masses[LSP]

        #Redefine mass range if necessary to make sure the mass ordering is respected:
        for pids in self.canonicalOrder:
            if pid == pids[0] and not (pids[1] in frozen):
                if pids[1] in protomodel.masses:
                    maxMass = protomodel.masses[pids[1]] #Do not allow for masses above the heavier state
            elif pid == pids[1] and not (pids[0] in frozen):
                if pids[0] in protomodel.masses:
                    minMass = protomodel.masses[pids[0]] #Do not allow for masses below the ligher state

        if pid in self.walledpids:
            ## heed the wall!
            minMass = max ( self.walledpids[pid], minMass )

        # offshell = False
        if pid in [ 1000023, 1000024 ]:
            # for C1 and N2 we want a 10% chance to start in the offshell regime -> increase prob?
            p = np.random.uniform ( 0, 1 )
            if p < 0.1:
                # offshell = True
                self.log ( f"Unfreezing {self.namer.asciiName(pid)}, randomly chose to restrict to offshell mass!" )
                if pid == 1000023: maxMass = minMass + smMasses["Z"] + smWidths["Z"]
                else: maxMass = minMass + smMasses["W"] + smWidths["W"]

        m_random = float(np.random.uniform ( 0., 1. ))
        tmpMass = minMass + (maxMass-minMass)*m_random

        for degeneracy in protomodel.forced_degeneracies:
            if pid in degeneracy:
                for degen_pid in degeneracy:
                    if degen_pid != pid and degen_pid in unfrozen:
                        degen_mass = protomodel.masses[degen_pid]
                        tmpMass = degen_mass
                        if mass is not None and degen_mass != mass:
                            self.warning(f"Unfreeze {pid} at mass {mass} while it should be mass degenerate at {degen_mass} with {degen_pid}.")
                        break

        ctr = 0
        while pid in [ 1000006, 2000006 ] and self.inCorridorRegion ( tmpMass, protomodel.masses[LSP] ):
            # if in corridor region, redraw!
            tmpMass =  minMass + (maxMass-minMass)*m_random
            mstop2 = 2000.
            if 2000006 in protomodel.masses:
                mstop2 = protomodel.masses[2000006]
                protomodel.masses[2000006] = mstop2 + 20.
                if pid == 2000006:
                    tmpMass = protomodel.masses[2000006]
            if pid == 1000006:
                maxMass = mstop2 + 20.

        if mass is not None:
            tmpMass = mass
        protomodel.masses[pid] = tmpMass

        self.record ( f"Unfreeze mass of {pid}({self.namer.texName(pid,addDollars=True)}) to {tmpMass:.1f}" )
        self.log ( f"Unfreeze mass of {self.namer.asciiName(pid)} to {protomodel.masses[pid]:.1f}" )

        # Set branchings
        self.log(f"Initializing Branchings for {self.namer.asciiName(pid)}({pid})")
        initialized = self.initBranchings(pid, protomodel=protomodel)
        if not initialized and not pid in self.decaylessParticles:
            self.proposal_ratio['add_par']['q'] = 1.
            if pid in self.M.decaylessParticles:
                self.log(f"No decays for {self.namer.asciiName(pid)}({pid}) -- but it's marked as decayless")
            else:
                self.log(f"No decays for {self.namer.asciiName(pid)}({pid})")
                self.freezeParticles ( pid, force=True, protomodel=protomodel )
            return None

        #Add pid pair production and associated production to protomodel.ssmultipliers:
        self.log(f"Initializing Production Modes for {self.namer.asciiName(pid)}({pid})")
        self.initSSMFor(pid, protomodel=protomodel, cap_ssm=cap_ssm)

        return tmpMass

    def randomlyChangeMasses ( self, prob = 0.05, dx = 200.0 ):
        """ take a random step in mass space for a single unfrozen particle

        :param prob: Probability for changing the mass
        :param dx: Defines the interval for selecting the delta m (-dx,dx) """

        uMass = np.random.uniform ( 0., 1. )
        if uMass < (1-prob):
            self.log("Not changing masses")
            return 0

        unfrozen = self.M.unFrozenParticles()
        if len(unfrozen)==0:
            self.log("Error: No particles to change masses?")
            return 0

        pid = int(np.random.choice ( unfrozen ))

        #Define mass interval
        maxMass = self.M.maxMass
        minMass = self.M.masses[LSP]
        #In case the pid corresponds to a lighter or heavier state of a pair of particles,
        #make sure the mass ordering is respected:
        for pids in self.canonicalOrder:
            if pid == pids[0] and pids[1] in unfrozen:
                # Do not allow for masses above the heavier state
                maxMass = self.M.masses[pids[1]]
            elif pid == pids[1] and pids[0] in unfrozen:
                # Do not allow for masses below the ligher state
                minMass = self.M.masses[pids[0]]

        #If the particle is the LSP, relax the lower limit
        if pid == LSP:
            minMass = 10.0
            maxMass = 1500.0

        # an artificial wall because the maps are bounded from below
        if abs(pid) in self.walledpids:
            ## heed the wall!
            minMass = max ( self.walledpids[abs(pid)], minMass )

        ret = self.randomlyChangeMassOf ( pid, dx=dx, minMass=minMass, maxMass=maxMass )
        if pid in [ 1000023, 1000024 ] and pid in self.M.unFrozenParticles():
            # for C1 and N2, if one of the two gets changed, have a 10% chance that the other gets set to the same value
            p = np.random.uniform(0,1)
            offshell = self.checkIfOffshell(pid)
            if offshell: p = np.random.uniform(0,0.5)        #SN: check if this makes sense?
            if self.run_mcmc: p = 1.0
            if p < .1:
                mass = self.M.masses[pid]
                otherpid = 1000024 if pid == 1000023 else 1000023
                heavypid = 1000037 if pid == 1000023 else 1000025 #get N3 and C2 masses
                # remember the frozen particles, so we can check if we just unfroze this guy
                were_frozen = self.M.frozenParticles()
                was_offshell = False
                if otherpid not in were_frozen: was_offshell = self.checkIfOffshell(otherpid)
                self.M.masses[otherpid] = float(mass * np.random.uniform ( .99, 1.01 ))
                if heavypid not in were_frozen: #check for canonical ordering
                    if self.M.masses[otherpid] > self.M.masses[heavypid]: self.M.masses[otherpid] = self.M.masses[heavypid] - 10.
                self.log ( f"mass of {self.namer.asciiName(pid)} got changed to {mass:.1f}. hattrick, changing also for {self.namer.asciiName(otherpid)}!" )
                # If the particle was frozen before, we need to unfreeze
                if otherpid in were_frozen:
                    initialized = self.initBranchings(otherpid)
                    if initialized: self.initSSMFor(otherpid)
                #if otherpid was not offshell before but now is offshell and vice versa, initialize branchings
                else:
                    if self.checkIfOffshell(otherpid) != was_offshell: self.initBranchings(otherpid)
                if otherpid in self.M.unFrozenParticles(): #added check since sometimes after initBranchings, total br is 0 and particle is removed
                    self.record ( f"change mass of {self.namer.asciiName(otherpid)} to {self.M.masses[otherpid]}" )
                    ret+=1

        #Fix branching ratios and rescale signal strenghts, so other channels are not affected
        self.removeIllegalBRs(rescaleSSMs=True)

        return ret

    def inCorridorRegion ( self, mstop, mlsp ):
        """ are we in the top corridor region, i.e. (mstop - mlsp) \approx mtop?
            i.e. mstop < 280 and 150 < (mstop - mlsp) < 200
        :returns: true, if in corridor region
        """
        if mstop>280:
            return False
        return 150. < (mstop-mlsp) < 200.

    def randomlyChangeMassOf ( self, pid : int, dx : Union[float,None] = None,
            minMass : Union[float,None] = None,
            maxMass : Union[float,None] = None ) -> int:
        """ randomly change the mass of pid
        :param dx: the delta x to change. If none, then use a model-dependent
                   default
        :param minMass: minimum allowed mass for the particle.
                        If not defined, use the LSP mass.
        :param maxMass: maximum allowed mass for the particle.
                        If not defined, use the protomodel maxMass.

        :returns: number of masses that were changed
        """
        denom = 1.0

        if self.M.TL is not None and self.M.TL > 0:
            #short term fix -> discuss with Wg!
            denom = np.sqrt(self.M.TL) + 1.0

        step_size = 100 if self.M.masses[pid]<1000 else 500
        dx = (step_size)/denom
        if dx < 0.:
            self.highlight ( "info", f"dx={dx}<0. this should not happen. pid={pid} mass={self.M.masses[pid]} denom={denom}" )

        self.log(f"Current mass of {self.namer.asciiName(pid)}({pid}) = {self.M.masses[pid]:.3f} GeV, dx = {dx:.3f} GeV")

        if not minMass:
            minMass = self.M.masses[LSP]
        if not maxMass:
            maxMass = self.M.maxMass

        was_offshell, offshell = self.checkIfOffshell(pid), False
        if pid in [ 1000023, 1000024 ] and not was_offshell:
            # for C1 and N2 we want a 10% chance to move into the offshell region
            p = np.random.uniform ( 0, 1 )
            if self.run_mcmc: p = 1.0        #dont jump from onshell to offshell and vice-versa in mcmc walk
            if p < 0.1:
                offshell = True
                self.log ( f"randomly chose {self.namer.asciiName(pid)} to restrict to offshell mass!" )
                if pid == 1000023: maxMax = minMass + smMasses["Z"] + smWidths["Z"]
                else: maxMax = minMass + smMasses["W"] + smWidths["W"]

        massIsLegal = False
        ctIterations = 0
        while not massIsLegal:
            ctIterations += 1
            massIsLegal = True
            #tmpmass = float(self.M.masses[pid] + np.random.uniform(-dx,dx))
            tmpmass = float(norm.rvs(loc=self.M.masses[pid], scale=dx))
            if offshell: tmpMass = float(np.random.uniform ( minMass, maxMass ))
            # Enforce mass interval:
            if pid in [ 1000006, 2000006 ] and self.inCorridorRegion ( tmpmass, self.M.masses[LSP] ):
                massIsLegal = False
            if pid == LSP and 1000006 in self.M.masses and self.inCorridorRegion ( self.M.masses[1000006], tmpmass ):
                massIsLegal = False
            if pid == LSP and 2000006 in self.M.masses and self.inCorridorRegion ( self.M.masses[2000006], tmpmass ):
                massIsLegal = False
            if tmpmass > maxMass: ## check again if we are legal
                # tmpmass = maxMass-1.0
                # not ok, rerun
                massIsLegal = False
            if tmpmass < minMass: # check again if we are legal
                # not ok, rerun
                # tmpmass = minMass+1.
                massIsLegal = False
            if tmpmass in [ float("nan"), float("inf"), None ]:
                massIsLegal = False
                self.pprint ( f"huh? we have a tmpmass of {self.namer.asciiName(pid)} at {tmpmass}, was at {self.M.masses[pid]}, dx={dx}" )
            dx = dx * 1.2 ## to make sure we always get out of this
            if ctIterations > 20: # seems like we are in a super constrained situation
                self.pprint ( f"huh? we have a tmpmass of {pid} is {tmpmass} was at {self.M.masses[pid]} dx={dx} breaking off after {ctIterations} iterations" )
                tmpmass = self.M.masses[pid]
                break

        if pid == LSP:
            delta_mass = tmpmass - self.M.masses[LSP]
            self.M.masses[pid] = tmpmass
            self.log(f"Randomly changing LSP mass to {tmpmass:.1f}.")
            shifted=self.shiftAllMassesBy(delta_mass, lsp=False)
            return len(shifted)

        allpids = self.forcedMassDegeneratePids(pid)
        nchanges = 0
        self.log(f"Now randomly changing {pid}->{allpids}: masses {self.M.masses[pid]}->{tmpmass}" )

        for ipid in allpids:
            self.M.masses[ipid]=tmpmass
            if ipid in [ 1000023, 1000024 ]:
                if ipid == 1000023: is_offshell = (tmpmass - self.M.masses[LSP]) < ( smMasses["Z"] + smWidths["Z"])
                if ipid == 1000024: is_offshell = (tmpmass - self.M.masses[LSP]) < (smMasses["W"] + smWidths["W"])
                if was_offshell != is_offshell:     #initialize branchings
                    if self.run_mcmc:
                        self.log(f"Jumping from onshell to offshell mass or vice versa during mcmc walk. Not allowed. Dont change mass of {ipid}.")
                    else:
                        self.log ( f"randomly changing mass of {self.namer.asciiName ( ipid )}({ipid}) to {tmpmass:.1f}" )
                        self.record ( f"change mass of {self.namer.texName(ipid,addDollars=True)} to {tmpmass:.1f}" )
                        self.initBranchings(ipid)
                        nchanges += 1
            else:
                self.log ( f"randomly changing mass of {self.namer.asciiName(ipid)}({ipid}) to {tmpmass:.1f}" )
                self.record ( f"change mass of {self.namer.texName(ipid,addDollars=True)} to {tmpmass:.1f}" )
                nchanges += 1

        return nchanges

    def reassignPIDs(self):
        """ Check if a heavier mass eigenstate is present when the lighter one is
        not. If so, reassign the heavier eigenstate to the lighter one."""
        unfrozen = self.M.unFrozenParticles()
        frozen = self.M.frozenParticles()
        for pids in self.canonicalOrder:
            if pids[0] in frozen and pids[1] in unfrozen:
                self.log(f"{self.namer.asciiName(pids[0])} not present but {self.namer.asciiName(pids[1])} present. Reassigning {self.namer.asciiName(pids[1])} to {self.namer.asciiName(pids[0])}")
                self.M.masses[pids[0]] = self.M.masses[pids[1]]
                self.M.masses.pop(pids[1])
                self.M.decays[pids[0]] = self.M.decays[pids[1]]
                self.M.decays.pop(pids[1])

                newssms = {}
                oldssms = self.M.ssmultipliers
                lightpid = int(abs(pids[0])*abs(pids[1])/pids[1])    #get the charge right
                for pidpair, ssm in oldssms.items():
                    newpidpair = [int(lightpid*pid/abs(pid)) if abs(pids[1]) == abs(pid) else pid for pid in pidpair]
                    newpidpair = tuple(sorted(newpidpair))
                    newssms[newpidpair] = ssm
                self.M.ssmultipliers  = newssms
                self.removeIllegalBRs(rescaleSSMs=True)

    def simplifyModel ( self, dm= 200. ):
        """ Try to simplify model, merging pair of candidate particles with similar masses.
        :param dm: Maximum mass difference for merging
        :returns: None, if no mergable particle pair exists, else returns new protomodel
                  with all possible particles merged.
        """


        self.log ( "trying to simplify model" )
        #Make a copy of the model:
        newModel = self.M.copy()
        merged = True
        nMerges = 0
        #Keep merging candidates until no new merge is possible:
        while merged:
            merged = self.mergeParticles(dm,newModel)
            nMerges += merged #Count number of mergers

        if nMerges > 0:
            #Update cross-sections (if needed)
            newModel.getXsecs()
            return newModel
        else:
            return None

    def mergeParticles( self, dm : float = 200,
        protomodel : Union[ProtoModel,None] = None ):
        """ Look for pair of candidates with mass difference smaller than dm and
        merge them. If several particles can be merged, only merge the ones with
        the smallest mass difference. If protomodel is defined merge the particles
        of the given model, else merge particles in self.M

        :param dm: Maximum mass difference for merging
        :param protomodel: ProtoModel to be modified. If None, use self.M

        :returns: False if no merge was performed, else returns True
        """

        #Loop over candidates
        minDMass = dm
        pidA = None
        pidB = None
        if not protomodel:
            protomodel = self.M
        unfrozen = protomodel.unFrozenParticles()
        for pidGroup in self.mergerCandidates:
            pG = sorted(pidGroup) #Make sure the pids are ordered
            for pA,pB in itertools.product(pG,pG):
                if pA >= pB: continue #Only need to check for unique pairings
                if (not pA in unfrozen) or not (pB in unfrozen):
                    continue #Only need to consider unfrozen particles
                dmass = abs(protomodel.masses[pA]-protomodel.masses[pB])
                if dmass < minDMass:
                    minDMass = dmass
                    pidA = pA
                    pidB = pB

        if pidA is not None:
            if self.M.masses[pidA] > self.M.masses[pidB]:
                self.pprint(f"can not merge particles with wrong mass hierarchy ({int(pidA)} > {int(pidB)})")
                return False
            if pidA in [ 1000006, 2000006 ] and pidB in [ 1000006, 2000006 ]:
                ## merging stops. check if we would end up in corridor.
                avgM = self.computeAvgMass ( (pidA,pidB) )
                if self.inCorridorRegion ( avgM, self.M.masses[LSP] ):
                    self.pprint ( "wont merge the stops since we would end up in corridor region!" )
                    return False
            self.merge((pidA,pidB),protomodel)
            return True
        else:
            return False

    def merge ( self, pair : Tuple[int],
        protomodel : Union[ProtoModel,None] = None,
        strategy : str = "avg" ) -> ProtoModel:
        """ merge the particles with pidA and pidB in protomodel.

        :param pair: Pair of particle pids to be merged
        :param protomodel: ProtoModel to be modified. If None, use self.M
        :param strategy: strategy for resultant mass. Avg: take avg of
        two previous masses. Lower: take lower of the two

        :returns: Protomodel with the particles merged
        """

        if not protomodel:
            protomodel = self.M

        n_par_old = len(protomodel.unFrozenParticles())
        n_par_new = n_par_old - 1
        n_decays_old = sum([len(dc.keys())- 1 for par, dc in protomodel.decays.items() if par != LSP ])

        ## Store original decays
        olddecays = {}
        for mpid,decays in protomodel.decays.items():
            olddecays[mpid] = dict([[dpids,br] for dpids,br in decays.items()])
        ## Store orignal xsecs (needed for rescaling the SSMs)
        oldxsecs = None
        tmpx = protomodel.getXsecs()
        if len(tmpx)>0:
            oldxsecs = tmpx[0]

        pair = list(pair)
        pair.sort()
        p1,p2 = pair[0], pair[1]
        self.log(f"Merging {self.namer.asciiName(p1)} and {self.namer.asciiName(p2)}")
        self.log(f"Masses before merger: {protomodel.masses[p1]:.2f}, {protomodel.masses[p2]:.2f}")
        strategy = strategy.lower()
        assert strategy in [ "avg", "lower" ], "strategy has to be one of: avg, lower"
        if strategy == "avg":
            avgM = self.computeAvgMass ( pair )
            self.log(f"Avg mass for {str(pair)} is {avgM:.2f}")
            protomodel.masses[ p1 ] = avgM ## set this one to the avg mass
        else:
            m = min ( [ self.M.masses[x] for x in pair ] )

        #Get p2 decays:
        p2decays = protomodel.decays[p2]
        #Get allowed decay channels for p1:
        openChannels = self.M.getOpenChannels(p1)
        ## add the decays from pid2 to pid1 if decay is allowed:
        for pids,br in p2decays.items():
            if not pids in openChannels:
                continue
            if not p1 in protomodel.decays:
                self.log(f"how is {p1} not in {protomodel.decays}?" )
                continue
            if pids in protomodel.decays[p1]:   #add to existing br
                if br > 0.001:
                    self.log(f"Add to decays {p1}/{pids}: {br:.2f}")
                protomodel.decays[p1][pids] += br
            else:                               #add new decay mode
                self.log(f"Set decays of {p1}/{pids} to {br:.2f}")
                protomodel.decays[p1][pids] = br

        self.log(f"Normalize branchings of {self.namer.asciiName(p1)} after merge" )
        self.normalizeBranchings ( p1, protomodel=protomodel )

        #Now replace all decays to p2 by decays to p1
        #(since the new (average) p1 mass is always smaller than the p2 mass,
        #there is no chance of running into offshell decays)
        for mpid,decays in olddecays.items():
            for dpids,br in decays.items():
                newpids = dpids
                #Replace any appearence of p2/-p2 in decays by p1/-p1:
                if isinstance(dpids,(list, tuple)):
                    newpids = [dpid if abs(dpid) != abs(p2) else abs(p1)*dpid/abs(dpid) for dpid in dpids]
                    newpids = tuple(newpids)
                elif isinstance(dpids,int) and abs(dpids) == abs(p2):
                    newpids = abs(p1)*dpids/abs(dpids)

                #If original channel did not contain p2, do nothing
                if newpids == dpids:
                    continue
                #Print log message for non-negligible BRs:s
                if br > 0.0001:
                    self.log ( f"redirecting decay of {int(mpid)} from {dpids} to {newpids}: br={br:.2f}" )

                #If new channel was already present, simply add to BR:
                if newpids in protomodel.decays[mpid]:
                    br += protomodel.decays[mpid][newpids]

                #Remove original decay:
                protomodel.decays[mpid].pop(dpids)
                #Add new channel:
                protomodel.decays[mpid][newpids]=br

        n_decays_new = sum([len(dc.keys())- 1 for par, dc in protomodel.decays.items() if par != LSP ])

        if oldxsecs != None:
            ## merge the signal strength multipliers:
            n_ssms_old, n_ssms_new = self.mergeSSMs( pair,
                    oldXsecs = oldxsecs, protomodel=protomodel )

        #Get proposal ratio for merge move
        a,b,c = 2,4,8
        shift_parameter = 10
        z_old = (n_par_old - shift_parameter)/a + n_decays_old/b + n_ssms_old/c
        z_new = (n_par_new - shift_parameter)/a + n_decays_new/b + n_ssms_new/c

        prob_rem = 1.0          #Removing degrees of freedom
        prob_add = (1 + np.exp(z_new))/(1 + np.exp(z_old))
        self.proposal_ratio['merge'] = {'q':prob_add}
        #print(f"prob  merge = {prob_add}")

        mpair = pair
        if mpair[0] < mpair[1]:
            mpair = [ pair[1], pair[0] ]
        ## finally freeze p2:
        self.freezeParticles(p2,protomodel=protomodel,merge=True,
                reassignPID = mpair )
        self.reassignPIDs()
        return protomodel

    def computeAvgMass ( self, pids ):
        """ compute the average mass
        :param merge_strategy: allow for different ways to merge
        :returns: mass, as scalar, in GeV
        """
        ret=0.
        for pid in pids:
            ret+=self.M.masses[pid]
        return ret / len(pids)

    def mergeSSMs ( self, pair : Tuple[int], oldXsecs : list,
            protomodel : Union[ProtoModel,None] = None ):
        """ merge signal strength multipliers for particles in pair.
        The cross-selections involving the merged particles are assumed to be
        added and the corresponding signal strengths are rescaled.

        :param pair: pair of particle PIDs being merged
        :param oldXsecs: cross-sections before the merge
        :param protomodel: protomodel to be modified. If not defined, use self.M

        :returns: tuple(number of old xsecs, number of new xsecs)
        """

        if not protomodel:
            protomodel = self.M

        #Get updated list of unfrozen particles
        unfrozen = protomodel.unFrozenParticles()
        pair = list(pair)
        pair.sort()
        p1,p2 = pair[0], pair[1]

        #Build dictionary with original cross-sections
        #(only select LO cross-sections at 13 TeV)
        oldxsecDict = dict([[xsec.pid,xsec.value.asNumber(fb)] for xsec in oldXsecs
                            if xsec.info.sqrts > 10.*TeV and xsec.info.order <= 0])

        #Find cross-sections PIDs containing p2 or -p2
        #and build the new PIDs (with p2 replaced by p1)
        p2Xsecs = {}
        procDict = {} ## this shows the mapping of production modes
        prodModes = protomodel.getAllowedProdModes()
        for pids,xsec in oldxsecDict.items():
            if not p2 in pids and not -p2 in pids:
                continue
            newpids = [pid if abs(pid) != abs(p2) else int(abs(p1)*pid/abs(pid)) for pid in pids ]
            newpids = tuple(sorted(newpids))
            #Skip processes containing frozen particles:
            if not all([abs(pid) in unfrozen for pid in newpids]):
                continue
            #Skip prod modes which are not allowed
            if newpids not in prodModes:continue
            p2Xsecs[pids] = xsec
            procDict[pids] = newpids

        #Now compute the new SSMs assuming that the cross-sections will be added:
        newSSMs= {}
        for oldpids,newpids in procDict.items():
           #Get xsec value for the process containing p2:
           value = p2Xsecs[oldpids]
           oldvalue = None
           #Check if the new process (with p2->p1) already existed
           if newpids in oldxsecDict:
               oldvalue = oldxsecDict[newpids]
               value += oldvalue #Combine xsec values

           #If the new process (with p2->p1) already existed, rescale SSM:
           if oldvalue:
               #Check if the SSM already existed (if not, take SSM = 1.0)
               oldssm = 1.0
               if newpids in protomodel.ssmultipliers:
                   oldssm = protomodel.ssmultipliers[newpids]
               #The new SSM is going to be the ratio of old and new cross-sections times the old SSM:
               if not newpids in newSSMs:
                   newSSMs[newpids] = 0
               newSSMs[newpids] += oldssm*(value/oldvalue)
           #If the new process does not exist take the SSM for the (old) process containing p2
           else:
               oldssm = 1.0
               if oldpids in protomodel.ssmultipliers:
                   oldssm = protomodel.ssmultipliers[oldpids]
               if not newpids in newSSMs:
                   newSSMs[newpids] = 0.
                   self.ssmultipliers[newpids] = 1.0
               # we dont know what the oldssm would roughly
               # correspond to with the newpids, so we use oldssm for newpids
               # FIXME we can improve here
               oldxsec = p2Xsecs[oldpids]
               newssm, newxsec = self.computeNewSSM ( newpids, oldpids, oldssm, p2Xsecs[oldpids] )
               newSSMs[newpids] += newssm

        #Now replace the SSMs in protomodel:
        for pid,ssm in newSSMs.items():
           protomodel.ssmultipliers[pid] = ssm

        n_ssm_old, n_ssm_new = len(oldxsecDict.keys()), len(newSSMs.keys())
        return n_ssm_old, n_ssm_new

    def computeNewSSM ( self, newpids, oldpids, oldssm, oldxsec ):
        """ ok we need the translation of oldssm to newpids,
        but we are missing an xsec """
        tmpSLHA = tempfile.mktemp( prefix=f".{self.walkerid}_xsecfile",
                                   suffix=".slha",dir=self.SLHATEMPDIR )
        tmpSLHA = self.createSLHAFile(tmpSLHA, addXsecs = False)
        self.computer.compute( 13, tmpSLHA, ssmultipliers = self.ssmultipliers )
        for xsec in self.computer.xsecs:
            if xsec.info.sqrts == 13*TeV and xsec.pid == newpids:
                newxsec = xsec.value.asNumber(fb)
                newssm = oldxsec / newxsec
                if os.path.exists ( tmpSLHA ): ## remove
                    os.remove( tmpSLHA )
                return newssm, newssm*newxsec
        return None, None

    def simplifyMasses ( self ):
        """ return the masses only of the unfrozen particles """
        ret ={}
        unfrozen = self.M.unFrozenParticles()
        for pid in unfrozen:
            ret[pid]=self.M.masses[pid]
        return ret

    def simplifyXSecs ( self, fbmin=.001*fb ):
        """ return the xsecs above a threshold only """

        xsecs={ 8:{}, 13:{} }
        modelXSecs = self.M.getXsecs()[0]
        for xsec in modelXSecs:
            if xsec.value < fbmin:
                continue
            sqrts = xsec.info.sqrts.asNumber(TeV)
            if not xsec.pid in xsecs[sqrts]:
                xsecs[sqrts][xsec.pid]=xsec
            else:
                if xsecs[sqrts][xsec.pid].info.order < xsec.info.order:
                    xsecs[sqrts][xsec.pid]=xsec

        return xsecs

    def printXSecs ( self, fbmin=.001*fb, useParticleNames : bool = False ):
        """ print the cross sections in a human-readable way
        :param useParticleNames: if true, use names for particles not pids
        """
        xsecs = self.simplifyXSecs( fbmin )
        for sqrts in xsecs.keys():
            pidss = list ( xsecs[sqrts].keys() ) # list of list of pids
            pidss.sort()
            print ( f"{sqrts} TeV:" )
            for pids in pidss:
                mpids = tuple ( [ x for x in pids if x != None ] )
                if len(mpids)==1:
                    mpids = mpids[0]
                if useParticleNames:
                    mpids = self.namer.asciiName ( mpids )
                else:
                    mpids = f"{str(mpids):>22s}"
                xsec = xsecs[sqrts][pids]
                label = ""
                comment = ""
                if hasattr ( xsec, "comment" ):
                    comment = xsec.comment
                if "dict" in xsec.info.label:
                    label = " (from dict)"
                print ( f" {mpids}: {xsec.value.asNumber(fb):.2f} fb{label} {comment}" )

    def simplifyDecays ( self, protomodel=None):
        """ return the decays only of the unfrozen particles,
            only != 0 """

        if protomodel is None:
            protomodel = self.M
        ret ={}
        unfrozen = protomodel.unFrozenParticles()
        for mpid,decays in protomodel.decays.items():
            if mpid not in unfrozen:
                continue
            d = {}
            for dpid,dbr in decays.items():
                if dbr > 1e-5:
                    d[dpid]=dbr
            ret[mpid]=d
        return ret

    def allXSecsAbove ( self, threshold=.01*fb, sqrts=13*TeV, order=LO ):
        """ return list of all cross sections above threshold.
        :returns: list of tuples of pids, cross sections (that had the SSM applied),
                          and SSMs that *were* applied.
        """
        if type(threshold)==float and threshold>0.:
            self.pprint ( "note: interpreting threshold as fb" )
            threshold = threshold * fb
        ret = []
        modelXSecs = self.M.getXsecs()[0]
        for xsec in modelXSecs:
            if xsec.info.order != order:
                continue
            if abs (( xsec.info.sqrts - sqrts ).asNumber(TeV)) > .1:
                continue
            xs = xsec.value
            ssm = 1.
            if xsec.pid in self.M.ssmultipliers:
                ssm = self.M.ssmultipliers[xsec.pid]
            ret.append ( (xsec.pid, xs, ssm) )
        ret.sort( key = lambda x: x[1], reverse = True )
        return ret

    """
    def xsecsFor ( self, pids : Tuple , sqrts=13*TeV, min_order=LO ) -> Tuple:
        # return the cross sections for pids.

        :param pids: tuple of two pids
        :param min_order: report only if order is above or equal min_order
        :returns: cross section (that had the SSM applied),
        and SSM that *was* applied.
        #

        ssm = 1.
        if pids[1] < pids[0]:
            pids = ( pids[1], pids[0] )
        if pids in self.M.ssmultipliers:
            ssm = self.M.ssmultipliers[pids]
        xs = 0. * fb
        modelXSecs = self.M.getXsecs()[0]
        order = -1
        for xsec in modelXSecs:
            if xsec.info.order < order:
                continue
            if abs ( ( xsec.info.sqrts - sqrts ).asNumber(TeV) ) > .1:
                continue
            if xsec.pid != pids:
                continue
            # we want the highest order
            if xsec.info.order > order:
                xs = xsec.value
                order = xsec.info.order
        return xs,ssm
    """

    def simplifySSMs ( self, removeOnes=False, removeZeroes=False,
                       threshold=0.001*fb, store = False ):
        """ return only SSMs for unfrozen particles
        :param removeOnes: if True, remove ssms == 1.
        :param removeZeroes: if True, remove ssms == 0.
        :param threshold: remove the SSMs for cross sections smaller
                                          than the given threshold (13TeV, LO).
        :param store: if True, overwrite original ssms with ours
        :returns: dictionary of SSMs
        """
        if type(threshold)==float and threshold>0.:
            self.pprint ( "note: interpreting threshold as fb" )
            threshold = threshold * fb
        ret = {}
        frozen = self.M.frozenParticles()
        modelXSecs = self.M.getXsecs()[0]

        for pids,v in self.M.ssmultipliers.items():
            if removeOnes and abs(v-1.)<1e-5:
                continue
            if removeZeroes and v<1e-7:
                continue
            xsecBigEnough = False
            if threshold > 0.*fb:
                for xsec in modelXSecs:
                    if xsec.info.sqrts.asNumber(TeV)<10:
                        continue
                    if pids == xsec.pid: # they are always sorted
                        sigma = xsec.value
                        if sigma > threshold:
                            xsecBigEnough = True
            if not xsecBigEnough:
                continue
            isFrozen = False
            for pid in pids:
                if pid in frozen or -pid in frozen:
                    isFrozen = True
            if isFrozen:
                continue
            ret[pids]=v
        if store:
            self.M.ssmultipliers = ret
        return ret

    def getSSMsNotInBestCombo ( self ):
        """Freeze Production Modes which do not contribute to the best combination"""
        if not self.M.bestCombo: #no combination
            return None

        nfrozen = 0
        prod_modes = []
        from tester.combiner import Combiner
        c = Combiner ( self.walkerid )
        prodmodes_combo = c.getAllSSMsOfCombo ( self.M.bestCombo )
        prodmodes_model = list(self.M.ssmultipliers.keys())[:]
        for pm in prodmodes_model:
            if pm not in prodmodes_combo:
                self.log(f"Not rescaling {pm} production as it does not contribute to the best combination.")
                prod_modes.append(pm)
                '''
                #Get proposal ratio for freezing ssms:
                #q(rem) = 1.0 since we force removal
                #q(add) get from function z_k
                #Total proposal ratio for move = q(add)*(0.7/len(allowed production modes)) (See how we add new production modes)
                q_12, q_21 = self.z_model(old_protomodel, self.M, force_move=True)
                len_prod_modes = len(self.M.getAllowedProdModes())
                q_mov = (q_21/q_12)*(0.7/len_prod_modes)
                self.proposal_ratio['q_total'] *= q_mov
                nfrozen += 1
                '''

        return prod_modes

    def getAllPidsOfBestCombo ( self ) -> Set:
        """ get all the particle ids of BSM particles in
        the best combo. """
        if not self.M.bestCombo: #no combination
            return None
        from tester.combiner import Combiner
        c = Combiner ( self.walkerid )
        return c.getAllPidsOfCombo ( self.M.bestCombo )

    def freezePidsNotInBestCombo ( self ):
        """ Freeze pids that are present in the protomodel, but do not contribute to the best combination"""

        okPids = self.getAllPidsOfBestCombo()
        if len(okPids)==0: ## means we dont have a best combo
            return 0
        unfrozen = self.M.unFrozenParticles( withLSP=False )
        import copy
        olddecays = copy.deepcopy(self.M.decays)
        nfrozen = 0
        for pid in unfrozen:
            if not pid in okPids:
                freeze = True
                #pid could be in the decay of another particle contributing to the bestCombo
                for par,dec in olddecays.items():
                    if par == pid: continue
                    dpids = dec.keys()
                    pidpresent = [True if pid in dp else False for dp in dpids]
                    if True in pidpresent:
                        self.log(f"{self.namer.asciiName(pid)} not in bestCombo but in decay {self.namer.asciiName(par)} in bestCombo. Not taking out {self.namer.asciiName(pid)}.")
                        freeze = False
                        break
                if freeze:
                    self.log(f"{self.namer.asciiName(pid)} does not contribute to bestCombo. Taking out {self.namer.asciiName(pid)}.")
                    old_protomodel = self.M.copy()
                    frozen_pids = self.freezeParticles ( pid, force=True )
                    if frozen_pids: nfrozen += len(frozen_pids)

        return nfrozen

    def backupModel ( self ):
        """ backup the current state """
        print(self.M)
        self._backup = { "llhd": self.M.llhd, "letters": self.M.letters, "TL": self.M.TL,
                         "K": self.M.K, "muhat": self.M.muhat,
                         "description": self.M.description,
                         "ul_critic_tpList": copy.deepcopy(self.M.ul_critic_tpList),
                         "bestCombo": copy.deepcopy(self.M.bestCombo),
                         "masses": copy.deepcopy(self.M.masses),
                         "ssmultipliers": copy.deepcopy(self.M.ssmultipliers),
                         "decays": copy.deepcopy(self.M.decays),
                         "environ": copy.deepcopy(self.M.environ),
                         "rvalues": copy.deepcopy(self.M.rvalues),
                         "_stored_xsecs" : copy.deepcopy(self.M._stored_xsecs),
                         "_xsecMasses" : copy.deepcopy(self.M._xsecMasses),
                         "_xsecSSMs" : copy.deepcopy(self.M._xsecSSMs),
                        }
        if hasattr ( self.M, "ul_critic" ): self._backup["ul_critic"]=self.M.ul_critic
        if hasattr ( self.M, "llhd_critic" ): self._backup["llhd_critic"]=self.M.llhd_critic


    def restoreModel ( self, reportReversion : bool = False ):
        """ restore from the backup

        :param reportReversion: if true, then call "self.record"
        """
        if not hasattr ( self, "_backup" ):
            raise Exception ( "no backup available" )
        if reportReversion:
            self.record ( "revert step" )
        for k,v in self._backup.items(): ## do not!! shallow copy here
            setattr ( self.M, k, copy.deepcopy(v) )

    def delBackup ( self ):
        """ delete protomodel backup dictionary"""
        # if all and hasattr ( self, "_backup" ):
        if hasattr ( self, "_backup" ):
            del self._backup

    def __getattr__(self, name : str ):
        """ delegate everything else to the protomodel you own """
        return getattr(self.M, name )

    def __dir__(self):
        """ for delegation, to make tab completion and introspection work """
        return list(set(super().__dir__() + dir(self.M)))

if __name__ == "__main__":
    import pickle
    f=open("hiscores.cache","rb" )
    protomodels = pickle.load(f)
    f.close()
    ma = Manipulator ( protomodels[0], verbose=True )
    print ( ma.getAllPidsOfBestCombo() )
    #ma.merge ( ( 1000001, 1000003 ), force_merge = True )
    #import IPython
    #IPython.embed()
