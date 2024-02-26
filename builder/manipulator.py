#!/usr/bin/env python3

""" Class that encapsulates the manipulations we perform on the protomodels,
    so that the protomodel class is a data-centric class, and this one
    an algorithm-centric class.
"""

__all__ = [ "Manipulator" ]

from ptools.sparticleNames import SParticleNames
from builder.protomodel import ProtoModel
from smodels.tools.physicsUnits import fb, TeV
from smodels.theory.crossSection import LO
from smodels.theory.theoryPrediction import TheoryPrediction
import smodels
import copy, numpy, time, os, sys, itertools, colorama, random
from colorama import Fore as ansi
from typing import Union, Dict, List, Tuple
from builder.loggerbase import LoggerBase
from os import PathLike

class Manipulator ( LoggerBase ):
    """ contains the protomodel manipulation algorithms. """

    # walledpids are particle ids that have a minimum mass requirement,
    # given as the values in the dictionary
    walledpids = { 1000001 : 310, 1000002 : 310, 1000003 : 310, 1000004 : 310, 
                   1000021 : 310, 1000023 : 20, 1000024 : 100, 1000037 : 100,
                   1000025 : 20 }
    ## forbiddenparticles are particle ids that we do not touch in this run
    forbiddenparticles = []

    def __init__ ( self, protomodel : Union[ProtoModel,Dict,PathLike],
            strategy: str = "aggressive", verbose : bool = False,
            do_record : bool = False, seed : Union[bool,int] = None,
            nth : int = 0 ):
        """
        :param protomodel: is either a protomodel, or a hiscore dictionary,
        or a path to a protomodel
        :param strategy: what combination strategy?
        currently we only have "aggressive"
        :param do_record: do record actions taken
        :param seed: random seed
        :param nth: if initialisation from hiscores dict file, initialise
        from nth entry in that dict file

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
        super(Manipulator, self).__init__ ( 0 )
        if type ( protomodel ) == ProtoModel:
            ## make sure we log correctly asap
            self.walkerid = protomodel.walkerid
        self.namer = SParticleNames ( False )
        self.M = protomodel
        if type(protomodel) == dict:
            self.M = ProtoModel ( )
            self.initFromDict ( protomodel )
        if type(protomodel) == str:
            self.M = ProtoModel ( )
            if protomodel.endswith ( ".dict" ):
                self.initFromDictFile ( protomodel, nth = nth )
        self.walkerid = self.M.walkerid
        self.seed = seed
        self.strategy = strategy
        self.verbose = verbose
        #Store a canonical order for the masses. So the ordering in each
        #tuple is enforced
        self.canonicalOrder =  [ ( 1000006, 2000006 ), ( 1000005, 2000005 ),
                          ( 1000023, 1000025 ), ( 1000024, 1000037 ),
                          ( 1000025, 1000035 ), ( 1000023, 1000025 ),
                          ( 1000001, 1000002 ), (1000001, 1000003) ]
        #Define groups of particles to be merged if their mass difference is below
        #a given mass gap
        self.mergerCandidates =   [ (1000001, 1000002, 1000003, 1000004 ),
                                    ( 1000005, 2000005 ), (1000006, 2000006),
                                    ( 1000024, 1000037 ), ( 1000023, 1000025 )  ]
        self.do_record = do_record ## if True, then record changes
        self.recording = [] ## just in case

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
            self.pprint ( "when trying to teleport, found %s. cancel teleportation." % e )
            # can happen if it is just being written. in this case dont teleport
            return
        ith = 0
        choices = []
        f = 1
        for i in range(len(dicts)-1,-1,-1):
            choices += [i]*f
            f=f*2
        ith = random.choice ( choices )
        self.pprint ( "teleporting, we have %d dicts" % len(dicts) )
        self.pprint ( "choosing the %dth entry, it has a K of %.2f" % \
                      ( ith, dicts[ith]["K"] ) )
        step = self.M.step
        nth = "%dth" % ith
        longforms = { 1: "first", 2: "second", 3: "third" }
        if ith in longforms:
            nth = longforms[ith]
        self.record ( "teleporting to %s hiscore" % nth )
        self.initFromDict ( dicts[ith], initTestStats=True )
        self.M.step = step ## continue counting!
        self.M.bestCombo = None

    def shiftAllMassesBy ( self, dm : float ):
        """ shift the masses of all unfrozen particles by dm [GeV]
        (simple covenience function)

        :param dm: the shift of all masses, in GeV
        """
        for pid,m in self.M.masses.items():
            self.M.masses[pid]=m+dm

    def writeDictFile ( self, outfile : Union[str,None] = "pmodel.dict",
            cleanOut : bool = True, comment : str = "", appendMode : bool = False,
            ndecimals : int = 6 ) -> Dict:
        """ write out the dict file to outfile

        :param outfile: output file, but replacing %t with int(time.time()). If None,
                        then dont write file, just create dictionary object
        :param cleanOut: clean the dictionary from defaults
        :param comment: add a comment field
        :param ndecimals: number of digits after decimal when rounding
        :param appendMode: if true, append to file, and add comma after dictionary.
                           if false, overwrite, and no comma at the end.
        :returns: the dictionary with the protomodel
        """
        if type(self.M) == type(None):
            ## there is nothing to write
            return
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
        import time
        if hasattr ( self, "seed" ) and self.seed != None:
            D["seed"]=self.seed
        D["timestamp"]=time.asctime()
        D["Z"]=round(self.M.Z,ndecimals)
        D["K"]=round(self.M.K,ndecimals)
        D["walkerid"]=self.M.walkerid
        D["step"]=self.M.step
        D["codever"]=self.M.codeversion
        D["smodelsver"]=smodels.installation.version()
        D["dbver"]=self.M.dbversion
        D["description"]=self.M.description
        if hasattr ( self.M, "critic_description" ):
            D["critic_description"]=self.M.critic_description
        if len(comment)>0:
            D["comment"]=comment
        if outfile == None:
            return D
        fname = outfile.replace("%t", str(int(time.time())) )
        if not appendMode:
            self.pprint ( f"writing model to {fname}" )
        mode,comma = "wt",""
        if appendMode:
            mode,comma = "at",","
        with open ( fname, mode ) as f:
            f.write ( f"{D}{comma}\n" )
            f.close()
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
        :param initTestStats: if True, set also test statistics K and Z
        :param nth: if we find a list of models, pick the nth. 0 = 1st. If nth
                    does not exist, return False
        :returns: true, if successful
        """
        if not os.path.exists ( filename ):
            self.pprint ( f"filename {filename} does not exist!" )
            return False
        with open ( filename, "rt" ) as f:
            D = eval ( f.read() )
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


    def initFromDict ( self, D : Dict, filename : str = "",
            initTestStats : bool = False ):
        """ setup the protomodel from dictionary D.
        :param D: dictionary, as defined in pmodel*.dict files.
        :param filename: name of origin. not necessary, only for logging.
        :param initTestStats: if True, set also test statistics K and Z
        """
        scom = ""
        if "comment" in D:
                scom = ": " + D["comment"]
        if filename == "":
            line = "initializing from dictionary: "
            for k,v in D["masses"].items():
                line += self.namer.asciiName(k)+ ", "
            line = line[:-2]
            self.pprint ( line )
        else:
            self.highlight ( "info", f"starting with {os.getcwd()}/{filename}{scom}" )
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
        if "walkerid" in D:
            self.M.walkerid = D["walkerid"]
        if initTestStats:
            if "Z" in D:
                self.M.Z = D["Z"]
            if "K" in D:
                self.M.K = D["K"]
        if "xsecs[fb]" in D:
            tmp = D["xsecs[fb]"]
            xsecs = []
            from smodels.theory.crossSection import XSection
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

    def cheat ( self, mode = 0 ):
        ## cheating, i.e. starting with models that are known to work well
        if mode == 0: ## no cheating
            return
        filename = f"pmodel{mode}.dict"
        if not os.path.exists ( filename ):
            self.highlight ( "red", f"cheat mode {mode} started, but no {os.getcwd()}/{filename} found" )
            sys.exit(-1)
        # scom = ""
        with open ( filename, "rt" ) as f:
            m = eval ( f.read() )
        self.initFromDict ( m, filename )

    def checkForNans ( self ):
        """ check protomodel for NaNs, for debugging only """
        for pid,m in self.M.masses.items():
            if numpy.isnan ( m ):
                self.pprint ( "checking for nans: mass of %d is nan" % pid )

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
            print ( f" `- {ansi.GREEN}{i.analysisId()}:{dId}: {txns}{ansi.RESET}" )
            line = "               "
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
                        line += f" Z={ansi.RED}{(dI.observedN - eBG ) / toterr :.1f}*sigma{ansi.RESET}"
                    print ( line )
                if dtype in [ "upperLimit", "combined" ]:
                    eUL = i.getUpperLimit ( expected = True ).asNumber(fb)
                    oUL = i.getUpperLimit ( expected = False ).asNumber(fb)
                    sigma_exp = eUL / 1.96 # the expected scale, sigma
                    Z = ( oUL - eUL ) / sigma_exp
                    line += f"obs={oUL:.1f}*fb exp={eUL:.1f}*fb Z={ansi.RED}{Z:.1f}*sigma{ansi.RESET}"
                    print ( line )

            for pids in i.PIDs[:2]:
                s = str(pids)
                if len(s) > 80:
                    s=s[:76]+" ..."
                print ( f"              {s}" )
            if len(i.PIDs)>3:
                print ( "               ..." )

    def printAllTheoryPredictions ( self, detailed : bool = False ):
        """ pretty print all theory predictions for the model 
        :param detailed: if true, give more details
        """
        print ( "theory predictions:" )
        combo = self.M.tpList
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
            for pids in i.PIDs[:2]:
                s = f"{self.namer.asciiName(pids)}"
                if len(s) > 80:
                    s=s[:76]+" ..."
                print ( f"              {s}" )
            if len(i.PIDs)>3:
                print ( "               ..." )

    def removeAllOffshell ( self, rescaleSSMs=False, protomodel = None ):
        """ remove all offshell decays and decays of frozen particles. Renormalize all branchings """

        if protomodel is None:
            protomodel = self.M

        for pid in protomodel.frozenParticles():
            if pid in protomodel.decays:
                protomodel.decays.pop(pid)

        #Loop over all decays:
        for pid in protomodel.decays:
            #Get allowed decay channels:
            openChannels = protomodel.getOpenChannels(pid)
            #Check if any of the existing decays are forbidden:
            delDecays = [dpid for dpid in protomodel.decays[pid]
                            if not  dpid in openChannels]
            for dpid in delDecays:
                protomodel.decays[pid].pop(dpid)

            #Make sure to normalize the branchings
            self.normalizeBranchings(pid, rescaleSSMs=rescaleSSMs,
                                        protomodel=protomodel)

    def initBranchings ( self, pid, protomodel=None):
        """ Intialize BRs to diffrent open decay channels for pid.
            Either assign 'democratic' BRs or assign random BRs.
        """

        if not protomodel:
            protomodel = self.M

        #Do not modify the LSP decays
        if pid == protomodel.LSP:
            return

        #Erase BRs (if any has been stored) for offshell too?
        protomodel.decays[pid] = {}
        
        #Get the allowed decay channels:
        openChannels = self.M.getOpenChannels(pid)
        dkeys = set()
        for dpid in openChannels:
            dk = self.M.decay_keys[pid][dpid]
            dkeys.add(dk)
        dkeys = list(dkeys)
         
        
        nitems = len(openChannels)
        
        offshell = False
        if pid == 1000024 and (protomodel.masses[pid] - protomodel.masses[1000022]) < 80.377: offshell = True
        elif pid == 1000023 and (protomodel.masses[pid] - protomodel.masses[1000022]) < 91.1876: offshell = True
        
        for dk in dkeys:
            decay_chan = [key for key,value in self.M.decay_keys[pid].items() if value == dk]
            br = random.gauss ( 1. / nitems, numpy.sqrt ( .5 / nitems )  )
            br = max ( 0., br )
            
            for dpid in decay_chan:
                if offshell:
                    if pid == 1000023:
                        if '11' in dk: protomodel.decays[pid][dpid] = 1.0/18.0
                        elif '2' in dk: protomodel.decays[pid][dpid] = 12.0/18.0
                        elif '5' in dk: protomodel.decays[pid][dpid] = 3.0/18.0
                        else: protomodel.decays[pid][dpid] = 0.0
                    elif pid == 1000024:
                        if '11' in dk: protomodel.decays[pid][dpid] = 1.0/9.0
                        elif '2' in dk: protomodel.decays[pid][dpid] = 6.0/9.0
                        else: protomodel.decays[pid][dpid] = 0.0
                else:
                    if len(dpid) == 3: protomodel.decays[pid][dpid] = 0.0  #turn off offshell 3 body decays
                    else: protomodel.decays[pid][dpid] = br

        #Make sure there is at least one open channel:
        BRtot = sum(self.M.decays[pid].values())
        if BRtot == 0.0 and len(openChannels)>0:
            dk = random.choice(dkeys)
            decay_chan = [key for key,value in self.M.decay_keys[pid].items() if value == dk]
            br = 1.0/len(decay_chan)
            self.M.decays[pid] = {}
            for dpid in decay_chan:
                self.record ( f"change decay of {self.namer.texName(pid,addDollars=True)} -> {self.namer.texName(dpid,addDollars=True)} to {br:.2f}" )
                self.log ( f"changed decay of {self.namer.asciiName(pid)} -> {self.namer.asciiName(dpid)} to {br:.2f}" )
                self.M.decays[pid].update({dpid: br})
            BRtot = 1.0

        #Make sure to normalize the branchings:
        self.normalizeBranchings(pid, protomodel=protomodel)

    def normalizeBranchings ( self, pid, rescaleSSMs=False, protomodel=None ):
        """ normalize branchings of a particle if the total BR is differs from 1.0.

        :param pid: Particle to have their branchings normalized. If pid = None, normalize all decays.
        :param rescaleSSMs: if True, rescale the corresponding signal strength multipliers,
                         so that sigma x br stays the same.
        """

        if not protomodel:
            protomodel = self.M

        if not pid in protomodel.decays:
            protomodel.pprint ( "when attempting to normalize: %d not in decays" % pid )
            return

        BRtot = sum(protomodel.decays[pid].values())
        if BRtot == 0:
            if pid != protomodel.LSP:
                protomodel.pprint ( "when attempting to normalize: total BR (%d) is zero" % pid )
            return

        if abs(BRtot-1.0) < 1e-4:
            #BRs are already normalized.
            return
        self.log ( f"normalized branchings of {self.namer.asciiName(pid)} by {BRtot:.2f}" )

        for dpid in protomodel.decays[pid]:
            protomodel.decays[pid][dpid] *= 1/BRtot

        ## adjust the signal strength multipliers to keep everything else
        ## as it was
        if not rescaleSSMs:
            return

        for pidpair,ssm in protomodel.ssmultipliers.items():
            if pidpair in [ (pid,pid),(-pid,-pid),(-pid,pid),(pid,-pid) ]:
                newssm = min(1e5,ssm*BRtot*BRtot) #Rescale pair production by BRtot^2
            elif (pid in pidpair) or (-pid in pidpair):
                newssm = min(1e5,ssm*BRtot) #Rescale associated production by BRtot
            else:
                continue
            protomodel.ssmultipliers[pidpair]=newssm

    def initSSMFor(self, pid):
        """ Set SSM multipliers to 1 (for pair production of particle/anti-particle):
        """
        unfrozen = self.M.unFrozenParticles()
        pBlist = [pid]
        if self.M.hasAntiParticle(pid):
            pBlist.append(-pid)
        for pidpair in itertools.product(pBlist,pBlist):
            ppair = tuple(sorted(pidpair))
            if not ppair in self.M.ssmultipliers:
                # self.record ( f"change ssm of {self.namer.texName(ppair,addDollars=True)} to 1.0" )
                self.M.ssmultipliers[ppair] = 1.0

        #Set SSM multipliers to 1 for all associated productions with pid
        for pA in unfrozen:
            if abs(pA) == abs(pid):
                continue
            pAlist = [pA]
            if self.M.hasAntiParticle(pA):
                pAlist.append(-pA)
            for pidpair in itertools.product(pAlist,pBlist):
                ppair = tuple(sorted(pidpair))
                if not ppair in self.M.ssmultipliers:
                    # self.record ( f"change ssm of {self.namer.texName(ppair,addDollars=True)} to 1.0" )
                    self.M.ssmultipliers[ppair] = 1.0

    def describe ( self, allTheoryPredictions : bool = False ):
        """ lengthy description of protomodel
        :param allTheoryPredictions: if true, list all theory preds
        """
        sK, sZ = str(self.M.K), str(self.M.Z)
        try:
            sK="%1.2f" % self.M.K
        except:
            pass
        try:
            sZ="%1.2f" % self.M.Z
        except:
            pass
        print( f'\nK = {sK}, Z = {sZ}, muhat = {self.M.muhat:1.2f}, mumax={self.M.mumax:1.3g}' )
        print('  * Best Combo:')
        for tp in self.M.bestCombo:
            txns = ",".join ( set ( map ( str, tp.txnames ) ) )
            eUL = "no ULexp"
            anaId = tp.expResult.globalInfo.id
            dt = tp.dataType( short=True )
            fullId = anaId+":"+dt
            if hasattr ( tp, "expectedUL" ) and type(tp.expectedUL) != type(None):
                eUL = "UL_exp=%1.2f" % tp.expectedUL.asNumber(fb)
            if dt in [ "comb", "combined" ]:
                eUL = f"{tp.getUpperLimit ( expected=True ).asNumber(fb):1.2g}*fb"
            if dt in [ "em", "efficiencyMap" ]:
                dI = tp.dataset.dataInfo
                pred = f"{float ( (tp.xsection.value * tp.dataset.globalInfo.lumi).asNumber() ):1.2g}"
                fullId = anaId+":"+tp.dataset.dataInfo.dataId
                print(f'     - {fullId} [{txns}] obsN={dI.observedN} expBG={dI.expectedBG}+/-{dI.bgError} pred={pred}' )
            else:
                pred=f"{tp.xsection.value.asNumber(fb):.2g}*fb"
                UL = f"{tp.getUpperLimit().asNumber(fb):.2g}*fb"
                print(f'     - {fullId} [{txns}] pred={pred} UL={UL} eUL={eUL}' )

        print ( )
        print('  * Constraints:')
        for tp in sorted( self.M.tpList, key = lambda x: x['robs'], reverse=True ):
            if not allTheoryPredictions and tp['robs'] < 1.0:
                # if not all theory predictions are asked for, only do r>=1
                continue
            txns = ",".join ( set ( map ( str, tp[2].txnames ) ) )
            eUL = ""
            # eUL = ", no ULexp"
            anaId = tp[2].expResult.globalInfo.id
            anaId+=":"+tp[2].dataType()
            dt = tp[2].dataType( short=True )
            if hasattr ( tp[2], "expectedUL" ) and type(tp[2].expectedUL) != type(None):
                eUL = f", UL_exp={tp[2].expectedUL.asNumber(fb):1.2g}*fb"
            if dt in [ "comb", "combined" ]:
                eUL = f", UL_exp={tp[2].getUpperLimit ( expected=True ).asNumber(fb):1.2g}*fb"
            UL=f"{tp[2].upperLimit.asNumber(fb):1.2g}*fb"
            print( f'     - r={tp[2].getRValue():1.2f} {anaId} [{txns}] pred={tp[2].xsection.value.asNumber(fb):1.2g}*fb, UL={UL}{eUL}' )

    def rescaleSignalBy ( self, s : float ):
        """ multiply the signal strength multipliers with s """

        if s == 0.:
            self.pprint ( "rescaling by zero? Ignore." )
            return
        if abs ( s - 1.0 ) < 1e-5:
            return
        self.log ( "rescaling signal by muhat of %.2f" % s )
        self.M.rvalues = [r*s for r in self.M.rvalues[:]]
        self.M.muhat *= 1./s
        self.M.mumax *= 1./s
        self.M.rescaleXSecsBy(s)

        if hasattr(self.M,'tpList'):
            for i,tp in enumerate(self.M.tpList[:]):
                rnew = tp['robs']*s
                if tp['rexp']:
                    rexpnew = tp['rexp']*s
                else:
                    rexpnew = tp['rexp']
                tpNew = tp['tp']
                tpNew.xsection.value *= s #rescale theoryPrediction
                #Remove likelihood and chi2, since they are no longer valid
                #if hasattr(tpNew,'likelihood'):
                #    del tpNew.likelihood
                #if hasattr(tpNew,'chi2'):
                #    del tpNew.chi2
                self.M.tpList[i] = { "robs": rnew,"rexp": rexpnew, "tp": tpNew }
        if hasattr(self.M,'bestCombo'):
            for tp in self.M.bestCombo:
                tp.xsection.value *= s #rescale theoryPrediction
                #Remove likelihood and chi2, since they are no longer valid
                #if hasattr(tp,'likelihood'):
                #    del tp.likelihood
                #if hasattr(tp,'chi2'):
                #    del tp.chi2

    def randomlyChangeModel(self,sigmaUnFreeze : float = 0.5, probBR : float = 0.2,
            probSS : float = 0.25, probSSingle : float = 0.8, ssmSigma : float = 0.1,
            probMerge : float = 0.05, sigmaFreeze : float = 0.5,
            probMassive : float = 0.3, probMass : float = 0.05, dx : float = 200):
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

        nChanges = 0
        nChanges += self.randomlyUnfreezeParticle(sigma=sigmaUnFreeze)
        nChanges += self.randomlyChangeBranchings(prob=probBR)
        nChanges += self.randomlyChangeSignalStrengths(prob = probSS,
                                       probSingle = probSSingle, ssmSigma = ssmSigma)
        nChanges+=self.randomlyFreezeParticle(sigma= sigmaFreeze, probMassive = probMassive)
        if not nChanges: #If nothing has changed, force a random change of masses
            nChanges+=self.randomlyChangeMasses(prob=1.0, dx = dx)
        else: #Change masses with 5% probability
            nChanges+=self.randomlyChangeMasses(prob = probMass, dx = dx)

        #Update cross-sections (if needed)
        self.M.getXsecs()

    def randomlyUnfreezeParticle ( self, sigma=0.5, force = False ) -> int:
        """ Unfreezes a (random) frozen particle according to gaussian distribution
            with a width of <sigma>.

        :param sigma: Width of the gaussian distribution
        :param force: If True force the unfreezing.

        :returns: 1 if a particle got unfrozen, 0 if not.
        """

        #Decide whether to unfreeze according to the number of active particles
        #(always unfreeze if the model only has one particle)
        nUnfrozen = len( self.M.unFrozenParticles() )
        if (not force) and nUnfrozen > 1:
            nTotal = len ( self.M.particles )
            denom = self.M.Z+1.
            if denom < 1.:
                denom = 1.
            mu = 1. - .7 / denom ## make it more unlikely when Z is high
            uUnfreeze = random.gauss( mu ,sigma)
            if uUnfreeze < nUnfrozen/float(nTotal):
                return 0

        self.log ( "unfreeze random particle" )
        # Randomly select the pid:
        frozen = self.M.frozenParticles()
        if len(frozen)==0:
            return 0
        pid = random.choice ( frozen )

        if pid in self.forbiddenparticles:
            self.log ( f"wanted to unfreeze {self.namer.asciiName(pid)} but its forbidden" )
            return 0

        #Check for canonical ordering.
        #If pid matches the heavier state and the lighter state is frozen,
        #unfreeze the lighter state instead
        for pids in self.canonicalOrder:
            if pid == pids[1] and (pids[0] in frozen):
                pid = pids[0] #Unfreeze the lighter state
                break

        self.log ( f"Unfreezing {self.namer.asciiName(pid)}" )
        return self.unFreezeParticle(pid)

    def randomlyChangeBranchings ( self, prob=0.2, zeroBRprob = 0.05, singleBRprob = 0.05 ):
        """ randomly change the branchings of a single particle

        :param prob: Probability for changing a branching ratio
        :param zeroBRprob: With zeroBRprob probability, close decay channel
        :param singleBRprob: With probability singleBRprob, keep only one decay channel

        """

        uBranch = random.uniform(0,1)
        if uBranch < (1-prob):
            return 0

        self.log ( "randomly change branchings" )
        unfrozenparticles = self.M.unFrozenParticles( withLSP=False )
        if len(unfrozenparticles)<2:
            self.pprint ( "not enough unfrozen particles to change random branching" )
            return 0
        p = random.choice ( unfrozenparticles )
        if not p in self.M.decays.keys():
            self.highlight ( "error", "why is %d not in decays?? %s" % ( p, self.M.decays.keys() ) )
            # we dont know about this decay? we initialize with the default!

        return self.randomlyChangeBranchingOfPid ( p, zeroBRprob, singleBRprob )

    def record ( self, change : str ):
        """ log the changes that have been performed on the model
        :param change: a string that describes my latest action
        """
        if not self.do_record:
            return
        self.recording.append ( change )
        if len(self.recording)>20: ## make sure this never explodes
            self.recording = self.recording[-20:]


    def randomlyChangeBranchingOfPid ( self, pid, zeroBRprob = 0.05, singleBRprob = 0.05):
        """ randomly change the branching a particle pid """

        openChannels = self.M.getOpenChannels(pid)
        dkeys = set()
        for dpid in openChannels:
            dk = self.M.decay_keys[pid][dpid]
            dkeys.add(dk)
        
        dkeys = list(dkeys)

        if len(openChannels) < 2:
            self.pprint ( f"number of open channels of {pid} is {len(openChannels)}: cannot change branchings." )
            # not enough channels open to tamper with branchings!
            return 0

        dx = 0.1/numpy.sqrt(len(openChannels)) ## maximum change per channel??

        #Keep only one channel (with probability singleBRprob)
        uSingle = random.uniform( 0., 1. )
        if uSingle < singleBRprob:
            #Choose random decay key:
            dk = random.choice(dkeys)
            decay_chan = [key for key,value in self.M.decay_keys[pid].items() if value == dk]  #get decay channel assocaiated with key, make sure all channels assocaited with same key get same branchings
            self.M.decays[pid] = {}
            br = 1.0/len(decay_chan)
            for dpid in decay_chan:
                self.record ( f"change decay of {self.namer.texName(pid,addDollars=True)} -> {self.namer.texName(dpid,addDollars=True)} to {br:.2f}" )
                self.log ( f"changed decay of {self.namer.asciiName(pid)} -> {self.namer.asciiName(dpid)} to {br:.2f}" )
                self.M.decays[pid].update({dpid: br})
            
            return 1

        #Otherwise randomly change each channel(s) (based on the current BR)
        for dk in dkeys:
            oldbr = 0.
            
            #Check if decay channel already existed:
            decay_chan = [key for key,value in self.M.decay_keys[pid].items() if value == dk]
            if decay_chan[0] in self.M.decays[pid]:
                oldbr = self.M.decays[pid][decay_chan[0]]

            #Close channel(s) (with zeroBRprob probability)
            if oldbr > 0:
                
                uZero = random.uniform( 0., 1. )
                if uZero < zeroBRprob:
                    for dpid in decay_chan:
                        self.record ( f"change branchings of {self.namer.texName(pid,addDollars=True)} -> {self.namer.texName(dpid,addDollars=True)} to 0." )
                        self.log ( f"changed  branchings of {self.namer.asciiName(pid)} -> {self.namer.asciiName(dpid)} from {oldbr:.2f} to 0" )
                        self.M.decays[pid][dpid] = 0.
                    continue

            #Randomly change BR around old value
            Min,Max = max(0.,oldbr-dx), min(oldbr+dx,1.)
            br = random.uniform( Min, Max )/len(decay_chan)
            for dpid in decay_chan:
                self.record ( f"change branchings of {self.namer.texName(pid,addDollars=True)} -> {self.namer.texName(dpid,addDollars=True)} to {br:.2f}" )
                self.log ( f"changed  branchings of {self.namer.asciiName(pid)} -> {self.namer.asciiName(dpid)} to {br:.2f}" )
                self.M.decays[pid][dpid] = br
        
        
        #Make sure there is at least one open channel:
        BRtot = sum(self.M.decays[pid].values())
        if BRtot == 0.0:
            dk = random.choice(dkeys)
            decay_chan = [key for key,value in self.M.decay_keys[pid].items() if value == dk]
            br = 1.0/len(decay_chan)
            self.M.decays[pid] = {}
            for dpid in decay_chan:
                self.record ( f"change decay of {self.namer.texName(pid,addDollars=True)} -> {self.namer.texName(dpid,addDollars=True)} to {br:.2f}" )
                self.log ( f"changed decay of {self.namer.asciiName(pid)} -> {self.namer.asciiName(dpid)} to {br:.2f}" )
                self.M.decays[pid].update({dpid: br})
           
        

        #Make sure BRs add up to 1:
        self.normalizeBranchings(pid)

        return 1

    def randomlyChangeSignalStrengths ( self, prob : float =0.25,
            probSingle : float =0.8, ssmSigma=0.1 ) -> int:
        """ randomly change one of the signal strengths according to a gaussian
        distribution centered around the original SSM.

        :param prob: Probability for changing the signal strengths
        :param probSingle: Probability for changing the signal strength of a
        single particle
        :param ssmSigma: Width for the gaussian

        :returns: 1 if something got changed, else 0
        """

        uSSM = random.uniform(0,1)
        if uSSM < (1-prob):
            return 0

        self.log ( "randomly change signal strengths" )
        if random.uniform(0,1) < probSingle:
            return self.randomlyChangeSSOfOneParticle()
        unfrozenparticles = self.M.unFrozenParticles( withLSP=False )
        if len(unfrozenparticles)<2:
            self.pprint ( "not enough unfrozen particles to change random signal strength" )
            return 0
        #Randomly choose which process pids to change:
        p = random.choice ( unfrozenparticles )
        q = random.choice ( unfrozenparticles )
        #Half of the time select the anti-particle:
        if self.M.hasAntiParticle(p) and random.uniform(0,1)<.5:
            p = -p
        if self.M.hasAntiParticle(q) and random.uniform(0,1)<.5:
            q = -q
        pair = self.M.toTuple(p,q)
        if not pair in self.M.ssmultipliers:
            self.record ( f"change ssm of {self.namer.texName(pair,addDollars=True)} to 1.0" )
            self.M.ssmultipliers[pair]=1.
        newSSM=self.M.ssmultipliers[pair]*random.gauss(1.,ssmSigma) + random.gauss(.1,.1)
        if newSSM < 0.:
            newSSM = 0.
        self.changeSSM(pair,newSSM)
        self.log ( "changing signal strength multiplier of %s,%s: %.2f." % \
                   ( self.namer.asciiName(pair[0]), self.namer.asciiName(pair[1]), newSSM ) )
        self.record ( "change ssm of %s,%s to %.2f." % \
                    ( self.namer.texName(pair[0]), self.namer.texName(pair[1]), newSSM ) )
        return 1

    def randomlyChangeSSOfOneParticle ( self, pid = None ):
        """ randomly change the SS's consistently for one pid
        :param pid: change for this pid. If None, change of a random pid.
        """
        unfrozenparticles = self.M.unFrozenParticles( withLSP=False )

        if len(unfrozenparticles)<2:
            self.pprint ( "not enough unfrozen particles to change random signal strength" )
            return 0
        p = random.choice ( unfrozenparticles )
        if pid != None:
            p = pid
        a = random.uniform ( 0., 1. )
        if a > .9: ## sometimes, just knock out a random SSM
            randomProd = random.choice ( list ( self.M.ssmultipliers.keys() ) )
            self.record ( f"change ssm of {self.namer.texName(randomProd,addDollars=True)} to 1e-5" )
            self.M.ssmultipliers[randomProd]=0.00001
            return 1
        if a < .1: ## sometimes, just try to set to 1.
            randomProd = random.choice ( list ( self.M.ssmultipliers.keys() ) )
            self.record ( f"change ssm of {self.namer.texName(randomProd,addDollars=True)} to 1." )
            self.M.ssmultipliers[randomProd]=1.
            return 1
        if .1 < a < .2: ## sometimes, just try to set to ssm of different particle
            randomProd = random.choice ( list ( self.M.ssmultipliers.keys() ) )
            v = random.choice ( list ( self.M.ssmultipliers.values() ) )
            self.record ( f"change ssm of {self.namer.texName(randomProd,addDollars=True)} to {v:.2f}" )
            self.M.ssmultipliers[randomProd]=v
            return 1
        f = random.uniform ( .8, 1.2 )
        self.log ( f"randomly changing ssms of {self.namer.asciiName ( p )} by a factor of {f:.2f}" )
        ssms = []
        for dpd,v in self.M.ssmultipliers.items():
            if p in dpd or -p in dpd:
                newssm = self.M.ssmultipliers[dpd]*f
                #if newssm > 10000.:
                #    newssm = 10000.
                # self.M.ssmultipliers[dpd]= newssm
                self.changeSSM ( dpd, newssm )
                ssms.append ( newssm )
        self.log ( " `- %s: ssms are now %.2f+/-%.2f" % \
                 ( self.namer.asciiName(p), numpy.mean ( ssms ), numpy.std ( ssms) ) )
        return 1

    def pidPairIsInSSMs ( self, pids : Tuple ) -> bool:
        """ is a given pid pair in SSMs? """
        if pids[1] < pids[0]:
            pids = ( pids[1], pids[0] )
        return pids in self.M.ssmultipliers.keys()

    def changeSSM ( self, pids : Tuple, newssm, recursive : bool = True,
                    verbose : bool = True ):
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
        if newssm > 10000.:
            newssm = 10000.
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


    def randomlyFreezeParticle ( self, sigma= 0.5, probMassive = 0.3):
        """ freezes a random unfrozen particle according to gaussian distribution with width sigma.

        :param sigma: Width of the gaussian distribution
        :param probMassive: Probability for freezing the most massive particle
        """

        nUnfrozen = len( self.M.unFrozenParticles() )
        #Always keep at least 2 particles
        if nUnfrozen <= 2:
            return 0

        nTotal = len ( self.M.particles )
        denom = self.M.Z+1.
        if denom < 1.:
            denom = 1.
        mu = .4 / denom ## make it more unlikely when Z is high
        uFreeze = random.gauss(mu,sigma)
        if uFreeze > nUnfrozen/float(nTotal):
            return 0

        # in every nth step freeze random particle
        if random.uniform(0,1) < probMassive:
            self.log ( "freeze most massive particle" )
            return self.freezeMostMassiveParticle()

        self.log ( "freeze random particle" )
        unfrozen = self.M.unFrozenParticles( withLSP = False )
        if len(unfrozen)<2:
            self.log ( "only two particles are unfrozen, so dont freeze anything" )
            return 0 ## freeze only if at least 3 unfrozen particles exist
        pid = random.choice ( unfrozen )

        self.freezeParticle ( pid )
        return 1

    def freezeMostMassiveParticle ( self, protomodel=None ):
        """ freezes the most massive unfrozen particle """

        if protomodel is None:
            protomodel = self.M

        unfrozen = protomodel.unFrozenParticles( withLSP=False )
        if len(unfrozen)<2:
            return 0 ## freeze only if at least 3 unfrozen particles exist
        pid,minmass=0,0
        for i in unfrozen:
            if protomodel.masses[i]>minmass:
                minmass = protomodel.masses[i]
                pid = i
        # p = random.choice ( unfrozen )
        protomodel.log ( f"Freezing most massive {self.namer.asciiName(pid)} ({minmass:.1f})" )
        self.freezeParticle ( pid, protomodel = protomodel )
        return 1

    def freezeParticle ( self, pid, force = False, protomodel = None ):
        """ freeze particle pid, take care of offshell removal, and
            branching normalization

        :param pid: PID to be frozen
        :param force: If False, will only freeze the particle if it does not violate
                      the canonical order (e.g. will not freeze stop1 if stop2 is unfrozen)
                      and the model contains at least 3 particles.
        :returns: 1 if freezing succesful, 0 if not
        """

        if protomodel is None:
            protomodel = self.M

        #Check for canonical ordering.
        unfrozen = protomodel.unFrozenParticles( withLSP=False )
        if not force:
            if len(unfrozen) < 2:
                return 0
            #If pid matches the lighter state and the heavier state is unfrozen,
            #do not freeze the particle
            for pids in self.canonicalOrder:
                if pid == pids[0] and pids[1] in unfrozen:
                    return 0
        protomodel.log ( f"Freezing {self.namer.asciiName(pid)}" )
        self.record ( f"freeze {self.namer.texName(pid,addDollars=True)}" )
        #Remove pid from masses, decays and signal multipliers:
        if  pid in protomodel.masses:
            protomodel.masses.pop(pid)
        if  pid in protomodel.decays:
            protomodel.decays.pop(pid)
        removeSSM = [pids for pids in protomodel.ssmultipliers if (pid in pids or -pid in pids)]
        for pids in removeSSM:
            protomodel.ssmultipliers.pop(pids)

        #Fix branching ratios and rescale signal strenghts, so other channels are not affected
        self.removeAllOffshell(rescaleSSMs=True, protomodel=protomodel)
        return 1

    def unFreezeParticle (self, pid, force = False, protomodel = None):
        """ unfreeze particle pid, assign masses, BRs and signal strength multipliers.

        :param pid: PID to be unfrozen
        :param force: If False, will only unfreeze the particle if it does not violate
                      the canonical order (e.g. will not unfreeze stop2 if stop1 is frozen).
        """

        if protomodel is None:
            protomodel = self.M

        #Check for canonical ordering.
        frozen = protomodel.frozenParticles( )
        if not force:
            #If pid matches the heavier state and the lighter state is frozen,
            #do not unfreeze the particle
            for pids in self.canonicalOrder:
                if pid == pids[1] and pids[0] in frozen:
                    return 0

        #Absolute mass range:
        maxMass = protomodel.maxMass
        minMass = protomodel.masses[protomodel.LSP]

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

        offshell = False
        if pid in [ 1000023, 1000024 ]:
            # for C1 and N2 we want a 10% chance to start in the offshell regime
            p = random.uniform ( 0, 1 )
            if p < 0.1:
                offshell = True
                self.log ( f"Unfreezing {self.namer.asciiName(pid)}, randomly chose to restrict to offshell mass!" )
                maxMax = minMass + 90.

        tmpMass = random.uniform ( minMass, maxMass )
        ctr = 0
        while pid in [ 1000006, 2000006 ] and \
                self.inCorridorRegion ( tmpMass, protomodel.masses[protomodel.LSP] ):
            # if in corridor region, redraw!
            tmpMass = random.uniform ( minMass, maxMass )
            ctr += 1
            if ctr > 5: ## seems like the air is too thin. make more space.
                mstop2 = 2000.
                if 2000006 in protomodel.masses:
                    mstop2 = protomodel.masses[2000006]
                protomodel.masses[2000006] = mstop2 + 20.
                if pid == 1000006:
                    maxMass = mstop2 + 20.

        self.record ( f"unfreeze {self.namer.texName(pid,addDollars=True)} to "\
                      f"{tmpMass:.1f}" )
        #Randomly select mass of unfrozen particle:
        protomodel.masses[pid] = tmpMass
        self.pprint ( f"Unfroze mass of {self.namer.asciiName(pid)} to "\
                      f"{protomodel.masses[pid]:.1f}" )

      
        # Set branchings
        self.initBranchings(pid)

        #Add pid pair production and associated production to protomodel.ssmultipliers:
        self.initSSMFor(pid)

        return 1

    def randomlyChangeMasses ( self, prob = 0.05, dx = 200.0 ):
        """ take a random step in mass space for a single unfrozen particle

        :param prob: Probability for changing the mass
        :param dx: Defines the interval for selecting the delta m (-dx,dx) """

        uMass = random.uniform ( 0., 1. )
        if uMass < (1-prob):
            return 0

        self.log ( "take random mass step" )
        unfrozen = self.M.unFrozenParticles()
        if len(unfrozen)==0:
            return 0

        pid = random.choice ( unfrozen )

        #Define mass interval
        maxMass = self.M.maxMass
        minMass = self.M.masses[self.M.LSP]
        #In case the pid corresponds to a lighter or heavier state of a pair of particles,
        #make sure the mass ordering is respected:
        for pids in self.canonicalOrder:
            if pid == pids[0] and pids[1] in unfrozen:
                # Do not allow for masses above the heavier state
                maxMass = self.M.masses[pids[1]]
            elif pid == pids[1] and pids[0] in unfrozen:
                # Do not allow for masses below the ligher state
                minMass = self.M.masses[pids[0]]

        #If the particle is the LSP, make sure its mass remains the lightest, but relax the lower limit
        if pid == self.M.LSP:
            minMass = 10.0
            if len(unfrozen) > 1:
                maxMass = min([self.M.masses[p] for p in unfrozen if p != self.M.LSP])
        # an artificial wall because the maps are bounded from below
        if abs(pid) in self.walledpids:
            ## heed the wall!
            minMass = max ( self.walledpids[abs(pid)], minMass )

        ret = self.randomlyChangeMassOf ( pid, dx=dx, minMass=minMass, maxMass=maxMass )
        if pid in [ 1000023, 1000024 ]:
            # for C1 and N2, if one of the two gets changed, have a 10% that
            # the other gets set to the same value
            p=random.uniform(0,1)
            if p < .1:
                mass = self.M.masses[pid]
                otherpid = 1000024 if pid == 1000023 else 1000023
                ## remember the frozen particles, so we can check
                # if we just unfroze this guy
                were_frozen = self.M.frozenParticles()
                self.M.masses[otherpid] = mass * random.uniform ( .99, 1.01 )
                self.log ( f"mass of {self.namer.asciiName(pid)} got changed to {mass:.1f}. hattrick, changing also for {self.namer.asciiName(otherpid)}!" )
                # FIXME if the particle was frozen before, we need to
                # unfreeze
                if otherpid in were_frozen:
                    self.initBranchings(otherpid)
                    self.initSSMFor(otherpid)
                self.record ( f"change mass of {self.namer.asciiName(otherpid)} to {self.M.masses[otherpid]}" )

        #Fix branching ratios and rescale signal strenghts, so other channels are not affected
        self.removeAllOffshell(rescaleSSMs=True)

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

        :returns: 1 for success
        """
        if dx == None:
            denom = self.M.Z + 1.
            if denom < 1. or denom == None:
                denom = 1.
            dx = 40. / numpy.sqrt ( len(self.M.unFrozenParticles() ) ) / denom
        if dx in [ float("nan"), float("inf"), None ] or dx > 200.:
            dx=40.

        if not minMass:
            minMass = self.M.masses[self.M.LSP]
        if not maxMass:
            maxMass = self.M.maxMass
        massIsLegal = False
        ctIterations = 0
        while not massIsLegal:
            ctIterations += 1
            massIsLegal = True
            tmpmass = self.M.masses[pid]+random.uniform(-dx,dx)
            # Enforce mass interval:
            if pid in [ 1000006, 2000006 ] and \
                    self.inCorridorRegion ( tmpmass, self.M.masses[self.M.LSP] ):
                massIsLegal = False
            if pid == self.M.LSP and 1000006 in self.M.masses and \
                    self.inCorridorRegion ( self.M.masses[1000006], tmpmass ):
                massIsLegal = False
            if pid == self.M.LSP and 2000006 in self.M.masses and \
                    self.inCorridorRegion ( self.M.masses[2000006], tmpmass ):
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
        self.pprint ( f"randomly changing mass of {self.namer.asciiName ( pid )} to {tmpmass:.1f}" )
        self.record ( f"change mass of {self.namer.texName(pid,addDollars=True)} to {tmpmass:.1f}" )
        self.M.masses[pid]=tmpmass

        return 1

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

    def mergeParticles(self,dm=200.,protomodel=None):
        """ Look for pair of candidates with mass difference smaller than dm and merge them.
            If several particles can be merged, only merge the ones with the smallest mass difference.
            If protomodel is defined merge the particles of the given model, else merge particles in self.M

            :param dm: Maximum mass difference for merging
            :param protomodel: ProtoModel to be modified. If None, use self.M

            :return: False if no merge was performed, else returns True
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
                self.pprint("can not merge particles with wrong mass hierarchy (%d > %d)" %(pidA,pidB))
                return False
            if pidA in [ 1000006, 2000006 ] and pidB in [ 1000006, 2000006 ]:
                ## merging stops. check if we would end up in corridor.
                avgM = self.computeAvgMass ( (pidA,pidB) )
                if self.inCorridorRegion ( avgM, self.M.masses[self.M.LSP] ):
                    self.pprint ( "wont merge the stops since we would end up in corridor region!" )
                    return False
            self.merge((pidA,pidB),protomodel)
            return True
        else:
            return False

    def merge ( self, pair, protomodel = None):
        """ merge the particles with pidA and pidB in protomodel.

        :param pair: Pair of particle pids to be merged
        :param protomodel: ProtoModel to be modified. If None, use self.M

        :return: Protomodel with the particles merged
        """

        if not protomodel:
            protomodel = self.M

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
        self.pprint ( "merging %s and %s" % \
                ( self.namer.asciiName(p1), self.namer.asciiName( p2 ) ) )
        self.log ( "masses before merger: %.2f, %.2f" % \
                   ( protomodel.masses[p1], protomodel.masses[p2] ) )
        avgM = self.computeAvgMass ( pair )
        self.log ( "avg mass for %s is %.1f" % ( str(pair), avgM ) )
        protomodel.masses[ p1 ] = avgM ## set this one to the avg mass

        #Get p2 decays:
        p2decays = protomodel.decays[p2]
        #Get allowed decay channels for p1:
        openChannels = self.M.getOpenChannels(p1)
        ## add the decays from pid2 to pid1 if decay is allowed:
        for pids,br in p2decays.items():
            if not pids in openChannels:
                continue
            if pids in protomodel.decays[p1]:
                if br > 0.001:
                    self.log ( "add to decays %s/%s: %.2f" % ( p1, pids, br ) )
                protomodel.decays[p1][pids] += br
            else:
                self.log ( "set decays of %s/%s to %.2f" % ( p1, pids, br ) )
                protomodel.decays[p1][pids] = br

        self.log ( f"now normalize branchings of {self.namer.asciiName(p1)}" )
        self.normalizeBranchings ( p1, protomodel=protomodel )

        #Now replace all decays to p2 by decays to p1
        #(since the new (average) p1 mass is always smaller than the p2 mass,
        #there is no chance of running into offshell decays)
        for mpid,decays in olddecays.items():
            for dpids,br in decays.items():
                newpids = dpids
                #Replace any appearence of p2/-p2 in decays by p1/-p1:
                if isinstance(dpids,(list, tuple)):
                    newpids = [dpid if abs(dpid) != abs(p2) else p1*dpid/abs(dpid)
                                for dpid in dpids]
                    newpids = tuple(newpids)
                elif isinstance(dpids,int) and abs(dpids) == abs(p2):
                    newpids = p1*dpids/abs(dpids)

                #If original channel did not contain p2, do nothing
                if newpids == dpids:
                    continue
                #Print log message for non-negligible BRs:s
                if br > 0.0001:
                    self.log ( "redirecting decay of %d from %s to %s: br=%.2f" % \
                               ( mpid, dpids, newpids, br ) )

                #If new channel was already present, simply add to BR:
                if newpids in protomodel.decays[mpid]:
                    br += protomodel.decays[mpid][newpids]

                #Remove original decay:
                protomodel.decays[mpid].pop ( dpids )
                #Add new channel:
                protomodel.decays[mpid][newpids]=br

        if oldxsecs != None:
            ## merge the signal strength multipliers:
            self.mergeSSMs( pair, oldXsecs = oldxsecs, protomodel=protomodel )

        ## finally freeze p2:
        self.freezeParticle(p2,protomodel=protomodel)

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

    def mergeSSMs ( self, pair, oldXsecs, protomodel=None ):
        """ merge signal strength multipliers for particles in pair. The cross-selections
            involving the merged particles are assumed to be added and the corresponding
            signal strengths are rescaled.

        :param pair: pair of particle PIDs being merged
        :param oldXsecs: cross-sections before the merge
        :param protomodel: protomodel to be modified. If not defined, use self.M
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
        procDict = {}
        for pids,xsec in oldxsecDict.items():
            if not p2 in pids and not -p2 in pids:
                continue
            newpids = [pid if abs(pid) != abs(p2) else p1*p2/abs(p2) for pid in pids ]
            newpids = tuple(sorted(newpids))
            #Skip processes containing frozen particles:
            if not all([abs(pid) in unfrozen for pid in newpids]):
                continue
            p2Xsecs[pids] = xsec
            procDict[pids] = newpids

        #Now compute the new SSMs assuming that the cross-sections will be added:
        newSSMs= {}
        for oldpid,newpid in procDict.items():
           #Get xsec value for the process containing p2:
           value = p2Xsecs[oldpid]
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
               newSSMs[newpids] = oldssm*(oldvalue/value)
           #If the new process does not exist take the SSM for the (old) process containing p2
           else:
               oldssm = 1.0
               if oldpid in protomodel.ssmultipliers:
                   oldssm = protomodel.ssmultipliers[oldpid]
               newSSMs[newpids] = oldssm

        #Now replace the SSMs in protomodel:
        for pid,ssm in newSSMs.items():
           protomodel.ssmultipliers[pid] = ssm

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

    def printXSecs ( self, fbmin=.001*fb ):
        """ print the cross sections in a human-readable way """
        xsecs = self.simplifyXSecs( fbmin )
        for sqrts in xsecs.keys():
            pidss = list ( xsecs[sqrts].keys() ) # list of list of pids
            pidss.sort()
            print ( f"{sqrts} TeV:" )
            for pids in pidss:
                xsec = xsecs[sqrts][pids]
                label = ""
                if "dict" in xsec.info.label:
                    label = " (from dict)"
                print ( f" {str(pids):>22s}: {xsec.value.asNumber(fb):.2f} fb{label}" )

    def simplifyDecays ( self ):
        """ return the decays only of the unfrozen particles,
            only != 0 """
        ret ={}
        unfrozen = self.M.unFrozenParticles()
        for mpid,decays in self.M.decays.items():
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

    def xsecsFor ( self, pids, sqrts=13*TeV, order=LO ):
        """ return the cross sections for pids.
        :param pids: tuple of two pids
        :returns: cross section (that had the SSM applied),
                  and SSM that *was* applied.
        """

        ssm = 1.
        if pids[1] < pids[0]:
            pids = ( pids[1], pids[0] )
        if pids in self.M.ssmultipliers:
            ssm = self.M.ssmultipliers[pids]
        xs = 0. * fb
        modelXSecs = self.M.getXsecs()[0]
        for xsec in modelXSecs:
            if xsec.info.order != order:
                continue
            if abs ( ( xsec.info.sqrts - sqrts ).asNumber(TeV) ) > .1:
                continue
            if xsec.pid != pids:
                continue
            xs = xsec.value
        return xs,ssm

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

    def getAllPidsOfBestCombo ( self ):
        """ get all pids that appear in the best combo """
        ret = set()
        if self.M.bestCombo is None:
            return ret
        for tp in self.M.bestCombo:
            for prod in tp.PIDs:
                for branch in prod:
                    for pid in branch:
                        if type(pid) == int:
                            ret.add ( abs(pid) )
                        if type(pid) in [ tuple, list ] and len(pid)>0:
                            ret.add ( abs(pid[0]) )

        return ret

    def freezePidsNotInBestCombo ( self ):
        """ all pids that arent in best combo but have
            unfrozen masses -- freeze them """
        okPids = self.getAllPidsOfBestCombo()
        if len(okPids)==0: ## means we dont have a best combo
            return 0
        unfrozen = self.M.unFrozenParticles( withLSP=False )
        nfrozen = 0
        for pid in unfrozen:
            if not pid in okPids:
                nfrozen += self.freezeParticle ( pid )
        return nfrozen

    def backupModel ( self ):
        """ backup the current state """

        self._backup = { "llhd": self.M.llhd, "letters": self.M.letters, "Z": self.M.Z,
                         "K": self.M.K, "muhat": self.M.muhat,
                         "description": self.M.description,
                         "tpList": copy.deepcopy(self.M.tpList),
                         "bestCombo": copy.deepcopy(self.M.bestCombo),
                         "masses": copy.deepcopy(self.M.masses),
                         "ssmultipliers": copy.deepcopy(self.M.ssmultipliers),
                         "decays": copy.deepcopy(self.M.decays),
                         "rvalues": copy.deepcopy(self.M.rvalues),
                         "_stored_xsecs" : copy.deepcopy(self.M._stored_xsecs),
                         "_xsecMasses" : copy.deepcopy(self.M._xsecMasses),
                         "_xsecSSMs" : copy.deepcopy(self.M._xsecSSMs),
                         }

    def restoreModel ( self, reportReversion=False ):
        """ restore from the backup """
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
