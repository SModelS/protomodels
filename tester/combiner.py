#!/usr/bin/env python3

""" Try out combinations from pickle file. """

__all__ = [ "Combiner", "loadUserPrior" ]

from smodels.decomposition import decomposer
from smodels.share.models.SMparticles import SMList
from pshare.model_spec import BSMList
from smodels.base.physicsUnits import fb
from smodels.statistics.basicStats import NllEvalType, observed, apriori
from smodels.matching.theoryPrediction import TheoryPrediction, \
         TheoryPredictionsCombiner
import sys, os, inspect, importlib.util
sys.path.insert(0,os.path.abspath ( os.path.dirname(__file__) ) )
from pbase.loggerbase import LoggerBase
import numpy, math, copy, sys, time
from colorama import Fore
from typing import List, Union, Tuple, Set, Callable
from ptools.helpers import getAllPidsOfTheoryPred, getAllProdModesOfTheoryPred
from ptools.bamCreator import selectMostSignificantSRs, bamAndWeights, \
         find_best_comb

def loadUserPrior ( path : str ) -> Callable:
    """ load the (optional) user-provided prior from a python file

    The file has to contain a function 'myprior' which takes the protomodel
    and the nll flag (as options, i.e. keyword arguments 'protomodel' and
    'nll') and returns a float. That float is then added to the standard
    prior (for nll=True) or multiplied to it (for nll=False).

    We do a few basic checks on the file and on the function:
    - the file has to exist
    - the file has to be importable (no exceptions while loading it)
    - it has to contain a callable 'myprior'
    - 'myprior' has to accept the protomodel and the nll flag

    :param path: path to the python file containing the 'myprior' function
    :returns: a function that takes (protomodel, nll) and returns a float
    :raises: Exception, if any of the checks above fails
    """
    if path is None:
        raise ValueError ( "no path to a user prior given" )
    if isinstance ( path, os.PathLike ):
        path = os.fspath ( path )
    if not isinstance ( path, str ):
        raise ValueError ( f"the path to the user prior has to be a string, but is {type(path)}" )
    fname = os.path.abspath ( os.path.expanduser ( path ) )
    if not os.path.exists ( fname ):
        raise FileNotFoundError ( f"the user prior file {fname} does not exist" )
    if not os.path.isfile ( fname ):
        raise FileNotFoundError ( f"the user prior file {fname} is not a file" )
    ## make the directory of the user file importable, so that the user can
    ## import helper modules from there
    dirname = os.path.dirname ( fname )
    if dirname not in sys.path:
        sys.path.insert ( 0, dirname )
    modname = "userprior_" + os.path.basename ( fname ).replace ( ".", "_" )
    try:
        spec = importlib.util.spec_from_file_location ( modname, fname )
    except Exception as e:
        raise Exception ( f"could not load the user prior file {fname}: {e}" )
    if spec is None or spec.loader is None:
        raise Exception ( f"could not load the user prior file {fname}" )
    module = importlib.util.module_from_spec ( spec )
    sys.modules [ modname ] = module
    try:
        spec.loader.exec_module ( module )
    except Exception as e:
        raise Exception ( f"error while executing the user prior file {fname}: {e}" )
    if not hasattr ( module, "myprior" ):
        raise Exception ( f"the user prior file {fname} does not contain a 'myprior' function" )
    func = getattr ( module, "myprior" )
    if not callable ( func ):
        raise Exception ( f"the 'myprior' in {fname} is not callable" )
    ## now check that we can call it with the protomodel and the nll flag
    try:
        sig = inspect.signature ( func )
    except ( TypeError, ValueError ) as e:
        raise Exception ( f"could not determine the signature of 'myprior' in {fname}: {e}" )
    params = list ( sig.parameters.values() )
    has_var_positional = any ( p.kind == p.VAR_POSITIONAL for p in params )
    has_var_keyword = any ( p.kind == p.VAR_KEYWORD for p in params )
    npositional = len ( [ p for p in params if p.kind in
                        ( p.POSITIONAL_ONLY, p.POSITIONAL_OR_KEYWORD ) ] )
    keyword_names = [ p.name for p in params if p.kind in
                     ( p.POSITIONAL_OR_KEYWORD, p.KEYWORD_ONLY ) ]
    use_keywords = ( ( "protomodel" in keyword_names or has_var_keyword ) and
                     ( "nll" in keyword_names or has_var_keyword ) )
    use_positional = npositional >= 2 or has_var_positional
    if not use_keywords and not use_positional:
        raise Exception ( f"the 'myprior' function in {fname} has to accept the protomodel "\
                          f"and the nll flag, but its signature is {sig}" )
    def userprior ( protomodel, nll : bool ) -> float:
        """ the wrapper that calls the user function """
        if use_keywords:
            return func ( protomodel=protomodel, nll=nll )
        return func ( protomodel, nll )
    userprior.__doc__ = f"user prior from {fname}"
    return userprior

class Combiner ( LoggerBase ):
    def __init__ ( self, walkerid: int=0,
                   userprior : Union[None,str,Callable] = None ):
        """ initialise the combiner

        :param walkerid: the id of the walker (for logging)
        :param userprior: optional path to a python file that contains a
        'myprior' function (taking the protomodel and the nll flag), which is
        then used to modify the standard prior. See loadUserPrior.
        """
        super ( Combiner, self ).__init__ ( walkerid )
        self.walkerid = walkerid
        self.userprior = None
        self.userpriorFile = None
        if userprior not in [ None, "", False ]:
            self.setUserPrior ( userprior )

    def setUserPrior ( self, userprior : Union[str,Callable] ) -> None:
        """ set the (optional) user-provided prior

        :param userprior: path to a python file containing a 'myprior'
        function, or the function itself.
        """
        if callable ( userprior ):
            self.userprior = userprior
            self.userpriorFile = None
            return
        self.userpriorFile = os.path.abspath ( os.path.expanduser ( userprior ) )
        self.userprior = loadUserPrior ( self.userpriorFile )
        self.pprint ( f"using the user prior from {self.userpriorFile}" )

    def getAllPidsOfCombo ( self, combo : List[TheoryPrediction] ) -> Set:
        """ get all particle IDs (PIDs) that make it into one combo """
        pids = set()
        for theoryPred in combo:
            pids = pids.union ( getAllPidsOfTheoryPred ( theoryPred ) )
        return pids

    def getAllSSMsOfCombo ( self, combo : List[TheoryPrediction] ) -> Set:
        """ get all production modes (PID1, PID2) that make it into one combo """
        prod_modes = set()
        for theoryPred in combo:
            prod_modes = prod_modes.union ( getAllProdModesOfTheoryPred ( theoryPred ) )
        return prod_modes

    def getAnaIdsWithPids ( self, combo, pids ):
        # from best combo, retrieve all ana ids that contain *all* pids
        anaIds = set()
        for theoryPred in combo:
            tpids = getAllPidsOfTheoryPred ( theoryPred )
            hasAllPids=True
            for pid in pids:
                if not pid in tpids:
                    hasAllPids=False
                    break
            if hasAllPids:
                anaIds.add ( theoryPred.analysisId() )
        return anaIds

    def printLLhds ( self, llhds ):
        keys = list ( llhds.keys() )
        keys.sort()
        for k in keys:
            v=llhds[k]
            self.pprint ( f"{k:.2f}: {v:.3g}" )

    def getLetters( self, predictions ):
        '''assign a letter to every prediction. for debugging'''

        letters={}
        if predictions is None:
            return letters

        letter=65                   #char = A
        # self.pprint ( "[combine] Letters assigned to results:" )
        for p in predictions:
            letters[p]=chr(letter)
            # self.pprint ( "[combine] Prediction %s: %s" % ( letters[p], p.expResult.globalInfo.id ) )
            letter+=1
            if letter == 91:        # skip the special characters
                letter = 97         #char = a
        return letters

    def getComboDescription ( self, combination ):
        def describe ( x ):
            return f"{x.analysisId()}({x.dataType().replace('upperLimit', 'ul').replace('efficiencyMap', 'em').replace('combined', 'comb')})"
        return ",".join( [ describe(x) for x in combination ] )

    def getLetterCode ( self, combination ):
        """ get the letter code of a combination """
        ret = ""
        for c in combination:
            ret += self.letters[c]
        return ret

    def priorForNDF ( self, nparticles : int, nbranchings : int, nssms : int, 
                    name : str ="expo1", verbose : bool = False, 
                      nll : bool =False ) -> float:
        """ Get the prior for this and this many degrees of freedom in the model.
        :param nparticles: number of unfrozen particles
        :param nbranchings: number of branchings > 0 and < 1
        :param nssms: number of signal strength multipliers > 0
        :param name: name of the prior, cause I will be defining quite a few.
        :param verbose: be verbose about computation
        :param nll: if true, compute nll of prior
        :returns: *proper* prior
        """
        if name == "flat":
            prior = 1.
        if name == "expo1":
            a,b,c = 2, 4, 8
            prior = numpy.exp ( -1 * ( nparticles/a + nbranchings/b + nssms/c ) )
        if name == "expo2":
            a,b,c = 3, 3.68*3, 5.7*3
            prior = numpy.exp ( -1 * ( nparticles/a + nbranchings/b + nssms/c ) )
        if name == "gauss1":
            a,b,c = 2, 8, 32 ## the "sigmas" of the Gaussians. Higher values means less punishment
            prior = numpy.exp ( -(1/2) * ( (nparticles/a)**2 + (nbranchings/b)**2 + (nssms/c)**2 ) )
        # verbose = True
        if verbose:
            self.log ( f"prior ``{name}'': {nparticles} particles, {nbranchings} branchings, {nssms:.1f} equivalent unique ssms: {prior}" )
        if nll:
            return - numpy.log ( prior )
        return prior

    def noSuchBranching ( self, branchings, br ):
        """ check if a branching ratio similar to br already exists in branchings """
        for cbr in branchings:
            if abs ( cbr - br ) / ( cbr + br ) < 0.01 / 2.: ## 1 percent rule
                return False
        return True

    def computePrior ( self, protomodel, nll : bool = False, 
                  verbose : bool =False, name : str ="expo1" ) -> float:
        """ compute the prior for protomodel, used to introduce regularization,
            i.e. penalizing for non-zero parameters, imposing sparsity.

        If a user prior has been provided (see self.setUserPrior), it is
        applied at the very end: the value returned by the user function is
        added to the prior (for nll=True) or multiplied to it (for nll=False).

        :param nll: if True, return negative log likelihood
        :param verbose: print how you get the prior.
        :param name: name of prior (expo1, gauss1, etc). See self.priorForNDF.

        :returns: a priori probability of the model
        """

        unfrozen = protomodel.unFrozenParticles ( withLSP=True )
        nUnfrozen = len ( unfrozen )
        nbr = 0
        ## every non-trivial branching costs something
        for mpid,decays in protomodel.decays.items():
            if not mpid in unfrozen or mpid == protomodel.LSP:
                continue ## frozen particles dont count
            ## memorize branchings, similar branchings count only once
            memBRs = set() 
            for dpid,br in decays.items():
                if br > 1e-5 and self.noSuchBranching ( memBRs, br ):
                    memBRs.add ( br )
            tmp = len ( memBRs ) - 1 ## subtract one cause they add up to 1.
            nbr += tmp

        ## every non-trivial signal strength multiplier costs something
        cssms = {}
        pun1 = 0. # punishment for ssm=1, we prefer zeroes!
        xsecs = protomodel.getXsecs()
        if len(xsecs)==0:
            self.pprint( f"could not get xsecs for protomodel with masses {protomodel.masses!s}" )
            if nll:
                return self.applyUserPrior ( -float("inf"), protomodel, nll )
            return self.applyUserPrior ( 0., protomodel, nll )
        modelXSecs = xsecs[0]
        for pids,ssm in protomodel.ssmultipliers.items():
            if (abs(pids[0]) not in unfrozen) or (abs(pids[1]) not in unfrozen):
                continue
            ## ignore non-existant xsecs
            hasXSec=False
            xsecv = 0.*fb
            for xsec in modelXSecs:
                if pids == xsec.pid:
                    hasXSec = True
                    xsecv = xsec.value
            ## ignore all very small xsecs
            if not hasXSec or xsecv < .001 * fb:
                continue
            ## every unique ssm > 0 costs a little, but only very little
            # ssmkey = int ( 100. * ssm )
            if ssm > 1e-4: # and abs ( ssm - 1. ) > .01: ssms of 1, are they special?
                foundSimilarValue=False
                for k,v in cssms.items():
                    dssm = 2. * abs ( ssm - k ) / ( ssm + k )
                    if dssm < .03: ## they are the same up to 3% percent?
                        cssms[k]+=1
                        foundSimilarValue=True
                        break
                if not foundSimilarValue:
                    cssms[ssm]=1
                """
                if not ssmkey in cssms:
                    cssms[ssmkey]=0
                cssms[ssmkey]+=1
                """
                # nssms += 1
            ## if we treat ones a special, see above, then we need this
            #if abs ( ssm - 1. ) < .01:
            #    pun1 += .1 ## we prefer zeroes over ones, so we punish the ones
        # print ( "cssms", cssms )
        pun1 += .1 * ( sum(cssms.values()) - len(cssms) ) ## small additional punishments for all non-zeros
        nssms = len( cssms )+pun1
        ret = self.priorForNDF ( nUnfrozen, nbr, nssms, name, verbose )
        if verbose:
            ssmstring = [ f"{x:.2f}" for x in cssms.keys() ]
            self.log ( f"           `- the unique ssms are: {', '.join(ssmstring)}" )
        if nll:
            ret = - math.log ( ret )

        ## susy_mode!
        if protomodel.environ.susy_mode == True:
            nll_penalty = 0.
            for pids, ssm in protomodel.ssmultipliers.items():
                # dnll = .1*float(abs(ssm-1.)) L1 penalty
                dnll = .1*float((ssm-1.)**2) # L2 penalty
                nll_penalty += dnll
            if nll:
                ret += nll_penalty
            else:
                ret *= numpy.exp ( - nll_penalty )
        pen_widths = self.penaltyForSillyWidths ( protomodel )
        if nll:
            ret += pen_widths
        else:
            ret *= numpy.exp ( pen_widths )
        ## and finally the (optional) user-provided prior
        #print ( f"@@2 computePrior {ret} nll {nll}" )
        ret = self.applyUserPrior ( ret, protomodel, nll )
        # print ( f"@@3 computePrior {ret}" )
        return float(ret)

    def applyUserPrior ( self, prior : float, protomodel,
                          nll : bool ) -> float:
        """ apply the (optional) user-provided prior to the standard prior

        :param prior: the standard prior (or nll, if nll=True)
        :param protomodel: the model the prior is computed for
        :param nll: if True, the user prior is added, else it is multiplied
        :returns: the modified prior
        """
        if self.userprior is None:
            return prior
        try:
            val = self.userprior ( protomodel, nll )
        except Exception as e:
            self.highlight ( "error", f"error in the user prior: {e}" )
            raise
        ## basic checks on the value returned by the user function
        if isinstance ( val, bool ) or not isinstance ( val, ( int, float,
                                                numpy.floating,
                                                numpy.integer ) ):
            raise TypeError ( f"the user prior has to return a float, "\
                              f"but it returned {type(val)}" )
        val = float ( val )
        if math.isnan ( val ) or math.isinf ( val ):
            raise ValueError ( f"the user prior has to return a finite float, "\
                               f"but it returned {val}" )
        if not nll and val < 0.:
            raise ValueError ( f"the user prior must not be negative when it is "\
                               f"multiplied to the prior, but it returned {val}" )
        if nll:
            return prior + val
        return prior * val

    def penaltyForSillyWidths ( self, protomodel )-> float:
        """ if widths dont fit with the mass gaps, penalize for this here
        :returns: penalty, as nll!!
        """
        if not hasattr ( protomodel, "widths" ):
            return 0.
        if protomodel.widths == None:
            return 0.
        if len ( protomodel.widths ) == 0:
            return 0.
        return 0.

    def penaltyForMissingResults ( self, predictions : List[TheoryPrediction] )\
                -> float:
        """ very simple hack for now, penalize if predictions are all from 
            the same experiment
        :returns: penalty -- 1e-3 if experiment is missing
        """
        if len(predictions)==0:
            return .01 # penalize but doesnt matter, there is no likelihood
        hasExperiment = { "ATLAS": False, "CMS": False }
        numExperiment = { "ATLAS": 0, "CMS": 0 }
        for p in predictions:
            for experiment in [ "CMS", "ATLAS" ]:
                if experiment in p.dataset.globalInfo.id:
                    hasExperiment[experiment]=True
                    numExperiment[experiment] += 1
        if hasExperiment["ATLAS"] and hasExperiment["CMS"]:
            return 1.
        elif numExperiment["ATLAS"] > 2 or numExperiment["CMS"] > 2:
            num = numExperiment["ATLAS"] if numExperiment["ATLAS"] > 0 else numExperiment["CMS"]
            self.log ( f"penalty! we only have {num} {','.join( k for k,v in hasExperiment.items() if v )} results" )
            return 0.8
        else:
            num = numExperiment["ATLAS"] if numExperiment["ATLAS"] > 0 else numExperiment["CMS"]
            self.log ( f"penalty! we only have {num} {','.join( k for k,v in hasExperiment.items() if v )} results" )
            return 0.5 # penalize!

    def hasNBodyDecays ( self, decays : dict , n : int = 3 ) -> bool:
        """ among the decays are there also n-body decays?
        :returns: true if pid decays via n-body decays
        """
        for decay,br in decays.items():
            if len(decay)==n:
                return True
        return False

    def hasHeavyBosonDecays ( self, decays : dict ) -> bool:
        """ amonst decays, are there also decays via heavy gauge bosons?
        :returns: true if pid decays via W or Z or higgs
        """
        for decay,br in decays.items():
            for g in [ 23, 24, 25 ]:
                if g in decay:
                    return True
        return False

    def penaltyForUndemocraticFlavors ( self, protomodel ) -> float:
        """ very simple hack for now, penalize for undemocratic flavor decays

        :returns: penalty -- 1 - delta(branching_ratios)/C for each flavor decay
        """
        ret = 1.
        for pid, decays in protomodel.decays.items():
            if not pid in [ 1000024, 1000023 ]:
                continue
            if len (decays) == 0:
                continue

            if all( [ len(decay) == 2 for decay in decays.keys() ] ):  # On-shell, continue
                continue
            if self.hasHeavyBosonDecays( decays ) and self.hasNBodyDecays ( decays ):
                #from termcolor import colored
                self.highlight("error", f"ERROR: pID {pid} of mass {protomodel.masses[pid]} has both on- AND off-shell decays: {decays}! Ignoring for now")
# self.highlight("error", f"ERROR: pID {pid} of mass {protomodel.masses[pid]} has both on- AND off-shell decays: {decays}! Returning a prior of 0.05")
                #import sys, IPython; IPython.embed( colors = "neutral" ); sys.exit()
                # from termcolor import colored
                # return 0.05
                continue
            if len ( decays ) == 1:
                if list(decays.keys())[0] in [(1000022, 2, 1), (1000022, 2, 2), (1000022, 5, 5)]: # If the open channel is for light quarks or bb, allow
                    continue

            lep_decays = 0
            lep_br = []
            C1_has_lep = False
            for decay in decays.keys():
                if pid == 1000024 and decay in [(1000022, 11, 12), (1000022, 13, 14), (1000022, 15, 16)] and decays[decay] != 0.:
                    lep_decays += 1
                    lep_br.append(decays[decay])
                    continue
                if pid == 1000023 and decay in [(1000022, 11, 11), (1000022, 13, 13), (1000022, 15, 15)]:
                    if decays[decay] == 0.:
                        C1_has_lep = True # Chargino 1 has leptonic decays but set to 0
                    else:
                        lep_decays += 1
                        lep_br.append(decays[decay])
                    continue
            if lep_decays == 0: # Multiple open channels but no leptonic one, ok for neutralino 2, disallow for chargino 1
                if pid == 1000023:
                    continue
                elif pid == 1000024 and not C1_has_lep:
                    self.pprint(f"Chargino 1 has multiple open channels but no leptonic one: {decays}! Returning a prior of 0.05")
                    return 0.05

            # Penalise for each missing leptonic channel
            mtau = 1.78 # GeV
            if pid == 1000024:
                mtau *= 2 # the chargino decays into two taus
            if (protomodel.masses[pid] - protomodel.masses[1000022]) >= mtau: # If tau decays are allowed
                if lep_decays != 3:
                    ret *= 1/(4-lep_decays)
                    continue
            else:
                if lep_decays != 2:
                    ret *= 1/2
                    continue

            if lep_br != []:
                delta_br = numpy.max( [numpy.max(lep_br)-numpy.mean(lep_br), numpy.mean(lep_br)-numpy.min(lep_br)] )
                if abs(delta_br - 0) > 5e-03: 
                    self.highlight("warning", f"Undemocratic lep flavor decay for {pid}, delta_br = {delta_br}, br = {lep_br}")
                l = 1. - delta_br / 3.
                ret *= l

        return ret

    def penaltyForExtremeSSMs ( self, protomodel ) -> float:
        """ very simple hack for now, penalize for unusual values of
        the signal strength multipliers, where unusual means orders of magnitudes
        different from unity.

        :returns: penalty -- 1 - log10(ssm)/20. per ssm
        """
        ret = 1.
        for k,v in protomodel.ssmultipliers.items():
            if v == 0. or 0.01 < v < 100.:
                continue
            if v < 1.:
                ## 0.01 is 0.9, 0.001 is 0.85, 1e-6 is 0.7
                l = 1 + numpy.log10 ( v ) / 10.
            if v > 1.:
                ## 100 is 0.9, 1000 is 0.85, 1e6 is 0.7
                l = 1 - numpy.log10 ( v ) / 10.
            if l < 1e-10:
                l = 1e-10
            ret *= l
        return ret

    def computeK ( self, TL : float, prior : float ) -> float:
        """ compute K from TL and prior (simple) """
        prior = max ( 1e-200, prior ) # avoid infs
        return TL + 2* numpy.log ( prior )

    def getPredictionID ( self, prediction : TheoryPrediction ):
        """ construct a unique id of a prediction from the analysis ID,
            pids (v2) or smsList (v3), data set name """
        # FIXME check that smsList works
        return f"{prediction.analysisId()}:{prediction.dataId()!s}:{'; '.join(map(str, prediction.smsList))}"
        # v2
        #return "%s:%s:%s" % ( prediction.analysisId(), str(prediction.dataId()),
        #                      "; ".join(map(str,prediction.PIDs) ) )

    def sortPredictions ( self, predictions : List ) -> List:
        """ filter the predictions and sort according to decreasing L_BSM/L_SM. Return the
            sorted list of filtered predictions

        :param predictions: list of theory predictions
        :returns: sorted list of filtered theory predictions
        """
        sorted_pred = {}

        for pred in predictions:
            l0 = pred.likelihood ( 0. )    #SM mu = 0
            l1 = pred.likelihood ( 1. )    #BSM mu = 1
            # llhd_ratio = -1.
            if type(l0)!=type(None) and type(l1)!=type(None):
                llhd_ratio = l1/l0
                diff = abs(llhd_ratio - 1.0)
                if diff < 1e-05: continue
            else: continue   #TP: If computation of l0 or l1 gave None, do not consider the SR
            sorted_pred[pred] = llhd_ratio

        #sort acc to decreasing llhd ratio
        sorted_pred = dict(sorted(sorted_pred.items(), key = lambda item:item[1], reverse=True))
        newpreds = list ( sorted_pred.keys() )

        return newpreds


    def getMostSignificantCombination(self, predictions : List[TheoryPrediction], use_pathfinder=True) -> Tuple:
        """
            Gets the most significant combination and its corresponding weight (-2 ln L0/L1 for the whole combination)
            given the list of theory predictions.
        """
        self.log(f"Finding most significant combination for {len(predictions)}")
        comb_dict = bamAndWeights(predictions, evaluationType=observed)        
        #get the true/false comb matrix, along with weights
        pred_dict = comb_dict['theoryPred']                     #a dict with tpId and correspond tpred

        if use_pathfinder: most_significant_comb_dict = find_best_comb(comb_dict)  #get the best combination given the matrix and weights

        else:
            from tester.alternate_pf import getBestComb
            most_significant_comb_dict = getBestComb(predictions, evaluationType = observed)

        comb_lbl, weight = most_significant_comb_dict['best'], most_significant_comb_dict['weight']

        #Removing Jamie's penalty for now
        if False and len(comb_lbl)>1:
            weight = weight / math.sqrt(len(comb_lbl) - 1) # Rescale to have all the combinations on the same footing

        #convert best labels to theorypreds
        tp_comb = []
        for lbl in comb_lbl:
            tp = pred_dict[lbl]
            tp_comb.append(tp)

        return tp_comb, weight

    def getMostSensitiveCombination(self, predictions : List[TheoryPrediction], use_pathfinder=True) -> Tuple:
        """
            Gets the most sensitive combination and its corresponding weight (-2 ln L1/L0 (expected likelihoods) for the whole combination)
            given the list of theory predictions.
        """
        self.log(f"Finding most sensitive combination for {len(predictions)}")
        comb_dict = bamAndWeights(predictions, evaluationType=apriori, 
                excl_mode=True)        
        # get the true/false comb matrix, along with weights
        pred_dict = comb_dict['theoryPred']                     
        # a dict with tpId and correspond tpred

        if use_pathfinder: most_sensitive_comb_dict = find_best_comb(comb_dict)  
        # get the best combination given the matrix and weights
        else:
            from tester.alternate_pf import getBestComb
            most_sensitive_comb_dict = getBestComb(predictions, evaluationType = apriori )

        comb_lbl, weight = most_sensitive_comb_dict['best'], most_sensitive_comb_dict['weight']

        #convert best labels to theorypreds
        tp_comb = []
        for lbl in comb_lbl:
            tp = pred_dict[lbl]
            tp_comb.append(tp)

        return tp_comb, weight

    def getMuhat(self, predictions : List[TheoryPrediction])->float:
        tpredcomb = TheoryPredictionsCombiner(predictions)
        return tpredcomb.muhat()

    def findHighestSignificance ( self, predictions : List[TheoryPrediction],
            expected : bool =False ) -> Tuple:
        """ for the given list of predictions find the combination with highest significance

        :param predictions: list of theory predictions
        :param expected: find the highest expected significance, not observed
        :returns: tuple with most significant combination (msc), TL of msc (-2log(L0/L1)), muhat of msc
        """

        #filter for most significant SRs
        # min_rel_weight = 0.05 ## default
        min_rel_weight = 0.10
        filtered_preds = selectMostSignificantSRs(predictions,
                        min_rel_weight=min_rel_weight)
        self.letters = self.getLetters ( filtered_preds )

        self.log(f"Filtered predictions from {len(predictions)} to {len(filtered_preds)} (min_rel_weight {min_rel_weight:.2f})")

        most_significant_comb, weight = self.getMostSignificantCombination(filtered_preds)
        import time
        start_time = time.time()
        muhat = self.getMuhat(most_significant_comb)
        if muhat: muhat = float(muhat)
        end_time = time.time()
        time_taken = end_time - start_time
        self.log(f"Computed muhat, taken {time_taken:.3f} seconds")
        TL = float(weight)

        return most_significant_comb,TL,muhat

    def removeDataFromBestCombo ( self, bestCombo ):
        """ remove the data from all theory predictions, we dont need them. """
        self.debug ( "removing Data from best Combo " )
        for ci,combo in enumerate(bestCombo):
            if hasattr ( combo, "elements" ):
                del bestCombo[ci].elements
            if hasattr ( combo, "avgElement" ):
                del bestCombo[ci].avgElement
            eR = bestCombo[ci].expResult
            for ds in eR.datasets:
                for tx in ds.txnameList:
                    if hasattr ( tx, "txnameData" ):
                        del tx.txnameData
                    if hasattr ( tx, "txnameDataExp" ):
                        del tx.txnameDataExp
        return bestCombo

    def removeDataFromTheoryPred ( self, tp ):
        """ remove unnecessary stuff from a theoryprediction object.
            for storage. """
        self.debug ( f"removing data from theory pred {tp.analysisId()}" )
        theorypred = copy.deepcopy( tp )
        if hasattr ( theorypred, "elements" ):
            del theorypred.elements
        if hasattr ( theorypred, "avgElement" ):
            del theorypred.avgElement
        eR = theorypred.expResult
        for ds in eR.datasets:
            for tx in ds.txnameList:
                if hasattr ( tx, "txnameData" ):
                    del tx.txnameData
                if hasattr ( tx, "txnameDataExp" ):
                    del tx.txnameDataExp
        return theorypred

def normalizePrior():
    c = Combiner()
    S=0.
    ctr,nmod=0,30
    control = 0
    for nparticles in range ( 1, 18 ):
        for nbr in range ( 0, 10*nparticles ):
            for nssms in range ( 1, 25*nparticles ):
                t = c.priorForNDF ( nparticles, nbr, nssms, 1. )
                ctr+=1
                control += c.priorForNDF ( nparticles, nbr, nssms, None )
                if ctr % nmod == 0:
                    print ( f"nparticles {int(nparticles)}, nbr {int(nbr)}, nssms {int(nssms)}, improper prior {t:.5f}" )
                    nmod=nmod*2
                S += t
    print ( f"The constant for normalizing the prior is {1.0 / S:.8f}" )
    print ( "With the current normalization we get", control )
    return 1./S

if __name__ == "__main__":
    if False:
        from smodels.base import runtime
        runtime._experimental = True
    import argparse
    argparser = argparse.ArgumentParser(
            description='combiner. if called from commandline, computes the highest Z' )
    argparser.add_argument ( '-f', '--slhafile',
            help='slha file to test [test.slha]',
            type=str, default="test.slha" )
    argparser.add_argument ( '-d', '--database',
            help='path to database [<rundir>/database.pcl]',
            type=str, default="<rundir>/database.pcl" )
    argparser.add_argument ( '-u', '--upper_limits',
            help='use only upper limits results', action='store_true' )
    argparser.add_argument ( '-e', '--efficiencyMaps',
            help='use only efficiency maps results', action='store_true' )
    argparser.add_argument ( '-E', '--expected',
            help='expected values, not observed', action='store_true' )
    argparser.add_argument ( '-P', '--prior',
            help='Compute normalization constant for prior, then quit', action='store_true' )
    args = argparser.parse_args()
    if args.prior:
        normalizePrior()
        sys.exit()
    if args.upper_limits and args.efficiencyMaps:
        print ( "[combiner] -u and -e are mutually exclusive" )
        sys.exit()
    from smodels.experiment.databaseObj import Database
    from smodels.decomposition import decomposer
    from smodels.share.models.SMparticles import SMList
    from smodels.base.model import Model
    from smodels.base.physicsUnits import fb
    from tester.combinationsmatrix import getYamlMatrix

    model = Model(BSMparticles=BSMList, SMparticles=SMList)
    model.updateParticles(inputFile=args.slhafile)

    combinationsmatrix, status = getYamlMatrix()
    if not combinationsmatrix or status != 0:
        sys.exit("Combination matrix not loaded correctly.")

    print ( "[combiner] loading database", args.database )
    db = Database ( args.database, combinationsmatrix = combinationsmatrix )
    print ( "[combiner] done loading database" )

    anaIds = [ "CMS-SUS-16-033" ]
    anaIds = [ "all" ]
    dts = [ "all" ]
    if args.upper_limits:
        dts = [ "upperLimit" ]
    if args.efficiencyMaps:
        dts = [ "efficiencyMap" ]
    listOfExpRes = db.getExpResults( analysisIDs = anaIds, dataTypes = dts,
                                     onlyWithExpected= True )
    smses = decomposer.decompose ( model, .01*fb )
    #print ( "[combiner] decomposed into %d topos" % len(smses) )
    from smodels.matching.theoryPrediction import theoryPredictionsFor
    combiner = Combiner()

    preds = theoryPredictionsFor ( db, smses )
    combo,globalTL, muhat = combiner.findHighestSignificance ( preds, expected=args.expected )
    print ( f"[combiner] global TL is {globalTL:.2f}: {combiner.getComboDescription(combo)} (muhat={muhat:.2f})" )

    for pred in preds:
        pred.computeStatistics()
        tpe = pred.dataType(True)
        tpe += f":{','.join(map(str, pred.txnames))}"
        print ( f"  `- llhd [{tpe}] SM={pred.likelihood(0.0, evaluationType=args.expected):.3g} BSM={pred.likelihood(1.0, evaluationType=args.expected):.3g}" )
    comb = Combiner()
