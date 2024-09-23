#!/usr/bin/env python3

""" Try out combinations from pickle file. """

__all__ = [ "Combiner" ]

from smodels.decomposition import decomposer
from smodels.share.models.SMparticles import SMList
from share.model_spec import BSMList
from smodels.base.physicsUnits import fb
from smodels.base.model import Model
from smodels.matching.theoryPrediction import TheoryPrediction, TheoryPredictionsCombiner
import sys, os
sys.path.insert(0,os.path.abspath ( os.path.dirname(__file__) ) )
from base.loggerbase import LoggerBase
from typing import Set

try:
    from tester import analysisCombiner
except:
    import analysisCombiner
import numpy, math, copy, sys, time
from colorama import Fore
from scipy import optimize, stats
from scipy.special import erf
from typing import List, Union, Tuple, Set
from ptools.helpers import getAllPidsOfTheoryPred
from ptools.bamCreator import selectMostSignificantSRs, bamAndWeights, find_best_comb
# import IPython

class Combiner ( LoggerBase ):
    def __init__ ( self, walkerid: int=0 ):
        super ( Combiner, self ).__init__ ( walkerid )
        self.walkerid = walkerid

    def getAllPidsOfCombo ( self, combo : List[TheoryPrediction] ) -> Set:
        """ get all PIDs that make it into one combo """
        pids = set()
        for theoryPred in combo:
            pids = pids.union ( getAllPidsOfTheoryPred ( theoryPred ) )
        return pids

    def getAllPidsOfTheoryPred ( self, theorypred ):
        print ( "FIXME the method combiner.getAllPidsOfTheoryPred is obsolete, call ptools.helpers.getAllPidsOfTheoryPred directly!!" )
        return getAllPidsOfTheoryPred (theorypred )

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

    def findCompatibles ( self, predA, predictions, strategy ):
        """ return list of all elements in predictions
            combinable with predA, under the given strategy """
        ret = []
        n=len(predictions)
        for ct,i in enumerate(predictions):
            if analysisCombiner.canCombine ( predA, i ):
                lpredA, li = predA, i
                if type(predA)!=list:
                    lpredA = [ predA ]
                if type(i)!=list:
                    li = [ i ]
                combo = lpredA + li
                ret.append ( combo )
                if ct < n:
                    deeper = self.findCompatibles ( combo, predictions[ct+1:], strategy )
                    for d in deeper:
                        ret.append ( d )
        return ret

    def removeDataType ( self, predictions, dataType ):
        """ remove from the predictions all the ones
        that match dataType """
        if predictions is None:
            return predictions
        tmp = []
        for pred in predictions:
            if pred.dataType() == dataType:
                continue
            tmp.append ( pred )
        self.pprint ( "removed %s, %d/%d remain" % \
                     ( dataType, len(tmp), len(predictions) ) )
        return tmp

    def findCombinations ( self, predictions : List, strategy : str ) -> List:
        """ finds all allowed combinations of predictions, for
            the given strategy
        :param predictions: list of predictions
        :returns: a list of combinations
        """
        self.log ( f"now find all combos of {len(predictions)} preds!" )
        combinables=[]
        """
        for iA,predA in enumerate(predictions):
            nexti = iA + 1
            compatibles = self.findCompatibles ( predA, predictions[nexti:], strategy )
            combinables += compatibles
        """
        from tester.combinationFinder import CombinationFinder
        comb_ob = CombinationFinder()
        combinables = comb_ob.getPossibleCombinations(predictions)
        #aids = [[a.analysisId() for a in comb] for comb in combinables]
        #print("possible combinations = ", aids)

        self.log ( f"found {len(combinables)} combos!" )
        return combinables

    def discussCombinations ( self, combinables ):
        """ simple method that writes some stats about a combination to the log file """
        count={}
        for i in combinables:
            n =len(i)
            if not n in count.keys():
                count[n]=0
            count[n]+=1
        npred = 0
        if 1 in count.keys():
            npred = count[1]
        self.debug ( "%d combinations from %d predictions" % \
                     (len(combinables),npred) )

    def getCombinedLikelihood ( self, combination, mu, expected=False, nll=False ):
        """ get the combined likelihood for a signal strength mu
        :param nll: compute the negative log likelihood
        """
        llhds = numpy.array ( [ c.likelihood(float(mu),expected=expected) for c in combination ], dtype=object )
        ret = numpy.prod ( llhds[llhds!=None] )
        if nll:
            if ret <= 0.:
                ret = 1e-70
            ret = - math.log ( ret )
        return ret

    def printLLhds ( self, llhds ):
        keys = list ( llhds.keys() )
        keys.sort()
        for k in keys:
            v=llhds[k]
            self.pprint ( "%.2f: %.3g" % ( k, v ) )

    def isSubset ( self, small, big ):
        """ is the small combo a subset of the big combo? """
        for s in small:
            if not s in big:
                return False
        return True

    def isSubsetOf ( self, small, combos ):
        """ is the small combo already a subset of any of the
            combos in 'combos'? """
        for c in combos:
            if self.isSubset ( small, c ):
                return True
        return False

    def sortOutSubsets ( self, combinations ):
        """ of all given combinations, sort out all those
            that are subsets of larger combinations.
        :returns: combinations without subsets
        """
        combinations.sort ( key=len, reverse=True ) ## sort them first be length
        ret = []
        for c in combinations:
            if self.isSubsetOf ( c, ret ):
                # self.pprint ( f"{getLetterCode(c)} is subset of bigger combo. skip." )
                continue
            ret.append ( c )
        self.pprint ( f"sorting out subsets, reduced {len(combinations)} -> {len(ret)} combinations." )
        return ret

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
            return "%s(%s)" % ( x.analysisId(), x.dataType().replace("upperLimit", "ul" ).replace ( "efficiencyMap", "em" ).replace ( "combined", "comb" ) )
        return ",".join( [ describe(x) for x in combination ] )

    def getSignificance ( self, combo, expected=False, mumax=None ):
        """ obtain the significance of this combo
        :param expected: get the expected significance, not observed
        :param mumax: maximum muhat before we run into exclusions
        :returns: Z (significance) and muhat ( signal strength multiplier that maximizes Z)
        """
        if len(combo)==0.:
            return 0.,0.

        muhat = self.findMuHat ( combo ) #Compute exactly
        """
        effCombo = [tp.dataset.dataInfo.expectedBG for tp in combo
                    if tp.dataset.dataInfo.dataType == 'efficiencyMap']

        minEvt = 0.
        if len(effCombo) > 0:
            #Check if gaussian approximation is valid for the likelihood:
            minEvt = min(effCombo)

        if minEvt > 15: #15 events keeps the approximation error ~under 20%
            muhat = self.findMuHatApprox ( combo ) #Use gaussian approximation
        else:
            muhat = self.findMuHat ( combo ) #Compute exactly
        """

        if mumax is None:
            mumax = float("inf")
        if muhat is None:
            return 0.,0.
        if muhat > mumax:
            self.debug ( "muhat(%.2f) > mumax(%.2f). use mumax" % ( muhat, mumax ) )
            muhat = mumax
        l0 = numpy.array ( [ c.likelihood(0.,expected=expected) for c in combo ], dtype=object )
        LH0 = numpy.prod ( l0[l0!=None] )
        l1 = numpy.array ( [ c.likelihood(muhat,expected=expected) for c in combo ], dtype=object )
        LH1 = numpy.prod ( l1[l1!=None] )
        if LH0 <= 0.:
            self.highlight ( "debug", f"l(mu=0) was {LH0:.2g}. Set to 1e-80" )
            LH0 = 1e-80
        if LH1 <= 0.:
            self.highlight ( "debug", f"l(muhat={muhat:.3g}) was {LH1:.2g}. Set to 1e-80." )
            LH1 = 1e-80
        chi2 = 2 * ( math.log ( LH1 ) - math.log ( LH0 ) ) ## chi2 with one degree of freedom
        if chi2 < 0.:
            chi2 = 0.
        Z = numpy.sqrt ( chi2 )
        return Z, muhat

    def _findLargestZ ( self, combinations, expected=False, mumax=None ):
        """ find the combo with the most significant deviation
        :param expected: find the combo with the most significant expected deviation
        :param mumax: Maximum muhat to allow before we run into an exclusion
        """
        self.log ( "now find largest Z" )
        combinations = self.sortOutSubsets ( combinations )
        # combinations.sort ( key=len, reverse=True ) ## sort them first by length
        # compute CLsb for all combinations
        highestZ,highest,muhat=None,"",mumax
        ## we will not look at combos that are subsets.
        doProgress=True
        try:
            import progressbar
        except ModuleNotFoundError:
            doProgress=False
        if len(combinations)<10:
            doProgress = False
        if doProgress:
            pb = progressbar.ProgressBar(widgets=["combination #",progressbar.Counter(),
                  "/%d " % len(combinations),
                  progressbar.Percentage(),
                  progressbar.Bar( marker=progressbar.RotatingMarker() ),
                  progressbar.AdaptiveETA()])
            pb.maxval = len(combinations)
            pb.start()
        for ctr,c in enumerate(combinations):
            if doProgress:
                pb.update(ctr)
            Z,muhat_ = self.getSignificance ( c, expected=expected, mumax=mumax )
            # self.log ( f"combination #{ctr}: {Z}" )
            if Z == None:
                continue
            # self.pprint ( "[combine] significance for %s is %.2f" % ( self.getLetterCode(c), Z ) )
            if highestZ is None or Z > highestZ:
                highestZ = Z
                highest = c
                muhat = muhat_
        if doProgress:
            pb.finish()
        return highest,highestZ,muhat

    def get95CL ( self, combination, expected ):
        """ compute the CLsb value for one specific combination
        :param expected: compute expected instead of observed value
        """
        llhds={}
        muhat = self.findMuHat ( combination )
        if muhat == None:
            return None
        Lmuhat = self.getCombinedLikelihood ( combination, muhat, expected=expected )
        # mumin=muhat/3.
        # Lmumin = getCombinedLikelihood ( combination, mumin, expected=expected )
        mumin = 0.
        mumax = muhat
        while True:
            mumax=2.*mumax+0.5
            Lmumax = self.getCombinedLikelihood ( combination, mumax, expected=expected )
            if Lmumax / Lmuhat < 1e-3: ## less than 1 permille? stop!
                break
        dmu = ( mumax - mumin ) / 30.
        for mu in numpy.arange(mumin,mumax,dmu): ## scan mu
            L = self.getCombinedLikelihood ( combination, mu, expected=expected )
            llhds[mu]=L
        # self.printLLhds ( llhds )
        Sm = sum ( llhds.values() )
        C = 0.
        xold = 0.0
        for x,v in llhds.items():
            Cold = C
            C+=v/Sm
            if C>.95:
                k = v/Sm / ( x - xold )
                d = C - k*x
                return ( 0.95 - d ) / k
                # return xold + ( x - xold ) * ( C - Cold )
            xold = x
        return 1.

    def root_func (self, x, a ):
        """Auxiliary function for finding mumax when computing events from limits"""

        return (erf((a-x))+erf(x)) / ( 1. + erf(x)) - .95

    def mmaxFit(self,a):
        """Auxiliary function which gives the root (x) to root_func"""

        if a < 1.38:
            return 0.0
        elif a < 1.5:
            c0,c1,c2 = -2.04738814,  1.50637722, -0.02005543
        elif a < 2.:
            c0,c1,c2 = -2.58981055, 2.28286693, -0.2972251
        else:
            c0,c1,c2 = -1.20,1.0,0

        return c0 + c1*a + c2*a**2

    def getEventsFromLimits(self, upperLimit, expectedUpperLimit, maxdiff = 0.4 ):
        """Using the gaussian limit, extract the (normalized) number of expected BG events,
        and the number of observed events. The normalization is such that the corresponding
        number of signal events is the total production cross-section."""

        dr = ( expectedUpperLimit - upperLimit ) / ( expectedUpperLimit + upperLimit )
        if abs(dr)>maxdiff:
            self.pprint ("asking for likelihood from limit but difference between oUL(%.2f) and eUL(%.2f) is too large (dr=%.2f)" % ( upperLimit, expectedUpperLimit, dr ) )
            return None

        sigma_exp = expectedUpperLimit / 1.96 # the expected scale, eq 3.24 in arXiv:1202.3415
        if upperLimit < expectedUpperLimit:
            ## underfluctuation. mumax = 0., nobs=nbg
            #Approximation for number of BG events (gaussian limit):
            nbg = sigma_exp**2
            nobs = nbg
            return (nobs,nbg)

        #Rescale everything by denominator:
        denominator = numpy.sqrt(2.) * sigma_exp
        a = upperLimit/denominator
        if numpy.sign(self.root_func(0.,a)*self.root_func(a,a)) > 0.:
            self.pprint ( "when computing likelihood: fA and fB have same sign")
            return None
        xmax = self.mmaxFit(a)

        #Approximation for number of signal events (gaussian limit) which maximizes the likelihood:
        mumax = xmax*denominator
        #Approximation for number of BG events (gaussian limit):
        nbg = sigma_exp**2
        #Approximation for number of observed events (gaussian limit):
        nobs = mumax + nbg

        return (nobs,nbg)

    #!remove function later as already in SModelS
    def findMuHatApprox ( self, combination ):
        """ find the maximum likelihood estimate for the signal strength mu in the gaussian limit"""

        nobs = []
        nbg = []
        bgerr = []
        ns = []
        for tp in combination:
            if tp.dataset.getType() != "combined" and tp.dataset.dataInfo.dataType == 'upperLimit':
                upperLimit = tp.upperLimit.asNumber(fb)
                expectedUL = tp.getUpperLimit ( expected=True ).asNumber(fb)
                n = self.getEventsFromLimits(upperLimit,expectedUL)
                if not n: continue #Could not get events
                observedN,expectedBG = n
                bgError = 0.0 #In this case we ignore systematics
                Nsig = tp.xsection.asNumber(fb) #Arbitrary normalization (does not affect muhat)
                nobs.append(observedN)
                nbg.append(expectedBG)
                bgerr.append(bgError)
                ns.append(Nsig)
            if tp.dataset.getType() != "combined" and tp.dataset.dataInfo.dataType == "efficiencyMap":
                observedN = tp.dataset.dataInfo.observedN
                expectedBG = tp.dataset.dataInfo.expectedBG
                bgError = tp.dataset.dataInfo.bgError
                Nsig = tp.xsection*tp.expResult.globalInfo.lumi
                Nsig = Nsig.asNumber()
                nobs.append(observedN)
                nbg.append(expectedBG)
                bgerr.append(bgError)
                ns.append(Nsig)
            if tp.dataset.getType() == "combined":
                for d in tp.dataset._datasets:
                    observedN = d.dataInfo.observedN
                    expectedBG = d.dataInfo.expectedBG
                    bgError = d.dataInfo.bgError
                    Nsig = tp.xsection*tp.expResult.globalInfo.lumi
                    Nsig = Nsig.asNumber()
                    nobs.append(observedN)
                    nbg.append(expectedBG)
                    bgerr.append(bgError)
                    ns.append(Nsig)

        if not ns:
            return 1.0
        ns = numpy.array(ns)
        nobs = numpy.array(nobs)
        nbg = numpy.array(nbg)
        bgerr = numpy.array(bgerr)
        sigma2 = nbg+bgerr**2 #total error^2 (systematics + statistics)

        num = numpy.sum(ns*(nobs-nbg)/sigma2)
        if num <= 0:
            return 0.0
        den = numpy.sum(ns*ns/sigma2)

        mu = num/den
        return mu
    #!remove function later as already in SModelS
    def findMuHat ( self, combination ):
        """ find the maximum likelihood estimate for the signal strength mu """
        def getNLL ( mu ):
            ret = self.getCombinedLikelihood ( combination, mu, nll=True )
            return ret
        #Quick guess for muhat:
        mustart = self.findMuHatApprox(combination)
        for start in [mustart,0.,1.0,0.1,10.,1e-2,1e-3]:
            ret = optimize.minimize ( getNLL, start, bounds=[(0.,None)] )
            # print ( "findMuHat combo %s start=%f, ret=%s" % ( combination, start, ret.fun ) )
            if ret.status==0:
                return ret.x[0]
        self.pprint ( "%serror finding mu hat for %s%s" % (Fore.RED, self.getLetterCode(combination), Fore.RESET ) )
        return None

    def getLetterCode ( self, combination ):
        """ get the letter code of a combination """
        ret = ""
        for c in combination:
            ret += self.letters[c]
        return ret

    def priorForNDF ( self, nparticles, nbranchings, nssms, name="expo1",
                      verbose=False, nll=False ):
        """ get the prior for this and this many degrees of freedom
            in the model.
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
            self.pprint ( f"prior ``{name}'': {nparticles} particles, {nbranchings} branchings, {nssms:.1f} equivalent unique ssms: {prior}" )
        if nll:
            return - numpy.log ( prior )
        return prior

    def noSuchBranching ( self, branchings, br ):
        """ check if a branching ratio similar to br already exists
            in branchings """
        for cbr in branchings:
            if abs ( cbr - br ) / ( cbr + br ) < 0.025: ## 5 percent rule
                return False
        return True

    def computePrior ( self, protomodel, nll : bool =False,
        verbose : bool =False, name : str ="expo1" ) -> float:
        """ compute the prior for protomodel, used to introduce regularization,
            i.e. penalizing for non-zero parameters, imposing sparsity.

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
            memBRs = set() ## memorize branchings, similar branchings count only once
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
            self.pprint( "could not get xsecs for protomodel with masses %s" % \
                         str(protomodel.masses) )
            if nll:
                return -float("inf")
            return 0.
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
            ssmstring = [ "%.2f" % x for x in cssms.keys() ]
            self.pprint ( "           `- the unique ssms are: %s" % ", ".join ( ssmstring ) )
        if nll:
            return - math.log ( ret )
        return ret

    def computeK ( self, TL : float, prior : float ) -> float:
        """ compute K from TL and prior (simple) """
        return TL + 2* numpy.log ( prior )

    def getPredictionID ( self, prediction : TheoryPrediction ):
        """ construct a unique id of a prediction from the analysis ID,
            pids (v2) or smsList (v3), data set name """
        # FIXME check that smsList works
        return "%s:%s:%s" % ( prediction.analysisId(), str(prediction.dataId()),
                              "; ".join(map(str,prediction.smsList) ) )
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

    '''
    def selectMostSignificantSR ( self, predictions : List ) -> List:
        """ given, the predictions, for any analysis and topology,
            return the most significant SR only. FILTER PREDS
        :param predictions: all predictions of all SRs
        :returns: list of predictions of most significant SR of each analysis
        """
        sortByAnaId = {}                            # first sort all by ana id + data Type
        for pred in predictions:
            Id = pred.analysisId()+":"+pred.dataType(True)
            if not Id in sortByAnaId:
                sortByAnaId[Id]=[]
            sortByAnaId[Id].append ( pred )         #keep all em-type ds of one analysis under one key
        ret = []
        keptThese = [] ## log the ana ids that we kept, for debugging only.
        for Id,preds in sortByAnaId.items():
            maxRatio, bestpred = 0., None
            for pred in preds:
                oul = pred.getUpperLimit(expected=False)
                eul = pred.getUpperLimit(expected=True)
                if oul is None or eul is None:
                    continue
                ratio = oul / eul
                if ratio > maxRatio:
                    maxRatio = ratio
                    bestpred = pred
            if maxRatio > 0. and bestpred != None:
                ret.append ( bestpred )
                keptThese.append ( self.getPredictionID ( bestpred ) )
        self.pprint ( f"selected predictions down via SRs from {len(predictions)}"\
                      f" to {len(ret)}." )
        debug = False ## print the selections in debug mode
        if debug:
            for ctr,i in enumerate(predictions):
                ul,eul=i.getUpperLimit(),i.getUpperLimit(expected=True)
                r=float("nan")
                if type(eul) != type(None) and eul.asNumber(fb) > 0.:
                    r = ( ul / eul ).asNumber()
                didwhat = Fore.RED + "discarded"
                pId = self.getPredictionID( i )
                if pId in keptThese:
                    didwhat = "kept     "
                print(f" `- {didwhat}: #{ctr} {pId}: r={r:.2f}{Fore.RESET}")
        return ret
    '''

    def penaltyForMissingResults ( self,
            predictions : List[TheoryPrediction] ) -> float:
        """ very simple hack for now, penalize if predictions are all
        from the same experiment

        :returns: penalty -- 1e-3 if experiment is missing
        """
        if len(predictions)==0:
            return .01 # penalize but doesnt matter, there is no likelihood
        hasExperiment = { "ATLAS": False, "CMS": False }
        for p in predictions:
            for experiment in [ "CMS", "ATLAS" ]:
                if experiment in p.dataset.globalInfo.id:
                    hasExperiment[experiment]=True
        if hasExperiment["ATLAS"] and hasExperiment["CMS"]:
            return 1.
        self.pprint ( f"penalty! we only have {','.join( k for k,v in hasExperiment.items() if v )}" )
        return .1 # penalize!

    def penaltyForUndemocraticFlavors ( self, protomodel ) -> float:
        """ very simple hack for now, penalize for undemocratic flavor decays

        :returns: penalty -- 1 - delta(branching_ratios)/C for each flavor decay
        """
        ret = 1.
        for pid, decays in protomodel.decays.items():
            if not pid in [ 1000024, 1000023 ]:
                continue
            if len ( decays ) == 1: # only one flavor? allow!
                continue
            ## is there an electron decay?
            if ( 1000022, 11 ) in decays and decays[(1000022,11)] not in [ 0., 1.]:
                if (1000022, 13) in decays:
                    delta_br = abs ( decays[(1000022,11)]-decays[(1000022,13)] )
                    l = 1. - delta_br / 5.
                    ret *= l
                if (1000022, 15) in decays:
                    delta_br = abs ( decays[(1000022,11)]-decays[(1000022,15)] )
                    l = 1. - delta_br / 20.
                    ret *= l
            elif ( 1000022, 13 ) in decays and decays[(1000022,13)] not in [ 0., 1.]:
                # no electron, but muon and tau!
                if (1000022, 15) in decays:
                    delta_br = abs ( decays[(1000022,13)]-decays[(1000022,15)] )
                    l = 1. - delta_br / 20.
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
                l = 1 + numpy.log10 ( v ) / 20.
            if v > 1.:
                ## 100 is 0.9, 1000 is 0.85, 1e6 is 0.7
                l = 1 - numpy.log10 ( v ) / 20.
            if l < 1e-10:
                l = 1e-10
            ret *= l
        return ret

    def getMostSignificantCombination(self, predictions : List[TheoryPrediction], use_pathfinder=True) -> Tuple:
        """
            Gets the most significant combination and its corresponding weight (-2 ln L0/L1 for the whole combination)
            given the list of theory predictions.
        """
        self.log(f"Finding most significant combination for {len(predictions)}")
        comb_dict = bamAndWeights(predictions, expected=False)        #get the true/false comb matrix, along with weights
        pred_dict = comb_dict['theoryPred']                     #a dict with tpId and correspond tpred
        
        if use_pathfinder: most_significant_comb_dict = find_best_comb(comb_dict)  #get the best combination given the matrix and weights
 
        else:
            from tester.alternate_pf import getBestComb
            most_significant_comb_dict = getBestComb(predictions, expected=False)

        comb_lbl, weight = most_significant_comb_dict['best'], most_significant_comb_dict['weight']
        
        #Removing Jamie's penalty for now
        #if len(comb_lbl)>1:
        #    weight = weight / math.sqrt(len(comb_lbl) - 1) # Rescale to have all the combinations on the same footing
        
        #from ptools.helpers import experimentalId
        #tpred_lbl = {experimentalId(tpred):tpred for tpred in predictions}

        #convert best labels to theorypreds
        tp_comb = []
        for lbl in comb_lbl:
            tp = pred_dict[lbl]
            tp_comb.append(tp)

        return tp_comb, weight

    def getMostSensitiveCombination(self, predictions : List[TheoryPrediction]) -> Tuple:
        """
            Gets the most sensitive combination and its corresponding weight (-2 ln L1/L0 (expected likelihoods) for the whole combination)
            given the list of theory predictions.
        """
        self.log(f"Finding most sensitive combination for {len(predictions)}")
        comb_dict = bamAndWeights(predictions, expected=True, excl_mode=True)        #get the true/false comb matrix, along with weights
        pred_dict = comb_dict['theoryPred']                     #a dict with tpId and correspond tpred

        most_sensitive_comb_dict = find_best_comb(comb_dict)  #get the best combination given the matrix and weights
        comb_lbl, weight = most_sensitive_comb_dict['best'], most_sensitive_comb_dict['weight']

        #from ptools.helpers import experimentalId
        #tpred_lbl = {experimentalId(tpred):tpred for tpred in predictions}

        #convert best labels to theorypreds
        tp_comb = []
        for lbl in comb_lbl:
            tp = pred_dict[lbl]
            tp_comb.append(tp)

        return tp_comb, weight

    def getMuhat(self, predictions : List[TheoryPrediction])->float:
        tpredcomb = TheoryPredictionsCombiner(predictions)
        return tpredcomb.muhat()

    def findHighestSignificance ( self, predictions : List[TheoryPrediction], expected : bool =False ) -> Tuple:
        """ for the given list of predictions and employing the given combination strategy,
        find the combination with highest significance, SORT PREDICTIONS!

        :param predictions: list of theory predictions
        :param strategy: the combination strategy to use
        :param expected: find the highest expected significance, not observed
        :param mumax: maximimal signal strength mu that is allowed before we run
        into an exclusion
        :returns: most significant combination (msc), TL of msc (-2log(L0/L1)), muhat of msc
        """

        #for now, keep all SRs, let Pathfinder find best combination of all SRs
        #predictions = self.selectMostSignificantSR ( predictions )         !old method
        #predictions = self.sortPredictions ( predictions )                 !old method

        filtered_preds = selectMostSignificantSRs(predictions)
        self.letters = self.getLetters ( filtered_preds )

        self.pprint(f"Filtered predictions from {len(predictions)} to {len(filtered_preds)}")

        most_significant_comb, weight = self.getMostSignificantCombination(filtered_preds)

        muhat = self.getMuhat(most_significant_comb)
        TL = weight

        #!old code
        #combinables = self.findCombinations ( predictions, strategy )      !old method
        # singlepreds = [ [x] for x in predictions ]
        ## optionally, add individual predictions ... nah!!
        ## combinables = singlepreds + combinables
        #self.discussCombinations ( most_significant_comb )
        #bestCombo,Z,muhat = self._findLargestZ ( combinables, expected=expected, mumax = mumax )
        #bestCombo = sorted( bestCombo, key = lambda tp: tp.expResult.globalInfo.id )
        ## compute a likelihood equivalent for Z
        #if Z is not None:
        #    llhd = stats.norm.pdf(Z)
        #else:
        #    llhd = None
        # self.pprint ( "bestCombo %s, %s, %s " % ( Z, llhd, muhat ) )

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
        self.debug ( "removing data from theory pred %s" % tp.analysisId() )
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
                    print ( "nparticles %d, nbr %d, nssms %d, improper prior %.5f" % \
                            ( nparticles, nbr, nssms, t ) )
                    nmod=nmod*2
                S += t
    print ( "The constant for normalizing the prior is %.8f" % (1./S) )
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
    combo,globalTL, muhat = combiner.findHighestSignificance ( preds, "aggressive", expected=args.expected )
    print ( "[combiner] global TL is %.2f: %s (muhat=%.2f)" % (globalTL, combiner.getComboDescription(combo),muhat ) )


    preds = theoryPredictionsFor ( db, smses )
    Z, muhat_ = combiner.getSignificance ( preds, expected=args.expected, mumax = None )
    #print ( "%s has %d predictions, local Z is %.2f" % ( expRes.globalInfo.id, len(preds), Z ) )
    for pred in preds:
        pred.computeStatistics()
        tpe = pred.dataType(True)
        tpe += ":" + ",".join ( map ( str, pred.txnames ) )
        print ( "  `- llhd [%s] SM=%.3g BSM=%.3g" % ( tpe, pred.likelihood(0.,expected=args.expected), pred.likelihood(1.,expected=args.expected) ) )
    comb = Combiner()
