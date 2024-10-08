#!/usr/bin/env python3

"""old combiner methods used in the previous run"""

from smodels.decomposition import decomposer
from smodels.share.models.SMparticles import SMList
from share.model_spec import BSMList
from smodels.base.physicsUnits import fb
from smodels.base.model import Model
from smodels.matching.theoryPrediction import TheoryPrediction, TheoryPredictionsCombiner
import sys, os
sys.path.insert(0,os.path.abspath ( os.path.dirname(__file__) ) )
from base.loggerbase import LoggerBase
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
    
    
def removeDataType ( self, predictions, dataType ):
    """ remove from the predictions all the ones that match dataType """
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

def findHighestSignificance ( self, predictions : List[TheoryPrediction], expected : bool =False ) -> Tuple:
    #All old code, might not work properly!
    predictions = self.selectMostSignificantSR ( predictions )
    predictions = self.sortPredictions ( predictions )
    #!old code
    combinables = self.findCombinations ( predictions, strategy )
    singlepreds = [ [x] for x in predictions ]
    #optionally, add individual predictions ... nah!!
    #combinables = singlepreds + combinables
    self.discussCombinations ( most_significant_comb )
    bestCombo,Z,muhat = self._findLargestZ ( combinables, expected=expected, mumax = mumax )
    bestCombo = sorted( bestCombo, key = lambda tp: tp.expResult.globalInfo.id )
    # compute a likelihood equivalent for Z
    if Z is not None:
        llhd = stats.norm.pdf(Z)
    else:
        llhd = None
    print ( "bestCombo %s, %s, %s " % ( Z, llhd, muhat ) )
