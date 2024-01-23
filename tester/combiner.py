#!/usr/bin/env python3

""" Try out combinations from pickle file. """

__all__ = [ "Combiner" ]

from smodels.theory import decomposer
from smodels.share.models.SMparticles import SMList
from smodels.particlesLoader import BSMList
from smodels.tools.physicsUnits import fb
from smodels.theory.model import Model
try:
    from tester import analysisCombiner
except:
    import analysisCombiner
import numpy, math, copy, sys, time
from colorama import Fore
from scipy import optimize, stats
from scipy.special import erf
from typing import List
# import IPython

class Combiner:
    def __init__ ( self, walkerid: int=0 ):
        self.walkerid = walkerid

    def getAllPidsOfCombo ( self, combo ):
        """ get all PIDs that make it into one combo """
        pids = set()
        for theoryPred in combo:
            pids = pids.union ( self.getAllPidsOfTheoryPred ( theoryPred ) )
        return pids

    def getAnaIdsWithPids ( self, combo, pids ):
        """ from best combo, retrieve all ana ids that contain *all* pids """
        anaIds = set()
        for theoryPred in combo:
            tpids = self.getAllPidsOfTheoryPred ( theoryPred )
            hasAllPids=True
            for pid in pids:
                if not pid in tpids:
                    hasAllPids=False
                    break
            if hasAllPids:
                anaIds.add ( theoryPred.analysisId() )
        return anaIds

    def getAllPidsOfTheoryPred ( self, pred ):
        """ get all pids that make it into a theory prediction """
        pids = set()
        for prod in pred.PIDs:
            for branch in prod:
                for pid in branch:
                    if type(pid) == list:
                        for p in pid:
                            pids.add ( abs(p) )
                    else:
                        pids.add ( abs(pid) )
        return pids

    def getTheoryPredsWithPids ( self, combo, pids ):
        """ from best combo, retrieve all theory preds that contain *all* pids """
        tpreds = set()
        for theoryPred in combo:
            tpids = self.getAllPidsOfTheoryPred ( theoryPred )
            hasAllPids=True
            for pid in pids:
                if not pid in tpids:
                    hasAllPids=False
                    break
            if hasAllPids:
                tpreds.add ( theoryPred )
        return tpreds

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
        for iA,predA in enumerate(predictions):
            nexti = iA + 1
            compatibles = self.findCompatibles ( predA, predictions[nexti:], strategy )
            combinables += compatibles
        self.log ( f"found {len(combinables)} combos!" )
        return combinables

    def highlight ( self, msgType : str = "info", *args ):
        """ logging, hilit """
        msgType = msgType.lower()
        if msgType == "debug":
            return
        col = Fore.GREEN
        if msgType == "error":
            col = Fore.RED
        print ( f"{col}[combine:{self.walkerid}] {' '.join(map(str,args))}{Fore.RESET}" )

    def error ( self, *args ):
        self.highlight ( "error", *args )

    def pprint ( self, *args ):
        """ logging """
        print ( "[combine:%d] %s" % (self.walkerid, " ".join(map(str,args))) )

    def log ( self, *args ):
        """ logging to file """
        with open( "walker%d.log" % self.walkerid, "a" ) as f:
            f.write ( "[combiner-%s] %s\n" % ( time.strftime("%H:%M:%S"), " ".join(map(str,args)) ) )

    def debug ( self, *args ):
        """ logging """
        self.highlight ( "debug", *args )

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
        letters={}
        if predictions is None:
            return letters
        ## assign a letter to every prediction. for debugging
        letter=65
        # self.pprint ( "[combine] Letters assigned to results:" )
        for p in predictions:
            letters[p]=chr(letter)
            # self.pprint ( "[combine] Prediction %s: %s" % ( letters[p], p.expResult.globalInfo.id ) )
            letter+=1
            if letter == 91: # skip the special characters
                letter = 97
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
            self.log ( f"combination #{ctr}: {Z}" )
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
                Nsig = tp.xsection.value.asNumber(fb) #Arbitrary normalization (does not affect muhat)
                nobs.append(observedN)
                nbg.append(expectedBG)
                bgerr.append(bgError)
                ns.append(Nsig)
            if tp.dataset.getType() != "combined" and tp.dataset.dataInfo.dataType == "efficiencyMap":
                observedN = tp.dataset.dataInfo.observedN
                expectedBG = tp.dataset.dataInfo.expectedBG
                bgError = tp.dataset.dataInfo.bgError
                Nsig = tp.xsection.value*tp.expResult.globalInfo.lumi
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
                    Nsig = tp.xsection.value*tp.expResult.globalInfo.lumi
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
        if verbose:
            self.pprint ( "prior ``%s'': %d particles, %d branchings, %.1f equivalent unique ssms: %.2f" % \
                      ( name, nparticles, nbranchings, nssms, prior ) )
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

    def computePrior ( self, protomodel, nll=False, verbose=False, name="expo1" ):
        """ compute the prior for protomodel, used to introduce regularization,
            i.e. punishing for non-zero parameters, imposing sparsity.
        :param nll: if True, return negative log likelihood
        :param verbose: print how you get the prior.
        :param name: name of prior (expo1, gauss1, etc). See self.priorForNDF.
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

    def computeK ( self, Z, prior ):
        """ compute K from Z and prior (simple) """
        return Z**2 + 2* numpy.log ( prior )

    def getPredictionID ( self, prediction ):
        """ construct a unique id of a prediction from the analysis ID,
            pids, data set name """
        return "%s:%s:%s" % ( prediction.analysisId(), str(prediction.dataId()),
                              "; ".join(map(str,prediction.PIDs) ) )

    def selectMostSignificantSR ( self, predictions ):
        """ given, the predictions, for any analysis and topology,
            return the most significant SR only.
        :param predictions: all predictions of all SRs
        :returns: filtered predictions
        """
        sortByAnaId = {} ## first sort all by ana id + data Type
        for i in predictions:
            Id = i.analysisId()+":"+i.dataType(True)
            if not Id in sortByAnaId:
                sortByAnaId[Id]=[]
            sortByAnaId[Id].append ( i )
        ret = []
        keptThese = [] ## log the ana ids that we kept, for debugging only.
        for Id,preds in sortByAnaId.items():
            maxR, bestpred = 0., None
            for pred in preds:
                eul = pred.getUpperLimit(expected=True)
                if eul is None:
                    continue
                r = pred.getUpperLimit() / eul
                if r > maxR:
                    maxR = r
                    bestpred = pred
            if maxR > 0. and bestpred != None:
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

    def findHighestSignificance ( self, predictions, strategy, expected=False,
                                  mumax = None ):
        """ for the given list of predictions and employing the given strategy,
        find the combo with highest significance

        :param expected: find the highest expected significance, not observed
        :param mumax: maximimal signal strength mu that is allowed before we run 
        into an exclusion
        :returns: best combination, significance (Z), likelihood equivalent
        """
        predictions = self.selectMostSignificantSR ( predictions )
        self.letters = self.getLetters ( predictions )
        combinables = self.findCombinations ( predictions, strategy )
        # singlepreds = [ [x] for x in predictions ]
        ## optionally, add individual predictions ... nah!!
        ## combinables = singlepreds + combinables
        self.discussCombinations ( combinables )
        bestCombo,Z,muhat = self._findLargestZ ( combinables, expected=expected,
                                                 mumax = mumax )
        bestCombo = sorted( bestCombo, key = lambda tp: tp.expResult.globalInfo.id )
        ## compute a likelihood equivalent for Z
        if Z is not None:
            llhd = stats.norm.pdf(Z)
        else:
            llhd = None
        # self.pprint ( "bestCombo %s, %s, %s " % ( Z, llhd, muhat ) )
        return bestCombo,Z,llhd,muhat

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
        from smodels.tools import runtime
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
    from smodels.theory import decomposer
    from smodels.particlesLoader import BSMList
    from smodels.share.models.SMparticles import SMList
    from smodels.theory.model import Model
    from smodels.tools.physicsUnits import fb
    model = Model(BSMparticles=BSMList, SMparticles=SMList)
    model.updateParticles(inputFile=args.slhafile)
    print ( "[combiner] loading database", args.database )
    db = Database ( args.database )
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
    from smodels.theory.theoryPrediction import theoryPredictionsFor
    combiner = Combiner()
    allps = []
    for expRes in listOfExpRes:
        preds = theoryPredictionsFor ( expRes, smses )
        if preds == None:
            continue
        for pred in preds:
            allps.append ( pred )
    combo,globalZ,llhd,muhat = combiner.findHighestSignificance ( allps, "aggressive", expected=args.expected )
    print ( "[combiner] global Z is %.2f: %s (muhat=%.2f)" % (globalZ, combiner.getComboDescription(combo),muhat ) )
    for expRes in listOfExpRes:
        preds = theoryPredictionsFor ( expRes, smses )
        if preds == None:
            continue
        Z, muhat_ = combiner.getSignificance ( preds, expected=args.expected, mumax = None )
        print ( "%s has %d predictions, local Z is %.2f" % ( expRes.globalInfo.id, len(preds), Z ) )
        for pred in preds:
            pred.computeStatistics()
            tpe = pred.dataType(True)
            tpe += ":" + ",".join ( map ( str, pred.txnames ) )
            print ( "  `- llhd [%s] SM=%.3g BSM=%.3g" % ( tpe, pred.getLikelihood(0.,expected=args.expected), pred.getLikelihood(1.,expected=args.expected) ) )
    comb = Combiner()
