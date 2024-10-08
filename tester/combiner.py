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
import numpy, math, copy, sys, time
from colorama import Fore
from typing import List, Union, Tuple, Set
from ptools.helpers import getAllPidsOfTheoryPred
from ptools.bamCreator import selectMostSignificantSRs, bamAndWeights, find_best_comb

class Combiner ( LoggerBase ):
    def __init__ ( self, walkerid: int=0 ):
        super ( Combiner, self ).__init__ ( walkerid )
        self.walkerid = walkerid

    def getAllPidsOfCombo ( self, combo : List[TheoryPrediction] ) -> Set:
        """ get all particle IDs (PIDs) that make it into one combo """
        pids = set()
        for theoryPred in combo:
            pids = pids.union ( getAllPidsOfTheoryPred ( theoryPred ) )
        return pids

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
            self.pprint ( "%.2f: %.3g" % ( k, v ) )

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

    def getLetterCode ( self, combination ):
        """ get the letter code of a combination """
        ret = ""
        for c in combination:
            ret += self.letters[c]
        return ret

    def priorForNDF ( self, nparticles, nbranchings, nssms, name="expo1", verbose=False, nll=False ):
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
            if abs ( cbr - br ) / ( cbr + br ) < 0.025: ## 5 percent rule
                return False
        return True

    def computePrior ( self, protomodel, nll : bool =False, verbose : bool =False, name : str ="expo1" ) -> float:
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
            self.log ( "           `- the unique ssms are: %s" % ", ".join ( ssmstring ) )
        if nll:
            return - math.log ( ret )
        return ret

    def penaltyForMissingResults ( self, predictions : List[TheoryPrediction] ) -> float:
        """ very simple hack for now, penalize if predictions are all from the same experiment
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
            if len (decays) == 0:
                continue

            if all( [ len(decay) == 2 for decay in decays.keys() ] ):  # On-shell, continue
                continue
            if not all( [ len(decay) == 3 for decay in decays.keys() ] ): # If both on- and off-shell decays, disallow
                from termcolor import colored
                self.pprint(colored(f"pID {pid} of mass {protomodel.masses[pid]} has both on- AND off-shell decays: {decays}! Returning a prior of 0.", "red"))
                return 0
            if len ( decays ) == 1:
                if list(decays.keys())[0] in [(1000022, 2, 1), (1000022, 2, 2), (1000022, 5, 5)]: # If the open channel is for light quarks or bb, allow
                    continue

            lep_decays = 0
            C1_has_lep = False
            for decay in decays.keys():
                if pid == 1000024 and decay in [(1000022, 11, 12), (1000022, 13, 14), (1000022, 15, 16)] and decays[decay] != 0.:
                    lep_decays += 1
                    continue
                if pid == 1000023 and decay in [(1000022, 11, 11), (1000022, 13, 13), (1000022, 15, 15)]:
                    if decays[decay] == 0.:
                        C1_has_lep = True # Chargino 1 has leptonic decays but set to 0
                    else:
                        lep_decays += 1
                    continue
            if lep_decays == 0: # Multiple open channels but no leptonic one, ok for neutralino 2, disallow for chargino 1
                if pid == 1000023:
                    continue
                elif pid == 1000024 and not C1_has_lep:
                    self.pprint(f"Chargino 1 has multiple open channels but no leptonic one: {decays}! Returning a prior of 0.")
                    return 0

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

            br = [decays[decay] for decay in decays]
            delta_br = numpy.max( [numpy.max(br)-numpy.mean(br), numpy.mean(br)-numpy.min(br)] )
            l = 1. - delta_br / 3.
            ret *= l

            # Old version
            # ## is there an electron decay?
            # if ( 1000022, 11 ) in decays and decays[(1000022,11)] not in [ 0., 1.]:
            #     if (1000022, 13) in decays:
            #         delta_br = abs ( decays[(1000022,11)]-decays[(1000022,13)] )
            #         l = 1. - delta_br / 5.
            #         ret *= l
            #     if (1000022, 15) in decays:
            #         delta_br = abs ( decays[(1000022,11)]-decays[(1000022,15)] )
            #         l = 1. - delta_br / 20.
            #         ret *= l
            # elif ( 1000022, 13 ) in decays and decays[(1000022,13)] not in [ 0., 1.]:
            #     # no electron, but muon and tau!
            #     if (1000022, 15) in decays:
            #         delta_br = abs ( decays[(1000022,13)]-decays[(1000022,15)] )
            #         l = 1. - delta_br / 20.
            #         ret *= l
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
        # if len(comb_lbl)>1:
        #    weight = weight / math.sqrt(len(comb_lbl) - 1) # Rescale to have all the combinations on the same footing

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
        comb_dict = bamAndWeights(predictions, expected=True, excl_mode=True)        #get the true/false comb matrix, along with weights
        pred_dict = comb_dict['theoryPred']                     #a dict with tpId and correspond tpred

        if use_pathfinder: most_sensitive_comb_dict = find_best_comb(comb_dict)  #get the best combination given the matrix and weights
        else:
            from tester.alternate_pf import getBestComb
            most_sensitive_comb_dict = getBestComb(predictions, expected=True)
            
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

    def findHighestSignificance ( self, predictions : List[TheoryPrediction], expected : bool =False ) -> Tuple:
        """ for the given list of predictions find the combination with highest significance

        :param predictions: list of theory predictions
        :param expected: find the highest expected significance, not observed
        :returns: most significant combination (msc), TL of msc (-2log(L0/L1)), muhat of msc
        """

        #filter for most significant SRs
        filtered_preds = selectMostSignificantSRs(predictions)
        self.letters = self.getLetters ( filtered_preds )

        self.log(f"Filtered predictions from {len(predictions)} to {len(filtered_preds)}")

        most_significant_comb, weight = self.getMostSignificantCombination(filtered_preds)

        muhat = self.getMuhat(most_significant_comb)
        TL = weight

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
    combo,globalTL, muhat = combiner.findHighestSignificance ( preds, expected=args.expected )
    print ( "[combiner] global TL is %.2f: %s (muhat=%.2f)" % (globalTL, combiner.getComboDescription(combo),muhat ) )

    for pred in preds:
        pred.computeStatistics()
        tpe = pred.dataType(True)
        tpe += ":" + ",".join ( map ( str, pred.txnames ) )
        print ( "  `- llhd [%s] SM=%.3g BSM=%.3g" % ( tpe, pred.likelihood(0.,expected=args.expected), pred.likelihood(1.,expected=args.expected) ) )
    comb = Combiner()
