#!/usr/bin/env python3

"""
the critic class, i.e. the class that tells if the
current model survives LHC constraints or not.
"""

__all__ = [ "Critic" ]

import pickle, time, os, sys
from smodels.decomposition import decomposer
from smodels.matching.theoryPrediction import theoryPredictionsFor, TheoryPrediction, TheoryPredictionList, TheoryPredictionsCombiner
sys.path.insert(0,"../")
from builder.protomodel import ProtoModel
from smodels.share.models.SMparticles import SMList
from share.model_spec import BSMList
from smodels.base.physicsUnits import fb, GeV, TeV
from smodels.experiment.databaseObj import Database
from smodels.base.model import Model
from smodels.base.exceptions import SModelSBaseError as SModelSError
from os import PathLike
from typing import List, Union
#from smodels.base.smodelsLogging import logger
from base.loggerbase import LoggerBase
from tester.combiner import Combiner
from tester.combinationsmatrix import getYamlMatrix

class Critic ( LoggerBase ):
    def __init__ ( self, walkerid : int, dbpath : PathLike = "official", expected : bool = False, select : str = "all", do_srcombine : bool = False ):
        #call the super class of the critic i.e Loggerbase
        super ( Critic, self ).__init__ ( walkerid )
        self.walkerid = walkerid
        self.do_srcombine = do_srcombine
        self.r_threshold = 1.38
        self.sensitivity_threshold = 0.7
        self.verbose = 1

        force_load = None
        if dbpath.endswith ( ".pcl" ):
            force_load = "pcl"

        combinationsmatrix, status = getYamlMatrix()
        if not combinationsmatrix or status != 0:
            sys.exit("Combination matrix not loaded correctly when instantiating Critic class.")

        self.database = Database( dbpath, force_load = force_load, combinationsmatrix = combinationsmatrix )
        self.combiner = Combiner(self.walkerid)

    # def getMaxAllowedMu(self, protomodel):
    #     """ Compute the maximum (global) signal strength normalization
    #         given the predictions.
    #     """
    #
    #     mumax = float("inf")
    #     if protomodel.ul_type_rvalues[0] > 0.:
    #         #Set mumax slightly below threshold, so the model is never excluded
    #         print(f"r value from UL: {protomodel.ul_type_rvalues[0]}")
    #         mumax_ul = 0.999*self.rthreshold / protomodel.ul_type_rvalues[0]
    #         #if 0 < protomodel.llhd_type_rvalue < mumax:
    #         print(f"r value from llhd: {protomodel.llhd_type_rvalue}")
    #         mumax_llhd = 0.999*self.rthreshold /protomodel.llhd_type_rvalue
    #         mumax = min(mumax_ul, mumax_llhd)
    #
    #     return mumax


    def updateModelPredictionsWithULPreds(self, protomodel, predictions, keep_predictions):
        """ Extract information from list of theory predictions and store list of dict with r_obs,
            r_exp and theory prediction(sorted according to decreasing r_obs values) in the protomodel.
            Also store description about the critic_tp in the protomodel.

            :param predictions: theory predictions for UL-type results (EM-type if no UL-type)
        """

        if not predictions: # If empty list
            protomodel.critic_description = "No dataset for fast critic."
            if keep_predictions:
                protomodel.ul_critic_robs = None
                protomodel.ul_critic_tpList = None
            return

        robs = [] # If there are no predictions set rmax and r2 to 0
        tpList = []

        for theorypred in predictions:
            r = theorypred.getRValue(expected=False)

            if r is None:
                self.highlight("warning","The computation of the observed r-value of the most sensitive combination gave None.")
                r = 20 # Something is wrong, we exclude

            robs.append(r)
            rexp = theorypred.getRValue(expected=True)
            tpList.append( { "robs": r, "rexp": rexp, "tp": theorypred } )


        # while len(robs)<2:
        #     robs.append(0.)

        robs.sort ( reverse = True )
        tpList.sort ( key = lambda x: x['robs'], reverse = True )
        srs = ", ".join ( [ f"{x:.2f}" for x in robs[:3] ] )

        if keep_predictions:
            protomodel.ul_critic_robs = robs
            protomodel.ul_critic_tpList = tpList[:]

        critic_description = []
        for tp in tpList[:3]:
            rtype = tp['tp'].dataType(short=True)
            tmp = f"{tp['tp'].analysisId()}({rtype}):[robs={tp['robs']:.2f},rexp={tp['rexp']:.2f}]"
            critic_description.append ( tmp )
        if len(tpList)>3:
            critic_description.append ( "..." )
        protomodel.critic_description += "Datasets: " + ",".join ( critic_description )

        return


    def updateModelPredictionsWithCombinedPreds(self, protomodel, mostSensiComb, robsComb):
        """ Extract information from list of theory predictions and store r_obs from
            the most sensitive combination of analyses in the protomodel.

            :param mostSensiComb: list of theory predictions making the most sensitive combination
            :param robsComb: r_obs of the most sensitive combination
        """
        from ptools.helpers import experimentalId

        protomodel.llhd_critic_robs = robsComb

        if mostSensiComb is None:
            protomodel.description += "; llhd-based critic has no theory prediction."
        else:
            protomodel.critic_description += "; llhd-based critic combined datasets:" + ",".join( [experimentalId(comb) for comb in mostSensiComb] ) + f"with r={robsComb}"

        return


    def runSModelS(self, inputFile : PathLike, combineSRs : bool, ULpreds: bool, sigmacut : float, maxcond : float = 0.2 ) -> List[TheoryPrediction]:
        """ run smodels proper.
        :param inputFile: the input slha file
        :param ULpreds: if true, also returns the list of theory predictions for UL-type results
        :param sigmacut: the cut on the topology weights, typically 0.02*fb
        :param maxcond: maximum relative violation of conditions for valid results

        :returns: list of all theory predictions
        """

        if not os.path.exists ( inputFile ):
            self.highlight("error", f"cannot find inputFile {inputFile}" )
            return []
        model = Model ( BSMList, SMList )
        try:
            model.updateParticles ( inputFile=inputFile )
        except SModelSError as e:
            if "No cross-sections found" in str(e):
                # okay, everything under control, we just return empty list
                # of theory predictions
                return []
            else:
                # no idea what that is. pass it on.
                raise e

        mingap=10*GeV
        topos = decomposer.decompose ( model, sigmacut, minmassgap=mingap )
        if False:
            from smodels.base import runtime
            runtime._experimental = True
        combinedIds = set() # the analysis ids of the combined
        datasetPreds = [] # the SR specific predictions
        predictions = []
        ulpreds = []

        theoryPredictions = theoryPredictionsFor ( self.database, topos, useBestDataset=True, combinedResults=combineSRs )
        preds = TheoryPredictionList(theoryPredictions, maxcond)

        if preds != None:
            for pred in preds:
                if pred.dataType() == 'upperLimit':
                    ulpreds.append ( pred )
                    continue
                datasetPreds.append ( pred )

        for pred in datasetPreds:
            if pred.dataType() == "combined":
                predictions.append ( pred )
                combinedIds.add ( pred.dataset.globalInfo.id )
        for pred in datasetPreds:
            if pred.dataset.globalInfo.id in combinedIds:
                continue
            predictions.append ( pred )
        if ULpreds:
            return ulpreds, predictions

        return predictions

    def predict_critic(self, protomodel : ProtoModel, sigmacut = 0.02*fb, keep_predictions : bool = True, keep_slhafile : bool = False ):
        """ Compute the critic predictions and statistical variables, for a protomodel.

        :param sigmacut: weight cut on the predict xsecs for theoryPredictions
        :param keep_predictions: if True, then keep *all* predictions -- not just the one that make it into the combination, store them as
                                 protomodel(self).predictions. Store the predictions for the critic in protomodel(self).critic_preds.
        :param keep_slhafile: if True, then keep the temporary slha file, print out its name

        :returns: False, if no combinations could be found, else True
        """
        # Create SLHA file (for running SModelS)
        slhafile = protomodel.createSLHAFile()

        # --- UL-based critic ---

        # Run SModelS to get for UL-type predictions, and best SR preditcions if no UL-type result.
        UL_preds, bestSR_preds = self.runSModelS( slhafile, combineSRs=False, ULpreds=True, sigmacut=sigmacut)

        # Use best SR preds only if no UL-type result.
        predictions = self.merge_preds(UL_preds,bestSR_preds)

        allowed_by_UL_critic = self.UL_critic(protomodel, predictions)

        # Extract the relevant prediction information and store in the protomodel:
        self.updateModelPredictionsWithULPreds(protomodel, predictions, keep_predictions)

        if not allowed_by_UL_critic:
            if keep_slhafile:
                self.info( f"Keeping {protomodel.currentSLHA}, as requested" )
            else:
                protomodel.delCurrentSLHA()
            return False

        self.info("Model allowed by UL-based critic. Starting llhd-based critic.")

        # --- llhd-based critic ---

        predictions = self.runSModelS( slhafile, combineSRs=True, ULpreds=False, sigmacut=sigmacut )

        allowed_by_llhd_critic, mostSensiComb, robsComb = self.llhd_critic(predictions, cut=0.1, keep_predictions=keep_predictions)

        # Extract the relevant prediction information and store in the protomodel:
        self.updateModelPredictionsWithCombinedPreds(protomodel, mostSensiComb, robsComb)

        if keep_slhafile:
            self.info(f"Keeping {protomodel.currentSLHA}, as requested" )
        else:
            protomodel.delCurrentSLHA()

        # Compute the maximum allowed (global) mu value given the r-values
        # stored in protomodel
        # protomodel.mumax = self.getMaxAllowedMu(protomodel)

        if allowed_by_llhd_critic:
            self.log(f"Model passed llhd-based critic with critic robs = {robsComb}.")
            return True
        else:
            self.info(f"Model failed llhd-based critic with critic robs = {robsComb}.")
            return False


    def merge_preds(self, pred_list_1, pred_list_2):
        """ Small script to merge the preds of pred_list_2 into pred_list_1 but
        only if no pred of pred_list_1 shares the same analysis ID """

        predictions = []
        set_ids = set ()

        for pred in pred_list_1:    #Ul preds
            predictions.append ( pred )
            set_ids.add ( pred.dataset.globalInfo.id )
        for pred in pred_list_2:    #best SR preds
            id = pred.dataset.globalInfo.id
            for ext in ['agg','ma5','ewk','strong','hino','multibin','exclusive','incl','adl','eff']:
                id = id.replace(f'-{ext}','')
            if id in set_ids:
                continue
            predictions.append ( pred )

        return predictions


    def UL_critic(self, protomodel, predictions):
        """ UL-based critic (can also use best SR results if no UL-type result available for a given analysis).

        :param predictions: list of theory predictions (UL-type and EM-type)

        :returns: False if the critic excludes the model, else True.
        """

        if not predictions: # If empty list
            return True     # the model is not excluded

        from scipy.stats import binom

        n_sensitive, n_excluding = 0, 0

        for pred in predictions:
            robs = pred.getRValue (expected=False)
            rexp = pred.getRValue (expected=True)

            if rexp is None:
                rexp = robs
            if rexp < self.sensitivity_threshold:
                continue

            n_sensitive += 1
            if robs > self.r_threshold:
                n_excluding += 1

        max_allowed = 0
        while binom.cdf(max_allowed,n_sensitive,0.05) <= 0.66:
            max_allowed += 1

        protomodel.critic_description = f"UL-based critic: n_sensitive={n_sensitive}, n_excluding={n_excluding}, max_allowed={max_allowed} => passes critic: {max_allowed >= n_excluding}. "
        self.log(f"UL-based critic: n_sensitive={n_sensitive}, n_excluding={n_excluding}, max_allowed={max_allowed} => passes critic: {max_allowed >= n_excluding}")
        return max_allowed >= n_excluding


    def llhd_critic(self, predictions, cut=0, keep_predictions=False):
        """ llhd-based critic.

        :param predictions: list of theory predictions (EM-type only). Combined dataset when available, best SR otherwise.
        :param cut: theory predictions giving an r_exp below this cut will not enter the combination

        :returns: False if the critic excludes the model, else True.
        """

        if not predictions:          # If empty list
            return True, None, None  # the model is not excluded

        EMpreds = []

        if cut > 0:
            for tpred in predictions:
                rexp = tpred.getRValue(expected = True)
                if rexp is not None and rexp >= cut:
                    EMpreds.append(tpred)
        else:
            EMpreds = predictions

        if keep_predictions:
            self.llhd_critic_preds = EMpreds
        self.log( f"Found {len(EMpreds)} llhd-based critic predictions." )
        r = None
        best_comb, _ = self.combiner.getMostSensitiveCombination(predictions)
        if best_comb:
            tpCombiner = TheoryPredictionsCombiner(best_comb)
            r = tpCombiner.getRValue(expected=False)

        if r is None:
            self.highlight("warning","The computation of the observed r-value of the most sensitive combination gave None.")
            return False, best_comb, None

        return r < 1, best_comb, r          #SN: r < r_threshold?
