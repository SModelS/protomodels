#!/usr/bin/env python3

"""
the predictor class, i.e. the class that computes the predictions,
finds the best combinations, and computes the final test statistics
"""

__all__ = [ "Predictor" ]

import pickle, time, os, sys
from smodels.decomposition import decomposer
from smodels.matching.theoryPrediction import theoryPredictionsFor, TheoryPrediction, TheoryPredictionList, TheoryPredictionsCombiner
sys.path.insert(0,"../")
from builder.protomodel import ProtoModel
from smodels.share.models.SMparticles import SMList
from smodels.share.models.mssm import BSMList
from smodels.base.physicsUnits import fb, GeV, TeV
from smodels.experiment.databaseObj import Database
from smodels.base.model import Model
from smodels.base.exceptions import SModelSBaseError as SModelSError
from os import PathLike
from typing import List, Union
from builder.loggerbase import LoggerBase
from tester.combiner import Combiner

class Critic ( LoggerBase ):
    def __init__ ( self, walkerid : int, dbpath : PathLike = "official", expected : bool = False, select : str = "all", do_srcombine : bool = False ):
        self.walkerid = walkerid
        self.do_srcombine = do_srcombine
        self.rthreshold = 1.38
        
        force_load = None
        if dbpath.endswith ( ".pcl" ):
            force_load = "pcl"
            
        self.database=Database( dbpath, force_load = force_load )
        self.combiner = Combiner(self.walkerid)
    
    def getMaxAllowedMu(self, protomodel):
        """ Compute the maximum (global) signal strength normalization
            given the predictions.
        """

        mumax = float("inf")
        if protomodel.ul_type_rvalues[0] > 0.:
            #Set mumax slightly below threshold, so the model is never excluded
            print(f"r value from UL: {protomodel.ul_type_rvalues[0]}")
            mumax_ul = 0.999*self.rthreshold / protomodel.ul_type_rvalues[0]
            #if 0 < protomodel.llhd_type_rvalue < mumax:
            print(f"r value from llhd: {protomodel.llhd_type_rvalue}")
            mumax_llhd = 0.999*self.rthreshold /protomodel.llhd_type_rvalue
            mumax = min(mumax_ul, mumax_llhd)

        return mumax
    
    def updateModelPredictionsWithULPreds(self, protomodel, ul_critic_preds):
        """ Extract information from list of theory predictions and store list of dict with r_obs,
            r_exp and theory prediction(sorted according to decreasing r_obs values) in the protomodel.
            Also store description about the critic_tp in the protomodel.

            :param ul_critic_preds: theory predictions for UL-type results
        """

        rvalues = [] #If there are no predictions set rmax and r2 to 0
        tpList = []
        for theorypred in ul_critic_preds:
            r = theorypred.getRValue(expected=False)
            if r == None:   #SN: !Fixme Error Message if r is None
                self.pprint ( "I received %s as r. What do I do with this?" % r )
                r = 23. #TP: Why high r and not 0? SN: if r is None, we have done sth wrong, so we dont want to condsider this protomodel
            rexp = theorypred.getRValue(expected=True)
            tpList.append( { "robs": r, "rexp": rexp, "tp": theorypred } )
            rvalues.append(r)
        while len(rvalues)<2:
            rvalues.append(0.)
        rvalues.sort(reverse = True )
        tpList.sort ( key = lambda x: x['robs'], reverse = True )
        srs = ", ".join ( [ f"{x:.2f}" for x in rvalues[:3] ] )
        #self.log ( f"top r values before rescaling are: {srs}" )
        protomodel.ul_type_rvalues = rvalues #Do not include initial zero values
        # protomodel.excluded = protomodel.rvalues[0] > self.rthreshold #The 0.99 deals with the case rmax = threshold
        protomodel.ul_type_tpList = tpList[:]
        critic_description = []
        for tp in tpList[:3]:
            rtype = tp['tp'].dataType(short=True)
            tmp = f"{tp['tp'].analysisId()}({rtype}):{tp['robs']:.2f}"
            critic_description.append ( tmp )
        if len(tpList)>3:
            critic_description.append ( "...")
        protomodel.critic_description = "Fast critic datasets:" + ",".join ( critic_description )

    def updateModelPredictionsWithCombinedPreds(self, protomodel, llhd_critic_preds, cut, keep_predictions = False):
        """ Extract information from list of theory predictions and store robs from
            the most sensitive combination of analyses in the protomodel.

            :param llhd_critic_preds: theory predictions for EM-type results
            :param cut: theory predictions giving an r_exp below this cut will not enter the combination
        """
        from ptools.helpers import experimentalId

        EMpreds = []

        if cut > 0:
            for tpred in llhd_critic_preds:
                r_exp = tpred.getRValue(expected = True)
                if r_exp >= cut:
                    EMpreds.append(tpred)
            if keep_predictions:
                self.llhd_critic_preds = EMpreds
        else:
            EMpreds = llhd_critic_preds

        best_comb, weight = self.combiner.getMostSensitiveCombination(llhd_critic_preds)

        if best_comb:
            tpCombiner = TheoryPredictionsCombiner(best_comb)
            r = tpCombiner.getRValue(expected=False)

        if r is None:
            logger.warning ("The computation of the observed r-value of the most sensitive combination gave None.")

        protomodel.llhd_type_rvalue = r
        protomodel.critic_description += "; llhd-based critic combined datasets:" + ",".join( [experimentalId(comb) for comb in best_comb] ) + f"with r={r}"

        return
            
    def runSModelS(self, inputFile : PathLike, sigmacut, allpreds : bool, ULpreds : bool ) -> List[TheoryPrediction]:
        """ run smodels proper.
        :param inputFile: the input slha file
        :param sigmacut: the cut on the topology weights, typically 0.02*fb
        :param allpreds: if true, return all predictions of analyses, else only best signal region
        :param llhdonly: if true, return only results with likelihoods

        :returns: list of all theory predictions
        """

        if not os.path.exists ( inputFile ):
            self.pprint ( f"error, cannot find inputFile {inputFile}" )
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

        # self.log ( "Now decomposing" )
        topos = decomposer.decompose ( model, sigmacut, minmassgap=mingap )
        #self.log ( f"decomposed model into {len(topos)} topologies." )


        if allpreds:
            bestDataSet=False
        else:
            bestDataSet=True

        # self.log ( "start getting preds" )
        if False:
            from smodels.base import runtime
            runtime._experimental = True
        combinedIds = set() # the analysis ids of the combined
        datasetPreds = [] # the SR specific predictions
        predictions = []
        ulpreds = []

        preds = theoryPredictionsFor ( self.database, topos, useBestDataset=bestDataSet, combinedResults=self.do_srcombine )
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

        sap = "best preds"
        if allpreds:
            sap = "all preds"
        sllhd = ""
        if ULpreds:
            sllhd = ", with UL preds"
        #self.log( f"returning {len(predictions)} predictions, {sap}{sllhd}" )

        if ULpreds:
            return ulpreds, predictions
        else:
            return predictions
    
    def predict_critic(self, protomodel : ProtoModel, sigmacut = 0.02*fb, keep_predictions : bool = False, keep_slhafile : bool = False ):
        """ Compute the crtic predictions and statistical variables, for a protomodel.

        :param sigmacut: weight cut on the predict xsecs for theoryPredictions
        :param keep_predictions: if True, then keep *all* predictions -- not just the one that make it into the combination, store them as
                                 predictor(self).predictions. Store the predictions for the critic in predictor(self).critic_preds.
        :param keep_slhafile: if True, then keep the temporary slha file, print out its name
        :returns: False, if no combinations could be found, else True
        """
        # Create SLHA file (for running SModelS)
        slhafile = protomodel.createSLHAFile()
            # First run SModelS using all results and considering only the best signal region.
        # thats the run for the critic
        ul_critic_preds, llhd_critic_preds = self.runSModelS( slhafile, sigmacut, allpreds=False, ULpreds=True )

        if keep_predictions:
            self.ul_critic_preds = ul_critic_preds
            self.llhd_critic_preds = llhd_critic_preds
        # Extract the relevant prediction information and store in the protomodel:
        self.updateModelPredictionsWithULPreds(protomodel,ul_critic_preds)
        self.updateModelPredictionsWithCombinedPreds(protomodel,llhd_critic_preds, cut=0.1, keep_predictions=keep_predictions)
        # self.log ( "model is excluded? %s" % str(protomodel.excluded) )

        # Compute the maximum allowed (global) mu value given the r-values
        # stored in protomodel
        protomodel.mumax = self.getMaxAllowedMu(protomodel)
