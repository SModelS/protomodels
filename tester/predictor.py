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
from builder.manipulator import Manipulator
from smodels.share.models.SMparticles import SMList
from share.model_spec import BSMList
from smodels.base.physicsUnits import fb, GeV, TeV
from smodels.experiment.databaseObj import Database
from smodels.base.model import Model
from smodels.base.exceptions import SModelSBaseError as SModelSError
from os import PathLike
from typing import List, Union
from base.loggerbase import LoggerBase
from tester.combinationsmatrix import getYamlMatrix

try:
    from tester.combiner import Combiner
except:
    from combiner import Combiner

class Predictor ( LoggerBase ):
    def __init__ ( self, walkerid : int, dbpath : PathLike = "official",
                   expected : bool = False, select : str = "all",
                   do_srcombine : bool = True ):
        """
        the predictor class, i.e. the class that computes the predictions,
        finds the best combinations, and computes the final test statistics

        :param do_srcombine: if True, then also use combined results,
                           both via simplified likelihoods and pyhf.
        """
        super ( Predictor, self ).__init__ ( walkerid )
        self.walkerid = walkerid
        self.do_srcombine = do_srcombine
        self.modifier = None
        self.select = select
        self.expected = expected

        if expected:
            from expResModifier import ExpResModifier
            self.modifier = ExpResModifier()
        force_load = None
        if dbpath.endswith ( ".pcl" ):
            force_load = "pcl"
        if "/" in dbpath:
            ntries = 0
            while not os.path.exists ( dbpath ):
                ## give it a few tries
                ntries += 1
                time.sleep ( ntries * 5 )
                if ntries > 5:
                    break

        combinationsmatrix, status = getYamlMatrix()
        if not combinationsmatrix or status != 0:
            sys.exit("Combination matrix not loaded correctly when instantiating Predictor class.")

        self.database=Database( dbpath, force_load = force_load, combinationsmatrix = combinationsmatrix )
        self.fetchResults()
        self.combiner = Combiner(self.walkerid)

    def filterForAnaIdsTopos ( self, anaIds, topo ):
        """ filter the list of expRes, keep only anaIds """
        keepExpRes = []
        nbefore = len(self.listOfExpRes)
        for er in self.listOfExpRes:
            eid = er.globalInfo.id
            if not eid in anaIds:
                continue
            txnames = [ x.txName for x in er.getTxNames() ]
            if not topo in txnames: ## can safely skip
                continue
            newDS = []
            for dataset in er.datasets:
                newTxNames = []
                for txName in dataset.txnameList:
                    if txName.txName != topo:
                        continue
                    newTxNames.append ( txName )
                if len(newTxNames)>0:
                    dataset.txnameList = newTxNames
                    newDS.append ( dataset )
            if len(newDS)>0:
                er.datasets = newDS
                keepExpRes.append ( er )
        nafter=len(keepExpRes)
        self.pprint ( f"filtered for {topo}, keeping {nafter}/{nbefore} expRes" )
        self.listOfExpRes = keepExpRes

    def matchesTopos ( self, topo, listOfTopos ):
        """ tell me if topo is in listOfTopos """
        if topo in listOfTopos:
            return True
        return False

    def filterForTopos ( self, topos : Union[List,str] ):
        """ filter the list of expRes, keep only the ones for topo
        :param topos: the topologies, either as list or as string with commas
        """
        keepExpRes = []
        nbefore = len(self.listOfExpRes)
        if type(topos) == str:
            from ptools import moreHelpers
            topos = topos.split(",")
        actualtopos = set()
        for topo in topos:
            names = moreHelpers.namesForSetsOfTopologies(topo)
            ts = names[0].split(",")
            for t in ts:
                actualtopos.add ( t )
        for er in self.listOfExpRes:
            txnames = [ x.txName for x in er.getTxNames() ]
            anyMatch = False
            for t in actualtopos:
                if self.matchesTopos ( t, txnames ):
                    anyMatch = True
            if not anyMatch: ## can safely skip
                continue
            newDS = []
            for dataset in er.datasets:
                newTxNames = []
                for txName in dataset.txnameList:
                    if not self.matchesTopos ( txName.txName, actualtopos ):
                        continue
                    newTxNames.append ( txName )
                if len(newTxNames)>0:
                    dataset.txnameList = newTxNames
                    newDS.append ( dataset )
            if len(newDS)>0:
                er.datasets = newDS
                keepExpRes.append ( er )
        self.pprint ( f"filtered for {','.join(actualtopos)}, keeping {len(keepExpRes)}/{nbefore} expRes" )
        self.listOfExpRes = keepExpRes

    def filterForSqrts ( self, keep : List = [ 13 ] ):
        """ filter listOfExpRes for sqrts.

        :param keep: list of sqrtses to *keep*
        """
        keepExpRes = []
        for er in self.listOfExpRes:
            sqrts = er.globalInfo.sqrts.asNumber(TeV)
            for k in keep:
                ds = abs ( sqrts - k )
                if ds < 1e-5:
                    keepExpRes.append ( er )
        self.listOfExpRes = keepExpRes

    def fetchResults ( self ):
        """ fetch the list of results, perform all selecting
            and modding """

        dataTypes = [ "all" ]
        if self.select == "em":
            dataTypes = [ "efficiencyMap" ]
        if self.select == "ul":
            dataTypes = [ "upperLimit" ]
        txnames = [ "all" ]
        if self.select.startswith("txnames:"):
            s = self.select.replace("txnames:","")
            from ptools.moreHelpers import namesForSetsOfTopologies
            txnames = namesForSetsOfTopologies ( s )[0]
            self.log ( f"I have been asked to select txnames for {txnames}" )
            txnames = txnames.split(",")

        listOfExpRes = self.database.getExpResults( dataTypes = dataTypes,
                                                    txnames = txnames,
                                                    useNonValidated=True )
        if self.modifier:
            listOfExpRes = self.modifier.modify ( listOfExpRes )
        self.log ( f"I will be working with {len(listOfExpRes)} results" )

        self.listOfExpRes = listOfExpRes
        if False:
            f=open("expresults.txt","wt")
            for expRes in self.listOfExpRes:
                f.write ( f"{expRes.id()} {expRes.datasets[0].dataInfo.dataId}\n" )
            f.close()

    # def removeRedundantULResults ( self,
    #         predictions : List[TheoryPrediction] ) -> List[TheoryPrediction]:
    #     """ from the given predictions, return UL results, if there is
    #     also a combined result """
    #     ret = []
    #     hasCombined = set()
    #     for p in predictions:
    #         if str(p.dataId()) == "(combined)":
    #             hasCombined.add ( p.analysisId() )
    #     for p in predictions:
    #         if p.dataType() == "upperLimit":
    #             if p.analysisId() in hasCombined:
    #                 continue
    #         ret.append ( p )
    #     return ret

    def predict ( self, protomodel : ProtoModel, sigmacut = 0.02*fb,
                  strategy : str = "aggressive",keep_predictions : bool = False, keep_slhafile : bool = False ) -> bool:
        """ Compute the predictions and statistical variables, for a
            protomodel.

        :param sigmacut: weight cut on the predict xsecs for theoryPredictions
        :param strategy: combination strategy, currently only aggressive is used
        :param keep_predictions: if True, then keep *all* predictions --
        not just the one that make it into the combination, store them as
        predictor(self).predictions. Store the predictions for the critic in
        predictor(self).critic_preds.
        :param keep_slhafile: if True, then keep the temporary slha file,
        print out its name
        :returns: False, if no combinations could be found, else True
        """


        if hasattr ( self, "predictions" ):
            del self.predictions ## make sure we dont accidentally use old preds
        self.walkerid = protomodel.walkerid ## set the walker ids, for debugging
        self.combiner.walkerid = protomodel.walkerid

        # Create SLHA file (for running SModelS)
        slhafile = protomodel.createSLHAFile()

        # now use all prediction with likelihood values to compute the TL of the model
        predictions = self.runSModelS( slhafile, sigmacut, allpreds=True, ULpreds=False )
        if not predictions: return False

        if keep_predictions:
            self.predictions = predictions

        # Compute significance and store in the model:
        self.computeSignificance( protomodel, predictions, strategy )

        if protomodel.TL is None:
            self.log ( f"done with prediction. Could not find combinations (TL={protomodel.TL})" )
            protomodel.delCurrentSLHA()
            return False
        else:
            self.log ( f"done with prediction. best TL={protomodel.TL:.2f} (muhat={protomodel.muhat:.2f})" )

        protomodel.cleanBestCombo()

        #Recompute predictions with higher accuracy for high score models:
        ## FIXME not a good idea! should probably remove!
        if False and protomodel.TL > 4.1 and protomodel.nevents < 55000:
            self.log ( f"TL {protomodel.TL:.2f}>2.7, repeat with higher stats!" )
            protomodel.nevents = 100000
            protomodel.computeXSecs()
            self.predict(protomodel,sigmacut=sigmacut, strategy= strategy,
                    keep_predictions = keep_predictions )

        if keep_slhafile:
            self.log ( f"keeping {protomodel.currentSLHA}, as requested" )
        else:
            protomodel.delCurrentSLHA()
        # we keep track of the database version, when predicting
        protomodel.dbversion = self.database.databaseVersion
        return True

    def runSModelS(self, inputFile : PathLike, sigmacut : float,
            allpreds : bool, ULpreds : bool, maxcond : float = 0.2 ) -> List[TheoryPrediction]:
        """ run smodels proper.
        :param inputFile: the input slha file
        :param sigmacut: the cut on the topology weights, typically 0.02*fb
        :param allpreds: if true, return all predictions of analyses, else
                         only best signal region
        :param ULpreds: if true, also returns the list of theory predictions for UL-type results
        :param maxcond: maximum relative violation of conditions for valid results

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
                self.log(f"No cross-sections found for protomodel point")
                return []
            else:
                # no idea what that is. pass it on.
                raise e

        mingap=10*GeV

        # self.log ( "Now decomposing" )
        topos = decomposer.decompose ( model, sigmacut, minmassgap=mingap )
        self.log ( f"decomposed model into {len(topos)} topologies." )

        if allpreds:
            bestDataSet=False
        else:
            bestDataSet=True

        # self.log ( "start getting preds" )
        if False:
            from smodels.base import runtime
            runtime._experimental = True
        combinedIds = set() # the analysis ids of the combined
        EMpreds = [] # the SR specific predictions
        predictions = []
        ulpreds = []
        self.log("start computing preds")
        import time
        start_time = time.time()
        theoryPredictions = theoryPredictionsFor ( self.database, topos, useBestDataset=bestDataSet, combinedResults=self.do_srcombine )
        preds = TheoryPredictionList(theoryPredictions, maxcond)

        end_time = time.time()
        time_taken = end_time - start_time
        self.log(f"finished computing preds, taken {time_taken:.6f} seconds")
        if preds != None:
            for pred in preds:
                if pred.dataType() == 'upperLimit':
                    ulpreds.append ( pred )
                    continue
                EMpreds.append ( pred )

        for pred in EMpreds:
            if pred.dataType() == "combined":
                predictions.append ( pred )
                combinedIds.add ( pred.dataset.globalInfo.id )
        for pred in EMpreds:
            if pred.dataset.globalInfo.id in combinedIds:
                continue
            predictions.append ( pred )

        sap = "best preds"
        if allpreds:
            sap = "all preds"
        sllhd = ""
        if ULpreds:
            sllhd = ", with UL preds"
        self.log( f"returning {len(predictions)} predictions, {sap}{sllhd}" )

        if ULpreds:
            return ulpreds, predictions
        else:
            return predictions

    def printPredictions ( self ):
        """ if self.predictions exists, pretty print them """
        if not hasattr ( self, "predictions" ):
            print ( "[predictor] no predictions. did you run .predict( ..., keep_predictions=True )?" )
        if hasattr ( self, "predictions" ):
            print ( f"[predictor] {len(self.predictions)} predictions for combiner:" )
            for p in self.predictions:
                dataId = "combined"
                if hasattr ( p.dataset, "dataInfo" ):
                    dataId = p.dataset.dataInfo.dataId
                if dataId == None:
                    dataId = "UL"
                txns = ",".join ( set ( map ( str, p.txnames ) ) )
                print ( f" - {p.analysisId()}:{dataId}: {txns}" )
        if hasattr ( self, "critic_preds" ):
            print ( )
            print ( f"[predictor] {len(self.critic_preds)} predictions for critic:" )
            for p in self.critic_preds:
                dataId = "combined"
                if hasattr ( p.dataset, "dataInfo" ):
                    dataId = p.dataset.dataInfo.dataId
                if dataId == None:
                    dataId = "UL"
                txns = ",".join ( set ( map ( str, p.txnames ) ) )
                print ( f" - {p.analysisId()}:{dataId}: {txns}" )


    def computeSignificance(self, protomodel, predictions, strategy):
        """ compute the K and TL values, and attach them to the protomodel """
        if len ( predictions ) == 0:
            protomodel.K = None
            protomodel.TL = None
            return

        self.log ( f"now find combo with highest TL given {len(predictions)} predictions" )
        ## find highest observed significance
        #(set mumax just slightly below its value, so muhat is always below)
        # mumax = protomodel.mumax
        bestCombo,TL,muhat = self.combiner.findHighestSignificance ( predictions, expected=False )

        ## DONT normalize here! rescaling by muhat is done in randomwalker predict function!
        '''
        ma = Manipulator ( protomodel )
        ma.rescaleSignalBy ( muhat )
        muhat = 1.
        '''

        if hasattr ( protomodel, "keep_meta" ) and protomodel.keep_meta:
            protomodel.bestCombo = bestCombo
        else:
            protomodel.bestCombo = self.combiner.removeDataFromBestCombo ( bestCombo )
        protomodel.letters = self.combiner.getLetterCode(protomodel.bestCombo)
        protomodel.description = self.combiner.getComboDescription(protomodel.bestCombo)

        protomodel.muhat = muhat

        protomodel.TL = TL

        if abs(muhat - 1.0) < 1e-02:
            prior = self.combiner.computePrior ( protomodel )
            ## temporary hack: penalize for missing experiment
            missingExpPenalty = self.combiner.penaltyForMissingResults ( predictions )
            extremeSSMs = self.combiner.penaltyForExtremeSSMs ( protomodel )
            undemocraticFlavors = self.combiner.penaltyForUndemocraticFlavors ( protomodel )
            oldprior = prior
            prior *= missingExpPenalty * extremeSSMs * undemocraticFlavors
            self.log ( f"prior={prior:.2f} before_penalties={oldprior:.2f} "\
                       f"missingExp={missingExpPenalty:.2f} "\
                       f"extremeSSMs={extremeSSMs:.2f} "\
                       f"undemocraticFlavors={undemocraticFlavors:.2f}" )

            if TL is None: # TL is None when no combination was found
                protomodel.K = None
            else:
                protomodel.K = self.combiner.computeK ( TL, prior )
            #protomodel.llhd = llhd


if __name__ == "__main__":
    import argparse
    argparser = argparse.ArgumentParser(
            description='calling the predictor from the command line' )
    argparser.add_argument ( '-d', '--database',
            help='Path to database pickle file [./default.pcl]',
            type=str, default="./default.pcl" )
    argparser.add_argument ( '-f', '--hiscore',
            help='input model, to be taken from hiscore file [./hiscores.cache]',
            type=str, default="./model.slha" )
    argparser.add_argument ( '-o', '--output',
            help='output pickle file [./predictions.pcl]',
            type=str, default="./predictions.pcl" )
    argparser.add_argument ( '-i', '--interactive',
            help='interactive shell',
            action="store_true" )
    args = argparser.parse_args()

    p = Predictor ( 0, args.database, do_srcombine=False )

    sys.path.insert(0,"../")
    from walker.hiscore import Hiscores
    hiscore = Hiscores ( 0, False, args.hiscore )
    hi = hiscore.hiscores[0]
    print ( f"Will scrutinize hiscore obtained with database {hi.dbversion} K={hi.K:.3f}" )
    oldK = hi.K

    hi.K=-3 # make sure it gets recomputed
    print ( "to be sure i reset K to %.2f" % hi.K )

    predictions = p.predict ( hi )

    if args.output not in [ "", "none", "None" ]:
        with open( args.output, "wb" ) as f:
            pickle.dump ( predictions, f )

    if args.interactive:
        import IPython
        IPython.embed ( using = False )
