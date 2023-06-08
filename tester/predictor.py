#!/usr/bin/env python3

""" store the theory predictions in pickle """

__all__ = [ "Predictor" ]

import pickle, time, os, sys
from smodels.theory import decomposer
from smodels.theory.theoryPrediction import theoryPredictionsFor, TheoryPrediction
from protomodels.builder.protomodel import ProtoModel
from smodels.share.models.SMparticles import SMList
from smodels.particlesLoader import BSMList
from smodels.tools.physicsUnits import fb, GeV
from smodels.experiment.databaseObj import Database
from smodels.theory.model import Model
from os import PathLike
from typing import List

try:
    from tester.combiner import Combiner
except:
    from combiner import Combiner

class Predictor:
    def __init__ ( self, walkerid : int, dbpath : PathLike = "./default.pcl",
                   expected : bool = False, select : str = "all",
                   do_combine : bool = False ):
        """
        :param do_combine: if True, then also use combined results,
                           both via simplified likelihoods and pyhf.
        """
        self.walkerid = walkerid
        self.do_combine = do_combine
        self.modifier = None
        self.select = select
        self.expected = expected
        self.rthreshold = 1.3 ## threshold for rmax
        if expected:
            from expResModifier import ExpResModifier
            self.modifier = ExpResModifier()
        force_load = None
        if dbpath.endswith ( ".pcl" ):
            force_load = "pcl"
        ntries = 0
        while not os.path.exists ( dbpath ):
            ## give it a few tries
            ntries += 1
            time.sleep ( ntries * 5 )
            if ntries > 5:
                break
        self.database=Database( dbpath, force_load = force_load )
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
        self.pprint ( "filtered for %s, keeping %d/%d expRes" % \
                      ( topo, len(keepExpRes), nbefore) )
        self.listOfExpRes = keepExpRes

    def filterForTopos ( self, topo ):
        """ filter the list of expRes, keep only the ones for topo """
        keepExpRes = []
        nbefore = len(self.listOfExpRes)
        for er in self.listOfExpRes:
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
        self.pprint ( "filtered for %s, keeping %d/%d expRes" % \
                      ( topo, len(keepExpRes), nbefore) )
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
            txnames = s.split(",")
            self.pprint ( "I have been asked to select txnames for %s" % s )

        listOfExpRes = self.database.getExpResults( dataTypes = dataTypes,
                                                    txnames = txnames,
                                                    useNonValidated=True )
        if self.modifier:
            listOfExpRes = self.modifier.modify ( listOfExpRes )

        self.listOfExpRes = listOfExpRes
        if False:
            f=open("expresults.txt","wt")
            for expRes in self.listOfExpRes:
                f.write ( "%s %s\n" % (expRes.id(), expRes.datasets[0] ) )
            f.close()

    def pprint ( self, *args ):
        """ logging """
        print ( "[predictor] %s" % ( " ".join(map(str,args))) )
        self.log ( *args )

    def log ( self, *args ):
        """ logging to file """
        with open( "walker%d.log" % self.walkerid, "a" ) as f:
            f.write ( "[predictor-%s] %s\n" % ( time.strftime("%H:%M:%S"), " ".join(map(str,args)) ) )

    def predict ( self, protomodel : ProtoModel, sigmacut = 0.02*fb,
                  strategy : str = "aggressive", keep_predictions : bool = False ):
        """ Compute the predictions and statistical variables, for a
            protomodel.

        :param sigmacut: weight cut on the predict xsecs for theoryPredictions
        :param strategy: combination strategy, currently only aggressive is used
        :param keep_predictions: if True, then keep all predictions (in self,
               not in protomodel!!)
        :returns: False, if no combinations could be found, else True
        """

        if hasattr ( self, "predictions" ):
            del self.predictions ## make sure we dont accidentally use old preds

        # Create SLHA file (for running SModelS)
        slhafile = protomodel.createSLHAFile()

        # First run SModelS using all results and considering only the best signal region.
        # thats the run for the critic
        critic_preds = self.runSModelS( slhafile, sigmacut,  allpreds=False,
                                           llhdonly=False )

        if keep_predictions:
            self.critic_preds = critic_preds
        # Extract the relevant prediction information and store in the protomodel:
        self.updateModelPredictions(protomodel,critic_preds)
        # self.log ( "model is excluded? %s" % str(protomodel.excluded) )

        # Compute the maximum allowed (global) mu value given the r-values 
        # stored in protomodel
        protomodel.mumax = self.getMaxAllowedMu(protomodel)

        # now use all prediction with likelihood values to compute the Z of the model
        predictions = self.runSModelS( slhafile, sigmacut, allpreds=True,
                                               llhdonly=True )

        if keep_predictions:
            self.predictions = predictions

        # Compute significance and store in the model:
        self.computeSignificance( protomodel, predictions, strategy )
        if protomodel.Z is None:
            self.log ( f"done with prediction. Could not find combinations (Z={protomodel.Z})" )
            protomodel.delCurrentSLHA()
            return False
        else:
            self.log ( f"done with prediction. best Z={protomodel.Z:.2f} (muhat={protomodel.muhat:.2f})" )

        protomodel.cleanBestCombo()

        #Recompute predictions with higher accuracy for high score models:
        if protomodel.Z > 4.1 and protomodel.nevents < 55000:
            self.log ( f"Z {protomodel.Z:.2f}>2.7, repeat with higher stats!" )
            protomodel.nevents = 100000
            protomodel.computeXSecs()
            self.predict(protomodel,sigmacut=sigmacut, strategy= strategy,
                    keep_predictions = keep_predictions )

        protomodel.delCurrentSLHA()
        return True

    def runSModelS(self, inputFile : PathLike, sigmacut, 
            allpreds : bool, llhdonly : bool ) -> List[TheoryPrediction]:
        """ run smodels proper.
        :param inputFile: the input slha file
        :param sigmacut: the cut on the topology weights, typically 0.02*fb
        :param allpreds: if true, return all predictions of analyses, else
                         only best signal region
        :param llhdonly: if true, return only results with likelihoods
        """

        if not os.path.exists ( inputFile ):
            self.pprint ( f"error, cannot find inputFile {inputFile}" )
            return []
        model = Model ( BSMList, SMList )
        model.updateParticles ( inputFile=inputFile )

        mingap=10*GeV

        # self.log ( "Now decomposing" )
        topos = decomposer.decompose ( model, sigmacut, minmassgap=mingap )
        self.log ( f"decomposed model into {len(topos)} topologies." )
            

        if allpreds:
            bestDataSet=False
            combinedRes=False
        else:
            bestDataSet=True
            combinedRes=self.do_combine

        # self.log ( "start getting preds" )
        from smodels.tools import runtime
        runtime._experimental = True
        combinedIds = set() # the analysis ids of the combined
        srpreds = [] # the SR specific predictions
        predictions = []
        for expRes in self.listOfExpRes:
            # get the SR specific predictions
            srpred = theoryPredictionsFor ( expRes, topos,
                                           useBestDataset=bestDataSet,
                                           combinedResults=combinedRes )
            if srpred != None:
                for p in srpred:
                    srpreds.append ( p )
            if allpreds:
                # get the SR-combined predictions
                cpreds = theoryPredictionsFor ( expRes, topos,
                                                useBestDataset=False,
                                                combinedResults=self.do_combine )
                if cpreds != None:
                    for c in cpreds:
                        anaId = c.dataset.globalInfo.id
                        combinedIds.add ( anaId ) # add with and without -agg
                        anaId = anaId.replace("-agg","")
                        combinedIds.add ( anaId )
                    for c in cpreds:
                        anaId = c.dataset.globalInfo.id
                        # if the aggregated version is in, then take this out
                        if anaId+"-agg" in combinedIds:
                            continue
                        # definitely add the combined predictions
                        c.computeStatistics()
                        predictions.append ( c )
        for srpred in srpreds:
            srId = srpred.dataset.globalInfo.id
            srId = srId.replace("-agg","")
            # add the others if not a combined result exists, and if we have llhds
            # (or we are just not asking specifically for llhds)
            if srId in combinedIds: # we continue only with the combination
                continue
            if (not llhdonly) or (srpred.likelihood != None):
                srpred.computeStatistics()
                predictions.append ( srpred )
        sap = "best preds"
        if allpreds:
            sap = "all preds"
        sllhd = ""
        if llhdonly:
            sllhd = ", llhds only"
        self.log( f"returning {len(predictions)} predictions, {sap}{sllhd}" )
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

    def updateModelPredictions(self, protomodel, predictions):
        """ Extract information from list of theory predictions and store in the protomodel.
        :param predictions: all theory predictions
        :returns: list of tuples with observed r values, r expected and
                  theory prediction info (sorted with highest r-value first)
        """

        rvalues = [] #If there are no predictions set rmax and r2 to 0
        tpList = []
        for theorypred in predictions:
            r = theorypred.getRValue(expected=False)
            if r == None:
                self.pprint ( "I received %s as r. What do I do with this?" % r )
                r = 23.
            rexp = theorypred.getRValue(expected=True)
            # tpList.append( (r, rexp, self.combiner.removeDataFromTheoryPred ( theorypred ) ) )
            tpList.append( (r, rexp, theorypred ) )
            rvalues.append(r)
        while len(rvalues)<2:
            rvalues.append(0.)
        rvalues.sort(reverse = True )
        srs = "%s" % ", ".join ( [ "%.2f" % x for x in rvalues[:3] ] )
        self.log ( "top r values before rescaling are: %s" % srs )
        protomodel.rvalues = rvalues #Do not include initial zero values
        # protomodel.excluded = protomodel.rvalues[0] > self.rthreshold #The 0.99 deals with the case rmax = threshold
        protomodel.tpList = tpList[:]

    def getMaxAllowedMu(self, protomodel):
        """ Compute the maximum (global) signal strength normalization
            given the predictions.
        """

        mumax = float("inf")
        if protomodel.rvalues[0] > 0.:
            #Set mumax slightly below threshold, so the model is never excluded
            mumax = 0.999*self.rthreshold / protomodel.rvalues[0]

        return mumax

    def computeSignificance(self, protomodel, predictions, strategy):
        """ compute the K and Z values, and attach them to the protomodel """

        self.log ( f"now find highest significance for {len(predictions)} predictions" )
        ## find highest observed significance
        #(set mumax just slightly below its value, so muhat is always below)
        mumax = protomodel.mumax
        bestCombo,Z,llhd,muhat = self.combiner.findHighestSignificance ( predictions,
                strategy, expected=False, mumax = mumax )
        prior = self.combiner.computePrior ( protomodel )
        if hasattr ( protomodel, "keep_meta" ) and protomodel.keep_meta:
            protomodel.bestCombo = bestCombo
        else:
            protomodel.bestCombo = self.combiner.removeDataFromBestCombo ( bestCombo )
        protomodel.Z = Z

        if Z is not None: # Z is None when no combination was found
            protomodel.K = self.combiner.computeK ( Z, prior )
        else:
            protomodel.K = None
        protomodel.llhd = llhd
        protomodel.muhat = muhat
        protomodel.letters = self.combiner.getLetterCode(protomodel.bestCombo)
        protomodel.description = self.combiner.getComboDescription(protomodel.bestCombo)

if __name__ == "__main__":
    import argparse
    argparser = argparse.ArgumentParser(
            description='calling the predictor from the command line' )
    argparser.add_argument ( '-d', '--database',
            help='Path to database pickle file [./default.pcl]',
            type=str, default="./default.pcl" )
    argparser.add_argument ( '-f', '--hiscore',
            help='input model, to be taken from hiscore file [./hiscore.hi]',
            type=str, default="./model.slha" )
    argparser.add_argument ( '-o', '--output',
            help='output pickle file [./predictions.pcl]',
            type=str, default="./predictions.pcl" )
    argparser.add_argument ( '-i', '--interactive',
            help='interactive shell',
            action="store_true" )
    args = argparser.parse_args()
    
    p = Predictor ( 0, args.database, do_combine=False ) 

    sys.path.insert(0,"../")
    from walker.hiscore import Hiscore
    hiscore = Hiscore ( 0, False, args.hiscore )
    hi = hiscore.hiscores[0]
    print ( "Will scrutinize hiscore obtained with database %s K=%.3f" % \
            ( hi.dbversion, hi.K ) )
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

