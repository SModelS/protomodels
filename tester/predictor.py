#!/usr/bin/env python3

""" store the theory predictions in pickle """

import pickle, time, os, sys
from smodels.theory import decomposer
from smodels.theory.theoryPrediction import theoryPredictionsFor
from smodels.share.models.SMparticles import SMList
from smodels.particlesLoader import BSMList
from smodels.tools.physicsUnits import fb, GeV
from smodels.experiment.databaseObj import Database
from smodels.theory.model import Model
from smodels.tools import runtime
runtime._cap_likelihoods = True

try:
    from tester.combiner import Combiner
except:
    from combiner import Combiner

class Predictor:
    def __init__ ( self, walkerid, dbpath = "./default.pcl",
                   expected = False, select = "all",
                   do_combine = False ):
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

    def predict ( self, protomodel, sigmacut = 0.02*fb,
                  strategy = "aggressive", keep_predictions = False ):
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
        bestpreds = self.runSModelS( slhafile, sigmacut,  allpreds=False,
                                           llhdonly=False )

        if keep_predictions:
            self.bestpreds = bestpreds
        # Extract the relevant prediction information and store in the protomodel:
        self.updateModelPredictions(protomodel,bestpreds)
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
            self.log ( "done with prediction. Could not find combinations (Z=%s)" % ( protomodel.Z) )
            protomodel.delCurrentSLHA()
            return False
        else:
            self.log ( "done with prediction. best Z=%.2f (muhat=%.2f)" % ( protomodel.Z, protomodel.muhat ) )

        protomodel.cleanBestCombo()

        #Recompute predictions with higher accuracy for high score models:
        if protomodel.Z > 2.7 and protomodel.nevents < 55000:
            protomodel.nevents = 100000
            protomodel.computeXSecs()
            self.predict(protomodel,sigmacut=sigmacut, strategy= strategy)

        protomodel.delCurrentSLHA()
        return True

    def runSModelS(self, inputFile, sigmacut, allpreds, llhdonly):
        """ run smodels proper.
        :param inputFile: the input slha file
        :param sigmacut: the cut on the topology weights, typically 0.02*fb
        :param allpreds: if true, return all predictions of analyses, else
                         only best signal region
        :param llhdonly: if true, return only results with likelihoods
        """

        if not os.path.exists ( inputFile ):
            self.pprint ( "error, cannot find inputFile %s" % inputFile )
            return []
        model = Model ( BSMList, SMList )
        model.updateParticles ( inputFile=inputFile )

        mingap=10*GeV

        # self.log ( "Now decomposing" )
        topos = decomposer.decompose ( model, sigmacut, minmassgap=mingap )
        self.log ( "decomposed model into %d topologies." % len(topos) )
            

        if allpreds:
            bestDataSet=False
            combinedRes=False
        else:
            bestDataSet=True
            combinedRes=self.do_combine

        preds = []
        # self.log ( "start getting preds" )
        from smodels.tools import runtime
        runtime._experimental = True
        for expRes in self.listOfExpRes:
            predictions = theoryPredictionsFor ( expRes, topos,
                                                 useBestDataset=bestDataSet,
                                                 combinedResults=combinedRes )
            if predictions == None:
                predictions = []
            if allpreds:
                combpreds = theoryPredictionsFor ( expRes, topos,
                                                   useBestDataset=False,
                                                   combinedResults=self.do_combine )
                if combpreds != None:
                    for c in combpreds:
                        predictions.append ( c )
            for prediction in predictions:
                prediction.computeStatistics()
                if (not llhdonly) or (prediction.likelihood != None):
                    preds.append ( prediction )
        sap = "best preds"
        if allpreds:
            sap = "all preds"
        sllhd = ""
        if llhdonly:
            sllhd = ", llhds only"
        self.log ( "returning %d predictions, %s%s" % \
                   (len(preds),sap, sllhd ) )
        return preds

    def printPredictions ( self ):
        """ if self.predictions exists, pretty print them """
        if hasattr ( self, "predictions" ):
            print ( "[predictor] all predictions (for combiner):" )
            for p in self.predictions:
                print ( " - %s %s, %s %s" % \
                        ( p.analysisId(), p.dataType(), p.dataset.dataInfo.dataId, p.txnames ) )
        if hasattr ( self, "bestpreds" ):
            print ( "[predictor] best SR predictions (for critic):" )
            for p in self.bestpreds:
                print ( " - %s %s, %s %s" % \
                        ( p.analysisId(), p.dataType(), p.dataset.dataInfo.dataId, p.txnames ) )

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

        self.log ( "now find highest significance for %d predictions" % len(predictions) )
        ## find highest observed significance
        #(set mumax just slightly below its value, so muhat is always below)
        mumax = protomodel.mumax
        combiner = self.combiner
        bestCombo,Z,llhd,muhat = combiner.findHighestSignificance ( predictions, strategy,
                                                expected=False, mumax = mumax )
        prior = combiner.computePrior ( protomodel )
        if hasattr ( protomodel, "keep_meta" ) and protomodel.keep_meta:
            protomodel.bestCombo = bestCombo
        else:
            protomodel.bestCombo = combiner.removeDataFromBestCombo ( bestCombo )
        protomodel.Z = Z

        if Z is not None: # Z is None when no combination was found
            protomodel.K = combiner.computeK ( Z, prior )
        else:
            protomodel.K = None
        protomodel.llhd = llhd
        protomodel.muhat = muhat
        protomodel.letters = combiner.getLetterCode(protomodel.bestCombo)
        protomodel.description = combiner.getComboDescription(protomodel.bestCombo)

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
    
    p = Predictor ( 0, args.database ) 

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

