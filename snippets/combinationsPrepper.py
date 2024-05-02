#!/usr/bin/env python3

""" 
.. module:: combinationsPrepper
   :synopsis: simple code snippet for jamie that takes a list of 
   TheoryPredictions and sets up everything for the actual path finder.

.. moduleauthor:: Wolfgang Waltenberger <wolfgang.waltenberger@gmail.com>

"""

from typing import Dict, List
from smodels.matching.theoryPrediction import TheoryPrediction 
import numpy as np

def bamAndWeights ( theorypredictions : List[TheoryPrediction] ) -> Dict:
    """ a simple function that takes a list of theory predictions,
    and from this compute a small binary acceptance matrix (bam) in the guise
    of a dictionary, returns the bam alongside with the dictionary of weights

    :returns: dictionary of bam and weights
    """
    def getTPName ( tpred : TheoryPrediction ) -> str:
        """ get the canonical name of a theory prediction: anaid:datasetid  """
        anaId = tpred.dataset.globalInfo.id
        dsId = "combined"
        if hasattr ( tpred.dataset, "dataInfo" ):
            dsId = tpred.dataset.dataInfo.dataId
        tpId = f"{anaId}:{dsId}"
        return tpId

    bam, weights = {}, {}
    for i,tpred in enumerate(theorypredictions):
        nll0 = tpred.lsm ( return_nll = True )
        nll1 = tpred.likelihood ( return_nll = True )
        w = float("nan")
        if nll0 != None and nll1 != None:
            w = 2 * ( nll1 - nll0 )
        tpId = getTPName ( tpred )
        weights[tpId]=w
        if not tpId in bam:
            bam[tpId]={}
        for tpred2 in theorypredictions[i+1:]:
            tpId2 = getTPName ( tpred )
            combinable = tpred.dataset.isCombinableWith ( tpred2.dataset )
            bam[tpId][tpId2] = combinable
            bam[tpId2][tpId] = combinable
    return { "weights": weights, "bam": bam }

if __name__ == "__main__":
    from smodels.experiment.databaseObj import Database
    from smodels.base.model import Model
    from smodels.decomposition import decomposer
    from smodels.particlesLoader import load
    from smodels.base.physicsUnits import fb, GeV, TeV
    from smodels.matching.theoryPrediction import theoryPredictionsFor
    from smodels.share.models.SMparticles import SMList
    database = Database("official")
    database.getExpResults ( dataTypes = [ "efficiencyMap" ] )
    BSMList = load()
    model = Model(BSMparticles=BSMList, SMparticles=SMList)
    slhafile = "./test.slha"
    model.updateParticles(inputFile=slhafile)
    topDict = decomposer.decompose(model, minmassgap=5*GeV, sigmacut=.1*fb)
    allPredictions = theoryPredictionsFor(database, topDict, combinedResults=True )
    for pred in allPredictions:
        print ( f"prediction {pred} {pred.dataset}" )
    ret = bamAndWeights ( allPredictions )
    print ( ret )
