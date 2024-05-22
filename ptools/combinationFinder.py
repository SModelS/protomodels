import numpy as np
try:
    import pathfinder as pf
except ImportError as e:
    # FIXME in the long run the line below should disappear
    import sys,os
    sys.path.insert(0, os.path.expanduser("~/PathFinder"))
    sys.path.insert(0, os.path.expanduser("~/git/PathFinder"))
    import pathfinder as pf
from typing import Iterable, Dict, List, Optional, Union
from numpy.typing import NDArray
from smodels.matching.theoryPrediction import theoryPredictionsFor, TheoryPrediction

__all__ = [ "selectMostSignificantSRs", "bamAndWeights", "find_best_comb"  ]

def selectMostSignificantSRs ( predictions: list[TheoryPrediction], bound: float = 5 ) -> List:
    """
    Given the predictions, return the "x" most significant SRs per analysis.

    :param predictions: all predictions of all SRs
    :param bound: an upper bound on the percentage of ratios of weights of SRs per analysis
    :returns: list of predictions of "x" most significant SRs of each analysis
    """
    sortByAnaId = {}                             # first sort all by ana id + data Type
    for pred in predictions:
        Id = pred.analysisId()+":"+pred.dataType(True)
        if not Id in sortByAnaId:
            sortByAnaId[Id]=[]
        sortByAnaId[Id].append ( pred )         #keep all em-type ds of one analysis under one key

    ret = []
    #keptThese = [] ## log the ana ids that we kept, for debugging only.
    for Id,preds in sortByAnaId.items():
        if len(preds) == 1:                 # If only 1 prediction, use it (could be combined result, or only one em result)
            ret.append( preds[0] )
            #keptThese.append ( preds[0].experimentalId() )   #self.getPredictionID ( pred )
            continue

        maxRatio, ratioList = 0., {}
        for pred in preds:
            nll0 = pred.likelihood(mu=0, expected=False, return_nll=True)
            nll1 = pred.likelihood(mu=1, expected=False, return_nll=True)
            if nll0 is None or nll1 is None: continue
            ratio = 2 * (nll0 - nll1)
            ratioList[pred] = ratio
            if ratio > maxRatio:
                maxRatio = ratio

        ratioList = {k:v for k,v in sorted(ratioList.items(), key=lambda item:item[1], reverse=True)}
        signPreds = [pred for pred, ratio in ratioList.items() if ratio/maxRatio >= (bound/100.) ]

        if len(signPreds) == 1:
            secondPred = list(ratioList.keys())[1]
            signPreds.append(secondPred)

        ret = ret + signPreds
        #keptThese = keptThese + [pred.experimentalId() for pred in signPreds]
        #self.pprint ( f"selected predictions down via SRs from {len(predictions)}"\f" to {len(ret)}." )

    return ret

def bamAndWeights(theorypredictions: list[TheoryPrediction], expected: bool = False, excl_mode: bool = False) -> Dict:
    """
    A simple function that takes a list of theory predictions,
    and from this compute a small binary acceptance matrix (bam) in the guise
    of a dictionary, returns the bam alongside with the dictionary of weights

    :returns: dictionary of bam and weights
    """
    from ptools.helpers import experimentalId

    bam, weights, theoryPred = {}, {}, {}

    for i, tpred in enumerate(theorypredictions):
        # nll0 = tpred.lsm(expected=expected, return_nll=True)
        nll0 = tpred.likelihood(mu=0, expected=expected, return_nll=True)
        nll1 = tpred.likelihood(mu=1, expected=expected, return_nll=True)
        w = np.NaN
        if nll0 is not None and nll1 is not None:
            # w = -2 * (ll0 - ll1) = 2 * (ll1 - ll0) = 2 * (-ll0 - (-ll1)) = 2 * (nll0 - nll1) for anamoly mode
            w = 2 * (nll0 - nll1)
            if excl_mode: w = -2 * (nll0 - nll1)

        tpId = experimentalId(tpred)
        weights[tpId] = w
        theoryPred[tpId] = tpred
        if tpId not in bam:
            bam[tpId] = set()
        for tpred2 in theorypredictions[i+1:]:
            tpId2 = experimentalId(tpred2)
            if tpred.dataset.isCombinableWith(tpred2.dataset):
                bam[tpId].add(tpId2)

    return {"weights": weights, "bam": bam, "theoryPred": theoryPred}

def get_bam_weight(over: Dict, weight: Dict) -> Dict[str, NDArray]:
    """
    Construct the binary_acceptance_matrix and weights from the the Dictionary type provided by SModelS.

    Args:
        over (Dict): SModelS correlation type dictionary in the form: {'label': [list of combinable labels]}
        weight (Dict): list of weights to chose combination in the form: {'label': weight (float)}

    Returns:
        Dict[str, NDArray, List]:
            bam (NDArray) -> Binary Acceptance Matrix: an N x N Symmetric True/False matrix that defines the allowable combinations
            weights (NDArray) -> 1 x N Weights that correspond to the binary acceptance matrix indices
            labels (List) -> 1 x N list of labels that match the binary acceptance and weight indices
    """
    columns_labels = list(weight.keys())
    bam = np.zeros((len(columns_labels), len(columns_labels)), dtype=bool)
    for i, key in enumerate(columns_labels):
        bam[i, :] = [True if sr in over[key] else False for sr in columns_labels]

    if not np.allclose(bam, bam.T):
        print("ERROR: Bam not symmetric!")      #move to loggerbase comment later on when loggerbase.py is moved

    bam |= np.triu(bam).T                       #ensure matrix is symmetric

    weight_array = np.array([item for _, item in weight.items()])
    order = np.argsort(weight_array)[::-1]
    return {'bam': bam[order, :][:, order],
            'weights': weight_array[order],
            'labels': [columns_labels[i] for i in order]}


def get_best_set(binary_acceptance_matrix: NDArray, weights: NDArray, sort_bam=False) -> Dict[str, NDArray]:
    """
    Get the highest sum weights that can be combined according to the binary acceptance matrix.

    Args:
        binary_acceptance_matrix (NDArray): N x N Symmetric matrix that defines the allowable combinations
        weights (NDArray): 1 x N Weights that correspond to the binary acceptance matrix indices
        sort_bam (bool, optional): Sort the binary acceptance matrix by weight (descending). Defaults to False.

    Returns:
        Dict[str, NDArray]: Containing the combination path indices and sum of weight sum.
    """
    weights -= 1                                                #check later when decided
    offset = 0.0
    if min(weights) < 0.0:
        offset = abs(min(weights)) + 1
    bam = pf.BinaryAcceptance(binary_acceptance_matrix, weights=weights + offset)
    results = {}
    if sort_bam:
        results['order'] = bam.sort_bam_by_weight()             #?? sorting because weights changed, how do use order again?
    whdfs = pf.WHDFS(bam, top=1, ignore_subset=True)
    whdfs.find_paths(verbose=False, runs=50)
    results['path'] = whdfs.best.path
    results['weight'] = whdfs.best.weight - (len(whdfs.best.path) * offset) + 1.0
    results['offset'] = offset
    return results


def find_best_comb(bam_weight_dict: Dict) -> Dict[str, float]:
    """
    Find the best combination of srs given the bam_weight_dict

    Args:
        bam_weight_dict (Dict): Correlation/overlap and weight information gathered from the SModelS API

    Returns:
        Dict[str, float, int]: Dictionary of list of best srs and corresponding weight for the given bam_weight_dict
    """

    bam_wgths = get_bam_weight(bam_weight_dict['bam'], bam_weight_dict['weights'])                #get bam and corresponfding weights
    temp_res = get_best_set(bam_wgths['bam'], bam_wgths['weights'])         #get best path for bam and weight
    best_labels = [bam_wgths['labels'][p] for p in temp_res['path']]        #get corresponding srs in best path

    result = {'best': best_labels, 'weight': temp_res['weight'], 'offset': temp_res['offset']}   #store labels and weight in dict
    return result
