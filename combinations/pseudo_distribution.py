import numpy as np
import pathfinder as pf
from multiprocessing import Process, Manager
from typing import Iterable, Dict, List
from numpy.typing import NDArray


def get_bam_weight(over: Dict, weight: Dict) -> Dict[str, NDArray]:
    """
    Construct the binary_acceptance_matrix and weights from the the Dictionary type provided by SModelS.

    Args:
        over (Dict): SModelS correlation type dictionary in the form: {'label': [list of combinable labels]}
        weight (Dict): list of weights to chose combination in the form: {'label': weight (float)}

    Returns:
        Dict[str, NDArray, List]:
            bam (NDArray) -> Binary Acceptance Matrix: an N x N Symmetric matrix that defines the allowable combinations
            weights (NDArray) -> 1 x N Weights that correspond to the binary acceptance matrix indices
            labels (List) -> 1 x N list of labels that match the binary acceptance and weight indices
    """
    columns_labels = list(weight.keys())
    bam = np.zeros((len(columns_labels), len(columns_labels)), dtype=bool)
    for i, key in enumerate(columns_labels):
        bam[i, :] = [True if sr in over[key] else False for sr in columns_labels]
    if not np.allclose(bam, bam.T):
        print('WARNING bam not symetric')
        bam |= np.triu(bam).T
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
    weights -= 1
    offset = 0.0
    if min(weights) < 0.0:
        offset = abs(min(weights)) + 1
    bam = pf.BinaryAcceptance(binary_acceptance_matrix, weights=weights + offset)
    results = {}
    if sort_bam:
        results['order'] = bam.sort_bam_by_weight()
    whdfs = pf.WHDFS(bam, top=1, ignore_subset=True)
    whdfs.find_paths(verbose=False, runs=50)
    results['path'] = whdfs.best.path
    results['weight'] = whdfs.best.weight - (len(whdfs.best.path) * offset) + 1.0
    results['offset'] = offset
    return results


def get_multi_best_set(pseudo_gen_dicts: List[Dict]) -> Dict[str, float]:
    """
    Iterate through a list of dictionaries containing the dictionaries of corelation and weight information
    gathered from the SModelS API

    Args:
        pseudo_gen_dicts (List[Dict]): Corelation/overlap and weight information gathered from the SModelS API

    Returns:
        Dict[str, float, int]: Dictionary containing
    """
    result = {}
    for i, item in enumerate(pseudo_gen_dicts):
        bam_wgths = get_bam_weight(item['bam'], item['weights'])
        temp_res = get_best_set(bam_wgths['bam'], bam_wgths['weights'])
        best_labels = [bam_wgths['labels'][p] for p in temp_res['path']]
        result[i] = {'best': best_labels, 'weight': temp_res['weight'], 'offset': temp_res['offset']}
    return result


def _best_set_worker(pseudo_gen_dicts: List[Dict], run_num: int, return_dict: Dict,) -> None:
    """
    Multi-processing worker for find_best_sets

    Args:
        pseudo_gen_dicts (List[Dict]): List of dictionaries containing a binary acceptance matrix
                                       and set of corresponding weights.
        run_num (int): Unique integer identifier for labeling return dictionary
        return_dict (Dict): DictProxy for Manager
    """
    for key, item in get_multi_best_set(pseudo_gen_dicts).items():
        idx = (run_num * len(pseudo_gen_dicts)) + key
        return_dict.update({idx: item})


def find_best_sets(pseudo_gen_dicts: List[Dict], num_cor: int = 1) -> Dict[int, Dict]:

    """
    Propagate the get_multi_best_set function over multiple CPU's

    Args:
        pseudo_gen_dicts (List[Dict]): List of dictionaries containing a binary acceptance matrix
                                       and set of corresponding weights.
        num_cor (int): Number of CPU cores to use

    Returns:
        Dict[Dict]: Enumerated result dictionaries containing the best sets with the corresponding weights
    """

    def split_list(list_in: List, nunber_of_chuncks: int) -> Iterable[List]:
        for i in range(0, nunber_of_chuncks):
            yield list_in[i::nunber_of_chuncks]

    if num_cor < 2:
        print(F"Starting job 1. Calculating {len(pseudo_gen_dicts)} best combinations")
        outputdict = get_multi_best_set(pseudo_gen_dicts)
    else:
        jobs = []
        manager = Manager()
        outputdict = manager.dict()
        bam_weights = split_list(pseudo_gen_dicts, num_cor)
        for i, bam_wgts in enumerate(bam_weights):
            p = Process(target=_best_set_worker, args=(bam_wgts, i, outputdict))
            jobs.append(p)
            p.start()
            print(F"Starting job {i+1}. Calculating {len(bam_wgts)} best combinations")
        for p in jobs:
            p.join()
    return dict(sorted(outputdict.items()))
