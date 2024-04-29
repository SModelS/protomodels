import os
import numpy as np
from typing import Optional
from pathlib import Path
from ptools.expResModifier import ExpResModifier
from smodels.theory.theoryPrediction import theoryPredictionsFor
from smodels.tools.physicsUnits import GeV, fb
from smodels.share.models.mssm import BSMList
from smodels.share.models.SMparticles import SMList
from smodels.theory.model import Model
from smodels.theory.decomposer import decompose
from multiprocessing import Process, Queue

pmodel_path = Path(__file__).parent / 'pmodels/'


def get_seeds(num: int, seed_seed: int = 65536):
    np.random.seed(int(seed_seed))
    return np.random.randint(0, int(1e9), size=num)


def get_pseudodata_args(database: str, tag: str, seed: float = None):
    args = {'dbpath': database,
            'max': 100,
            'rundir': os.getcwd(),
            'keep': False,
            'suffix': tag,
            'lognormal': False,
            'fixedsignals': False,
            'fixedbackgrounds': False,
            'seed': seed,
            'maxmassdist': 400.,
            'compute_ps': False
            }
    return args


def gen_pseudodata(database: str, model: Optional[list[str]] = None, datasetIDs: Optional[list[str]] = None,
                   analysisIDs: Optional[list[str]] = None, seed: Optional[float] = None,
                   bootstrap_num: int = 1) -> list[ExpResModifier]:
    """
    Generate pseudo Experimental results using "ExpResModifier"

    Args:
        database (str): SModels Database containing experimental data
        model (Optional[list[str]], optional): List of models (Topologies) provided to SModelS. Defaults to ['All'].
        datasetIDs (Optional[list[str]], optional): List of dataset IDs. Defaults to ['All'].
        analysisIDs (Optional[list[str]], optional):List of Analysis IDs. Defaults to ['All'].
        seed (float, optional): Random seed. Defaults to None.
        bootstrap_num (int, optional):  Number of bootstrapped sets. Defaults to 1.

    Returns:
        list[ExpResModifier]: List of bootstrapped Experimental results under the SM hypothesis
    """

    if analysisIDs is None:
        analysisIDs = ['all']
    if datasetIDs is None:
        datasetIDs = ['all']
    if model is None:
        model = ['all']

    expResArgs = dict(analysisIDs=analysisIDs, datasetIDs=datasetIDs, txnames=model)
    args = get_pseudodata_args(database, tag='fake', seed=seed)
    modifier = ExpResModifier(args)
    modifier.filter(expResArgs=expResArgs)
    updatedListOfExpRes = []
    for _ in range(bootstrap_num):
        listOfExpRes = modifier.removeEmpty(modifier.db.expResultList)
        updatedListOfExpRes.append(modifier.fakeBackgrounds(listOfExpRes))
    return updatedListOfExpRes


def gen_llr(database: str, slhafile: str, labels: list, model: Optional[list[str]] = None,
            datasetIDs: Optional[list[str]] = None, analysisIDs: Optional[list[str]] = None,
            seed: Optional[float] = None, bootstrap_num: int = 1, expected: bool = False,
            anomaly_mode: bool = True) -> list[dict[str, float]]:
    """
    Generate pseudo NLLR using "ExpResModifier"

    Args:
        database (str): SModels Database containing experimental data
        slhafile (str): Path to slha file (str)
        labels (list): List of SMoldels Unique SR ID's ({analysisId:}:{dataId})
        model (Optional[list[str]], optional): List of models (Topologies) provided to SModelS. Defaults to ['All'].
        datasetIDs (Optional[list[str]], optional): List of dataset IDs. Defaults to ['All'].
        analysisIDs (Optional[list[str]], optional):List of Analysis IDs. Defaults to ['All'].
        seed (float, optional): Random seed. Defaults to None.
        bootstrap_num (int, optional): Number of bootstrapped sets. Defaults to 1.
        expected (bool, optional): run with expected or observed data. Defaults to False.
        anomaly_mode (bool, optional): Run in aomaly mode. Defaults to True.

    Returns:
        list: List[Dict[str, float]] of bootstraped NLLR results under the SM hypothesis
    """

    if analysisIDs is None:
        analysisIDs = ['all']
    if datasetIDs is None:
        datasetIDs = ['all']
    if model is None:
        model = ['all']

    expResArgs = dict(analysisIDs=analysisIDs, datasetIDs=datasetIDs, txnames=model)
    args = get_pseudodata_args(database, tag='fake', seed=seed)
    modifier = ExpResModifier(args)
    modifier.filter(expResArgs=expResArgs)
    llr_dict = []
    for _ in range(bootstrap_num):
        listOfExpRes = modifier.removeEmpty(modifier.db.expResultList)
        exp_results = modifier.fakeBackgrounds(listOfExpRes)
        llr_dict.append(get_llr_at_point(slhafile, exp_results, labels, expected, anomaly_mode))
    return llr_dict


def get_llr_at_point(slhafile: str, exp_results: list, labels: list, expected: bool = False,
                     anomaly_mode: bool = True) -> dict[str, float]:

    """
    Calculate the log likelihood ratio for each SR as specified by model point.

    Args:
        slhafile (str): Str path to slha file
        exp_results (list): List SModels Experimental Results
        labels (list): List of identifying lables -> {analysisId:}:{dataId}
        expected (bool, optional): bool expected (True) or observed (False). Defaults to False.
        anomaly_mode (bool, optional): Run in aomaly mode, L(u=0)/L(u=1). Defaults to True.

    Returns:
        dict[str, float]: Dict of  NLLR results under the SM hypothesis
    """

    min_val = np.finfo(float).tiny
    model = Model(BSMparticles=BSMList, SMparticles=SMList)
    model.updateParticles(inputFile=slhafile)
    toplist = decompose(model, 0.005*fb, doCompress=True, doInvisible=True, minmassgap=5.*GeV)

    log_lilekihood_ratio = {}
    for e_r in exp_results:
        predictions = theoryPredictionsFor(e_r, toplist, useBestDataset=False, combinedResults=False)
        if predictions is not None:
            for pred in predictions:
                tag = f"{e_r.globalInfo.id}:{pred.dataset.dataInfo.dataId}"
                if tag not in labels:
                    continue
                if anomaly_mode:
                    pred.computeStatistics(expected=False)
                    l_Upper = pred.lsm(expected=False)
                    l_lower = pred.likelihood(expected=False)
                else:
                    pred.computeStatistics(expected=expected)
                    l_Upper = pred.likelihood(expected=expected)
                    l_lower = pred.lmax(expected=expected)
                if l_Upper is not None and l_lower is not None:
                    if l_Upper > min_val and l_lower > min_val:
                        log_lilekihood_ratio[tag] = -2 * np.log(l_Upper / l_lower)
                    elif l_Upper <= min_val and l_lower <= min_val:
                        log_lilekihood_ratio[tag] = 0.0
                    else:
                        log_lilekihood_ratio[tag] = np.NaN
    return log_lilekihood_ratio


def split_chunks(num: int, proc: int) -> list[int]:
    """
    Split N tasks (num) between P CPU's (proc) Parameters
    Args:
        num (int): number of tasks
        proc (int): number of CPU's

    Returns:
        list[int]: List of length N CPU's with M tasks per element (sum(M_i) = num)
    """

    if num > proc:
        run_chunks = proc * [num // proc]
        for i in range(num % proc):
            run_chunks[i] += 1
    else:
        run_chunks = num * [1]
    return run_chunks


# def strip_dict_list(list_of_dict: list, labels: list) -> NDArray:
#     """_summary_

#     Args:
#         list_of_dict (list): _description_
#         labels (list): _description_

#     Returns:
#         NDArray: _description_
#     """
#     ret = np.zeros((len(list_of_dict), len(list_of_dict[0])))
#     for i, element in enumerate(list_of_dict):
#         for key, item in element.items():
#             j = labels.index(key)
#             ret[i, j] = item
#     return ret


def get_pseudo_llr(slha_loc: str, data_base: str, labels: list, datasetIDs: Optional[list[str]] = None,
                   analysisIDs: Optional[list[str]] = None, num: int = 1, proc: int = 1) -> list:
    """_summary_

    Args:
        slha_loc (str): _description_
        data_base (str): _description_
        labels (list): _description_
        datasetIDs (Optional[list[str]], optional): _description_. Defaults to None.
        analysisIDs (Optional[list[str]], optional): _description_. Defaults to None.
        num (int, optional): _description_. Defaults to 1.
        proc (int, optional): _description_. Defaults to 1.

    Returns:
        list: _description_
    """

    if analysisIDs is None:
        analysisIDs = ['all']
    if datasetIDs is None:
        datasetIDs = ['all']

    def _llr_worker(args: dict, queue: Queue) -> None:
        """ Helper function to create queue
        Args:
            args (dict): Dictionary of arguments passed to gen_llr
            queue (Queue): Input Queue
        """
        queue.put(gen_llr(**args))

    args = dict(database=data_base, labels=labels, slhafile=slha_loc, datasetIDs=datasetIDs,
                analysisIDs=analysisIDs, seed=None, num=num)
    if proc < 2 or num == 1:
        args['seed'] = np.random.randint(0, 1e6)
        return gen_llr(**args)
    else:
        jobs = []
        input_queue = Queue()
        for item in split_chunks(num, proc):
            args['num'] = item
            args['seed'] = np.random.randint(0, 1e6)
            p = Process(target=_llr_worker, args=(args, input_queue))
            jobs.append(p)
            p.start()
        result_list = [input_queue.get() for _ in range(proc)]
        for p in jobs:
            p.join()
        return [item for sublist in result_list for item in sublist]
