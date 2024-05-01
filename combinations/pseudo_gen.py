import os
import numpy as np
from typing import Optional, Union
from pathlib import Path
from ptools.expResModifier import ExpResModifier
from smodels.experiment.databaseObj import Database
from smodels.matching.theoryPrediction import theoryPredictionsFor, TheoryPredictionList, TheoryPrediction
from smodels.base.physicsUnits import GeV, fb
from smodels.share.models.mssm import BSMList
from smodels.share.models.SMparticles import SMList
from smodels.base.model import Model
from smodels.decomposition.decomposer import decompose
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


def gen_llr(database: str, slhafile: str, labels: list, analysisIDs: list[str],
            model: Optional[list[str]] = None, seed: Optional[float] = None, bootstrap_num: int = 1,
            expected: bool = False, anomaly_mode: bool = True) -> list[dict[str, float]]:
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

    if model is None:
        model = ['all']

    # expResArgs = dict(analysisIDs=analysisIDs, datasetIDs=['all'], txnames=model)
    args = get_pseudodata_args(database, tag='fake', seed=seed)
    modifier = ExpResModifier(args)
    modifier.filter()
    llr_dict = []
    for _ in range(bootstrap_num):
        listOfExpRes = modifier.removeEmpty(modifier.db.expResultList)
        pseudo_databse = {'database': modifier.db, 'expResults': modifier.fakeBackgrounds(listOfExpRes)}
        llr_dict.append(get_llr_at_point(slhafile, analysisIDs, expected,
                                         anomaly_mode, pseudo_databse=pseudo_databse))
    return llr_dict


def get_prediction(pred: Union[TheoryPredictionList, TheoryPrediction], anomaly_mode: bool,
                   expected: bool, min_val: float = np.finfo(float).tiny) -> float:
    ratio = np.NaN
    if anomaly_mode:
        pred.computeStatistics(expected=False)
        likelihood = pred.lsm(expected=False)
        lmax = pred.likelihood(expected=False)
    else:
        pred.computeStatistics(expected=expected)
        likelihood = pred.likelihood(expected=expected)
        lmax = pred.lmax(expected=expected)

    if likelihood is not None and lmax is not None:

        if likelihood > min_val and lmax > min_val:
            ratio = -2 * np.log(likelihood / lmax)

        elif likelihood <= min_val and lmax <= min_val:
            ratio = 0.0

    return ratio


def get_llr_at_point(slhafile: Union[str, Path], analysisIDs: list[str], expected: bool = False,
                     anomaly_mode: bool = False, pseudo_databse: Optional[dict[str, Database, list]] = None) -> dict:

    srs_in_ids = any([":" in item for item in analysisIDs])
    if srs_in_ids:
        print("Signal regions found in analysis ID's.")
        print("Searching data ID's for possible matches.")
    model = Model(BSMparticles=BSMList, SMparticles=SMList)
    model.updateParticles(inputFile=slhafile)
    top_dict = decompose(model, sigmacut=0.005*fb, massCompress=True, invisibleCompress=True, minmassgap=5*GeV)
    if pseudo_databse is None:
        dbase = Database("official")
        expResults = dbase.getExpResults(analysisIDs=['all'], datasetIDs=['all'], dataTypes=['efficiencyMap'])
    else:
        dbase = pseudo_databse['database']
        expResults = pseudo_databse['expResults']
    allThPredictions = theoryPredictionsFor(dbase, top_dict)
    llr_dict = {}
    for expRes, thPreds in zip(expResults, allThPredictions):
        if not thPreds:
            continue
        aid_label = expRes.globalInfo.id
        if aid_label in analysisIDs:
            llr = get_prediction(thPreds, anomaly_mode=anomaly_mode, expected=expected)
            if not np.isnan(llr):
                llr_dict[aid_label] = llr
        elif srs_in_ids:
            for dataset in expRes.datasets:
                if hasattr(dataset.dataInfo, 'dataId'):
                    sms_label = f"{expRes.globalInfo.id}:{dataset.dataInfo.dataId}"
                    if sms_label in analysisIDs:
                        llr = get_prediction(thPreds, anomaly_mode=anomaly_mode, expected=expected)
                        if not np.isnan(llr):
                            llr_dict[sms_label] = llr
    return llr_dict


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


def strip_dict_list(list_of_dict: list, labels: list) -> np.ndarray:
    """_summary_

    Args:
        list_of_dict (list): _description_
        labels (list): _description_

    Returns:
        NDArray: _description_
    """
    ret = np.zeros((len(list_of_dict), len(labels)))
    for i, element in enumerate(list_of_dict):
        for key, item in element.items():
            if key in labels:
                j = labels.index(key)
                ret[i, j] = item
    return ret


def _llr_worker(args: dict, queue: Queue) -> None:
    """ Helper function to create queue
    Args:
        args (dict): Dictionary of arguments passed to gen_llr
        queue (Queue): Input Queue
    """
    queue.put(gen_llr(**args))


def get_pseudo_llr(slha_loc: str, data_base: str, labels: list, analysisIDs: list[str],
                   num: int = 1, proc: int = 1) -> list:
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

    args = dict(database=data_base, labels=labels, slhafile=slha_loc,
                analysisIDs=analysisIDs)
    if proc < 2 or num == 1:
        args['seed'] = np.random.randint(0, 1e6)
        return gen_llr(**args)
    else:
        jobs = []
        input_queue = Queue()
        for item in split_chunks(num, proc):
            # args['num'] = item
            args['seed'] = np.random.randint(0, 1e6)
            p = Process(target=_llr_worker, args=(args, input_queue))
            jobs.append(p)
            p.start()
        result_list = [input_queue.get() for _ in range(proc)]
        for p in jobs:
            p.join()
        return [item for sublist in result_list for item in sublist]
