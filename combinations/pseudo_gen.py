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
from multiprocessing import Process, Manager

pmodel_path = Path(__file__).parent / 'pmodels/'


def get_seeds(num: int, seed_seed: int = 65536):
    np.random.seed(int(seed_seed))
    return np.random.randint(0, int(1e9), size=num)


def get_pseudodata_args(database: str, tag: str, seed: float = None):
    args = {'dbpath': database,
            'max': 100,
            'rundir': os.getcwd(),
            'keep': False,
            'suffix': None,
            'lognormal': False,
            'fixedsignals': False,
            'fixedbackgrounds': False,
            'noupperlimits': True,
            'seed': seed,
            'maxmassdist': 400.,
            'compute_ps': False,
            'outfile': None,
            }
    return args


def gen_llr(database: str, slhafile: str, model: Optional[list[str]] = None, seed: Optional[float] = None,
            bootstrap_num: int = 1) -> list[dict[str, float]]:
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
        llr_dict.append(get_llr_at_point(slhafile, pseudo_databse=pseudo_databse))
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


def get_llr_at_point(slhafile: Union[str, Path], data_base: str = 'official',
                     pseudo_databse: Optional[dict[str, Database, list]] = None) -> dict:

    model = Model(BSMparticles=BSMList, SMparticles=SMList)
    model.updateParticles(inputFile=slhafile)
    top_dict = decompose(model, sigmacut=0.005*fb, massCompress=True, invisibleCompress=True, minmassgap=5*GeV)
    if pseudo_databse is None:
        dbase = Database(data_base)
        _ = dbase.getExpResults(analysisIDs=['all'], datasetIDs=['all'], dataTypes=['efficiencyMap'])
    else:
        dbase = pseudo_databse['database']
        _ = pseudo_databse['expResults']
    allThPredictions = theoryPredictionsFor(dbase, top_dict)
    return bamAndWeights(allThPredictions)


def bamAndWeights(theorypredictions: list[TheoryPrediction]) -> dict:
    """ a simple function that takes a list of theory predictions,
    and from this compute a small binary acceptance matrix (bam) in the guise
    of a dictionary, returns the bam alongside with the dictionary of weights

    :returns: dictionary of bam and weights
    """
    def getTPName(tpred: TheoryPrediction) -> str:
        """ get the canonical name of a theory prediction: anaid:datasetid  """
        anaId = tpred.dataset.globalInfo.id
        dsId = ""
        if hasattr(tpred.dataset, "dataInfo"):
            dsId = f":{tpred.dataset.dataInfo.dataId}"
        tpId = f"{anaId}{dsId}"
        return tpId

    bam, weights = {}, {}
    for i, tpred in enumerate(theorypredictions):
        nll0 = tpred.lsm(return_nll=True)
        nll1 = tpred.likelihood(return_nll=True)
        w = np.NaN
        if nll0 is not None and nll1 is not None:
            w = 2 * (nll1 - nll0)
        tpId = getTPName(tpred)
        weights[tpId] = w
        if tpId not in bam:
            bam[tpId] = set()
        for tpred2 in theorypredictions[i+1:]:
            tpId2 = getTPName(tpred2)
            if tpred.dataset.isCombinableWith(tpred2.dataset):
                bam[tpId].add(tpId2)
    return {"weights": weights, "bam": bam}


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


def _llr_worker(args: dict, outputlist: list) -> None:
    """ Helper function to create queue
    Args:
        args (dict): Dictionary of arguments passed to gen_llr
        queue (Queue): Input Queue
    """
    outputlist.extend(gen_llr(**args))


def get_pseudo_llr(slha_loc: str, data_base: str, bootstrap_num: int = 1, proc: int = 1) -> list[dict]:
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

    args = dict(database=data_base, slhafile=slha_loc)
    if proc < 2 or bootstrap_num == 1:
        args['seed'] = np.random.randint(0, 1e6)
        outputlist = [gen_llr(**args)]
    else:
        jobs = []
        manager = Manager()
        outputlist = manager.list()
        for item in split_chunks(bootstrap_num, proc):
            args['bootstrap_num'] = item
            args['seed'] = np.random.randint(0, 1e6)
            p = Process(target=_llr_worker, args=(args, outputlist))
            jobs.append(p)
            p.start()
        for p in jobs:
            p.join()
    return list(outputlist)
