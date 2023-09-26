import os
import numpy as np
from protomodels.ptools.expResModifier import ExpResModifier
from smodels.theory.theoryPrediction import theoryPredictionsFor
from smodels.tools.physicsUnits import GeV, pb, fb
from smodels.particlesLoader import BSMList
from smodels.share.models.SMparticles import SMList
from smodels.theory.model import Model
from smodels.theory.decomposer import decompose
from multiprocessing import Process, Queue

pmodel_path = 'Model_Data/pmodels/'

def get_seeds(num:int, seed_seed:int=65536):
    np.random.seed(int(seed_seed))
    return np.random.randint(0, int(1e9), size=num)

def get_pseudodata_args(database:str, tag:str, seed:float=None):
    args = {'dbpath':database, 
            'max': 100,
            'rundir':os.getcwd(), 
            'keep':False, 
            'suffix': tag,
            'lognormal': False,
            'fixedsignals': False,
            'fixedbackgrounds': False,
            'seed': seed,
            'maxmassdist':400.,
            'compute_ps': False
            }
    return args     


def gen_pseudodata(database:str, model:list=['all'], datasetIDs:list=['all'], 
                   analysisIDs:list=['all'], seed:float=None, num:int=1)->list:
    """
    Generate pseudo Experimental results using "ExpResModifier"

    Parameters
    ----------
    database    :   SModels Database containing experimental data, 
    model       :   List of models (Topologies) provided to SModelS -> default 'All'
    datasetIDs  :   List of dataset IDs -> Default 'All
    analysisIDs :   List of Analysis IDs -> Default 'All
    seed        :   Random seed
    num         :   Number of bootstrapped sets

    Returns
    -------
    pdf : updatedListOfExpRes
        List of bootstraped Experimental results under the SM hypothesis

    """
    expResArgs = dict(analysisIDs=analysisIDs, datasetIDs=datasetIDs, txnames=model)
    args = get_pseudodata_args(database, tag='fake', seed=seed) #,seed=np.random.randint(0, 1e6))
    modifier = ExpResModifier(args)
    modifier.filter (expResArgs=expResArgs)
    updatedListOfExpRes = []
    for _ in range(num):
        listOfExpRes = modifier.removeEmpty ( modifier.db.expResultList )
        updatedListOfExpRes.append(modifier.fakeBackgrounds (listOfExpRes))
    return updatedListOfExpRes

def gen_llr(database:str, slhafile:str, labels:list, model:list=['all'], datasetIDs:list=['all'], 
            analysisIDs:list=['all'], seed:float=None, num:int=1, expected:bool=False, anomaly_mode:bool=True) -> list:

    """
    Generate pseudo NLLR using "ExpResModifier"

    Parameters
    ----------
    database    :   SModels Database containing experimental data
    slhafile    :   Path to slha file (str)
    labels      :   List of SMoldels Unique SR ID's ({analysisId:}:{dataId})
    model       :   List of models (Topologies) provided to SModelS -> default 'All'
    datasetIDs  :   List of dataset IDs -> Default 'All
    analysisIDs :   List of Analysis IDs -> Default 'All
    seed        :   Random seed
    num         :   Number of bootstrapped sets

    Returns
    -------
    pdf : updatedListOfExpRes
        List[Dict] of bootstraped NLLR results under the SM hypothesis

    """
    expResArgs = dict(analysisIDs=analysisIDs, datasetIDs=datasetIDs, txnames=model)
    args = get_pseudodata_args(database, tag='fake', seed=seed) #,seed=np.random.randint(0, 1e6))
    modifier = ExpResModifier(args)
    modifier.filter (expResArgs=expResArgs)
    llr_dict = []
    for _ in range(num):
        listOfExpRes = modifier.removeEmpty ( modifier.db.expResultList )
        exp_results = modifier.fakeBackgrounds (listOfExpRes)
        llr_dict.append(get_llr_at_point(slhafile, exp_results, labels, expected, anomaly_mode))
    return llr_dict

def get_llr_at_point(slhafile:str, exp_results:list, labels:list, 
                     expected:bool=False, anomaly_mode:bool=True)-> dict:
    
    """
    Calculate the log likelihood ratio for each SR as specified by model point.

    Parameters
    ----------
    
        slhafile    :   Str path to slha file
        labels      :   List of identifying lables -> {analysisId:}:{dataId}
        exp_results :   List SModels Experimental Results 
        expected    :   bool expected (True) or observed (False)
        anomaly_mode:   L(u=0)/L(u=1)

    Returns
    -------
    pdf : updatedListOfExpRes
        Dict of  NLLR results under the SM hypothesis

    """

    min_val = np.finfo(float).tiny
    model = Model(BSMparticles=BSMList, SMparticles=SMList)
    model.updateParticles(inputFile=slhafile)
    toplist = decompose(model, 0.005*fb, doCompress=True,
                                doInvisible=True, minmassgap=5.*GeV)

    llr = {}
    for e_r in exp_results:
        predictions = theoryPredictionsFor(e_r, toplist, useBestDataset=False,
                                        combinedResults=False)
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
                    if l_Upper  > min_val and l_lower > min_val:
                        llr[tag] = -2 * np.log ( l_Upper / l_lower )
                    elif l_Upper <= min_val and l_lower <= min_val:
                        llr[tag] = 0.0
                    else:
                        llr[tag] = np.NaN
    return llr

def split_chunks(num:int, proc:int)->list:
    """
    split N tasks (num) between P CPU's (proc)
    Parameters
    ----------
    num     : int number of tasks
    proc    : int number of CPU's 

    Returns
    -------
    run_chunks : List of length N CPU's with M tasks per element (sum(M_i) = num)
    """
    if num > proc:
        run_chunks = proc * [num // proc]
        for i in range(num % proc):
            run_chunks[i] += 1
    else:
        run_chunks = num * [1]
    return run_chunks

def strip_dict_list(list_of_dict:list, labels:list)->np.ndarray:
    ret = np.zeros((len(list_of_dict), len(list_of_dict[0])))
    for i, element in enumerate(list_of_dict):
        for key, item in element.items():
            j = labels.index(key)
            ret[i, j] = item
    return ret

def llr_worker(args:dict, queue:Queue): #-> Connection:
    queue.put(gen_llr(**args))

def get_pseudo_llr(slha_loc:str, data_base:str, labels:list, analysisIDs:list=['all'], 
                 datasetIDs:list=['all'], num:int=1, proc:int=1)-> dict:
    args = dict(database=data_base, labels=labels, slhafile=slha_loc, 
                datasetIDs=datasetIDs, analysisIDs=analysisIDs,
                seed=0, num=num)
    if proc < 2 or num == 1:
        args['seed']=np.random.randint(0, 1e6)
        return gen_llr(**args)
    else:
        jobs = []
        input_queue = Queue()
        for item in split_chunks(num, proc):
            args['num'] = item
            args['seed']=np.random.randint(0, 1e6)
            p = Process(target=llr_worker, args=(args, input_queue))
            jobs.append(p)
            p.start()
        result_list = [input_queue.get() for _ in range(proc)]
        for p in jobs:
            p.join()
        return [item for sublist in result_list for item in sublist]
