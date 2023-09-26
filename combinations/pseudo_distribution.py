import numpy as np
import pathfinder as pf
from multiprocessing import Process, Queue
from scipy.stats import norm

from iminuit import Minuit
from iminuit.cost import UnbinnedNLL

def norm_max_pdf(data:np.ndarray, mu:float, sig:float, num:int)->np.ndarray:
    """
    Asymptotic distribution function of the given RV: 
    {Y_0, Y_1, ... Y_N} where Y = MAX{X_0, X_1, ... X_i}_0

    Where X_i is a normaly distributed random variable with

    Parameters
    ----------
    mu  :   location parameter, 
    sig :   scale parameter, 
    num :   sample size (N)

    Returns
    -------
    pdf : ndarray
     Distribution function evaluated at `y`

    """

    return num*norm.pdf(data, loc=mu, scale=sig) * (norm.cdf(data, loc=mu, scale=sig)**(num - 1))

def norm_max_cdf(data:np.ndarray, mu:float, sig:float, num:int):
    
    """
    Cumulative distribution function of the given RV: 
    {Y_0, Y_1, ... Y_N} where Y = MAX{X_0, X_1, ... X_i}_0

    Where X_i is a normaly distributed random variable with

    Parameters
    ----------
    mu  :   location parameter, 
    sig :   scale parameter, 
    num :   sample size (N)

    Returns
    -------
    cdf : ndarray
     Cumulative distribution function evaluated at `y`

    """
    return norm.cdf(data, loc=mu, scale=sig)**num

def data_fit(res:np.ndarray, num:int, mu:float=0, sig:float=0.5, fixnum:bool=True) -> Minuit:
    """
    Fit results to Asymptotic distribution (norm_max_pdf)

    Parameters
    ----------
    res     :   array of results sum of weights from best 
    num     :   sample size (N)
    mu      :   location parameter, 
    sig     :   scale parameter, 
    fixnum  : Fix the N parameter in fitting Default: True

    Returns
    -------
    im_fit : Minuit Fit results

    """
    cost = UnbinnedNLL(res, norm_max_pdf)
    im_fit = Minuit(cost, mu, sig, num)
    
    im_fit.fixed['mu'] = False
    im_fit.fixed['sig'] = False
    im_fit.fixed['num'] = fixnum
    im_fit.migrad()
    #im_fit.scan()
    #im_fit.simplex()
    im_fit.hesse()
    return im_fit

def fit_from_minuit(coefs:Minuit, xmin:float=-5, xmax:float=5):
    """
    converts the IMinuit fit into x and y values for plot
    """
    xval = np.linspace(xmin, xmax, 10000)
    yval = norm_max_pdf(xval, mu=coefs.values['mu'], sig=coefs.values['sig'], num=coefs.values['num'])
    return xval, yval

def get_best_set(binacc:np.ndarray, nllr:np.ndarray, )-> dict:
    """
    Finds best combination form binary acceptence matrix and list of weights
    Parameters
    ----------
    binacc  :   binary acceptance matrix, 2D Array[bool, bool] (N, N) 
    nllr :   NLLR values (weights), 1D array 

    Returns
    -------
    dict: path associated and accosoated weight
    """

    bam = pf.BinaryAcceptance(np.copy(binacc), weights=nllr)
    index_map = bam.sort_bam_by_weight()
    whdfs = pf.WHDFS(bam, top=1, ignore_subset=True)
    whdfs.find_paths(verbose=False, runs=50)
    return {'path': [index_map[i] for i in whdfs.best.path], 'weight': whdfs.best.weight}

def get_milti_bset_set(binacc:np.ndarray, nllr:np.ndarray)-> np.ndarray:
    
    """
    Itterate through 2D array of NLLR values finding the best combination for 
    each row using the corresponding binary acceptance matrix

    Parameters
    ----------
    binacc  :   binary acceptance matrix, 2D Array[bool, bool] (N, N) 
    nllr :   NLLR values (weights), 2D Array (N, M)  

    Returns
    -------
    results : ndarray
        list of Mx2 results, length and weighted sum of best combination from each NLLR set.

    """
    result = np.zeros((len(nllr), 2))
    for i, row in enumerate(nllr):
        res = get_best_set(binacc, row)
        result[i, 0] = res['weight']
        result[i, 1] = len(res['path'])
    return result
        

def best_set_worker(binacc:np.ndarray, nllr:np.ndarray, queue:Queue): #-> Connection:

    """
    Helper function to build queue 
    """
    queue.put(get_milti_bset_set(binacc, nllr))

def find_best_sets(bin_acc:np.ndarray, nllr_dat:np.ndarray, num_cor:int=1)-> np.ndarray:

    """
    Itterate through 2D array of NLLR values finding the best combination for 
    each row using the corresponding binary acceptance matrix

    Parameters
    ----------
    bin_acc  :   binary acceptance matrix, 2D Array[bool, bool] (N, N) 
    nllr_dat :   NLLR values (weights), 2D Array (N, M)  
    num_cor  :   Number of cores Default 1

    Returns
    -------
    results : ndarray
        list of M results, sum of best combination weights from each NLLR set.

    """
    if nllr_dat.ndim == 1:
        nllr_dat = nllr_dat[np.newaxis, ...]
    if num_cor < 2:
        print(F"Starting job 1. Calculating {len(nllr_dat)} best combinations")
        return get_milti_bset_set(bin_acc, nllr_dat)
    else:
        jobs = []
        input_queue = Queue(maxsize=num_cor)
        for i, chunk in enumerate(np.array_split(nllr_dat, num_cor, axis=0)):
            p = Process(target=best_set_worker, args=(bin_acc, chunk, input_queue))
            jobs.append(p)
            p.start()
            print(F"Starting job {i+1}. Calculating {len(chunk)} best combinations")
        result_list = [input_queue.get() for _ in range(num_cor)]
        for p in jobs:
            p.join()
        return np.array([item for sublist in result_list for item in sublist])