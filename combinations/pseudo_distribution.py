import numpy as np
import pathfinder as pf
from multiprocessing import Process, Queue
from scipy.stats import norm
from typing import Tuple
from numpy.typing import NDArray
from iminuit import Minuit
from iminuit.cost import UnbinnedNLL


def norm_max_pdf(data: NDArray, mu: float, sig: float, num: int) -> NDArray:
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


def norm_max_cdf(data: NDArray, mu: float, sig: float, num: int) -> NDArray:
    """Cumulative distribution function of the given RV:
    {Y_0, Y_1, ... Y_N} where Y = MAX{X_0, X_1, ... X_i}_0
    Where X_i is a normaly distributed random variable with

    Args:
        data (NDArray): array_like quantiles
        mu (float): Location parameter,
        sig (float): Scale parameter
        num (int): Sample size (N)

    Returns:
        NDArray: Cumulative distribution function evaluated at `y`
    """

    return norm.cdf(data, loc=mu, scale=sig)**num


def data_fit(res: NDArray, num: int, mu: float = 0, sig: float = 0.5, fixnum: bool = True) -> Minuit:
    """ Fit results to Asymptotic distribution (norm_max_pdf)

    Args:
        res (NDArray): array of results sum of weights from best
        num (int): sample size (N)_
        mu (float, optional): location parameter. Defaults to 0.
        sig (float, optional): scale parameter. Defaults to 0.5.
        fixnum (bool, optional): Fix the N parameter in fitting. Defaults to True.

    Returns:
        Minuit: Fit results
    """

    nll = UnbinnedNLL(res, norm_max_pdf)
    im_fit = Minuit(nll, mu, sig, num)

    im_fit.fixed['mu'] = False
    im_fit.fixed['sig'] = False
    im_fit.fixed['num'] = fixnum
    im_fit.scan()
    im_fit.hesse()
    return im_fit


def fit_from_minuit(coefs: Minuit, xmin: float = -5, xmax: float = 5) -> Tuple[NDArray, NDArray]:
    """
    Use Minuit fit to get X, Y plot points.
    Args:
        coefs (Minuit): Function minimizer used to fit data
        xmin (float, optional): x value minimum. Defaults to -5.
        xmax (float, optional): x value maximum. Defaults to 5.

    Returns:
        Tuple[NDArray, NDArray] : X and y values of the fit
    """
    xval = np.linspace(xmin, xmax, 10000)
    yval = norm_max_pdf(xval, mu=coefs.values['mu'], sig=coefs.values['sig'], num=coefs.values['num'])
    return xval, yval


def get_best_set(binacc: NDArray, nllr: NDArray) -> dict:
    """_summary_

    Args:
        binacc (NDArray): _description_
        nllr (NDArray): _description_

    Returns:
        dict: _description_
    """
    bam = pf.BinaryAcceptance(np.copy(binacc), weights=nllr)
    bam.sort_bam_by_weight()
    whdfs = pf.WHDFS(bam, top=1, ignore_subset=True)
    whdfs.find_paths(verbose=False, runs=50)
    return {'path': whdfs.best.path, 'weight': whdfs.best.weight}


def get_milti_bset_set(binacc: NDArray, nllr: NDArray) -> NDArray:
    result = np.zeros((len(nllr), 2))
    for i, row in enumerate(nllr):
        res = get_best_set(binacc, row)
        result[i, 0] = res['weight']
        result[i, 1] = len(res['path'])
    return result


def best_set_worker(binacc: NDArray, nllr: NDArray, queue: Queue) -> None:
    queue.put(get_milti_bset_set(binacc, nllr))


def find_best_sets(bin_acc: NDArray, nllr_dat: NDArray, num_cor: int = 1) -> NDArray:

    """
    Iterate through 2D array of NLLR values finding the best combination for
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
