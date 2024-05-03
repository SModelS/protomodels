import numpy as np
import pathfinder as pf
from multiprocessing import Process, Manager
from scipy.stats import norm
from typing import Tuple, Iterable
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


def get_best_set(binacc: NDArray, nllr: NDArray, sort_bam=False) -> dict[str, NDArray]:
    """_summary_

    Args:
        binacc (NDArray): _description_
        nllr (NDArray): _description_

    Returns:
        dict: _description_
    """
    bam = pf.BinaryAcceptance(np.copy(binacc), weights=nllr)
    results = {}
    if sort_bam:
        results['order'] = bam.sort_bam_by_weight()
    whdfs = pf.WHDFS(bam, top=1, ignore_subset=True)
    whdfs.find_paths(verbose=False, runs=50)
    results['path'] = whdfs.best.path
    results['weight'] = whdfs.best.weight
    return results


def get_bam_weight(over: dict, weight: dict) -> dict:
    columns_labels = list(weight.keys())
    bam = np.zeros((len(columns_labels), len(columns_labels)), dtype=bool)
    for i, key in enumerate(columns_labels):
        bam[i, :] = [True if sr in over[key] else False for sr in columns_labels]
    bam |= np.triu(bam).T
    assert np.allclose(bam, bam.T)
    weight_array = np.array([item for _, item in weight.items()])
    order = np.argsort(weight_array)[::-1]
    return {'bam': bam[order, :][:, order],
            'weights': weight_array[order],
            'labels': [columns_labels[i] for i in order]}


def split_list(list_in: list, nunber_of_chuncks: int) -> Iterable[list]:
    """
    Yield n number of striped chunks from l

    Args:
        list_in (list): _description_
        nunber_of_chuncks (int): _description_

    Returns:
        Iterable: _description_

    Yields:
        Iterator[Iterable]: _description_
    """
    for i in range(0, nunber_of_chuncks):
        yield list_in[i::nunber_of_chuncks]


def get_milti_bset_set(pseudo_gen_dicts: list[dict]) -> dict[str, float, int]:
    restut = {}
    for i, item in enumerate(pseudo_gen_dicts):
        bam_wgths = get_bam_weight(item['bam'], item['weights'])
        restut[i] = get_best_set(bam_wgths['bam'], bam_wgths['weights'])
    return restut


def best_set_worker(pseudo_gen_dicts: list[dict], run_num: int, return_dict: dict,) -> None:
    for key, item in get_milti_bset_set(pseudo_gen_dicts).items():
        idx = (run_num * len(pseudo_gen_dicts)) + key
        return_dict.update({idx: item})


def find_best_sets(pseudo_gen_dicts: list[dict[str, NDArray, list]], num_cor: int = 1) -> dict[dict]:

    """
    Iterate through 2D array of NLLR values finding the best combination for
    each row using the corresponding binary acceptance matrix

    Parameters
    ----------
    overlaps:
    nllr_dat :   NLLR values (weights), 2D Array (N, M)
    num_cor  :   Number of cores Default 1

    Returns
    -------
    results : ndarray
        list of M results, sum of best combination weights from each NLLR set.

    """
    if num_cor < 2:
        print(F"Starting job 1. Calculating {len(pseudo_gen_dicts)} best combinations")
        outputdict = get_milti_bset_set(pseudo_gen_dicts)
    else:
        jobs = []
        manager = Manager()
        outputdict = manager.dict()
        bam_weights = split_list(pseudo_gen_dicts, num_cor)
        for i, bam_wgts in enumerate(bam_weights):
            p = Process(target=best_set_worker, args=(bam_wgts, i, outputdict))
            jobs.append(p)
            p.start()
            print(F"Starting job {i+1}. Calculating {len(bam_wgts)} best combinations")
        for p in jobs:
            p.join()
    return dict(sorted(outputdict.items()))
