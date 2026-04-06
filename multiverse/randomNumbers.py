#!/usr/bin/env python3

""" simple numba methods to make random number generation ultra fast
round model. """

import numpy as np
from numba import njit, prange

# NORMAL
@njit(parallel=True)
def normal_parallel(n, mu, sigma):
    out = np.empty(n)
    for i in prange(n):
        out[i] = mu + sigma * np.random.normal()
    return out


# LOGNORMAL
@njit(parallel=True)
def lognorm_parallel(n, s, scale=1.0, loc=0.0):
    out = np.empty(n)
    for i in prange(n):
        z = np.random.normal()
        out[i] = loc + scale * np.exp(s * z)
    return out

# POISSON
@njit(parallel=True)
def poisson_parallel(n, lam):
    out = np.empty(n, dtype=np.int64)
    for i in prange(n):
        out[i] = np.random.poisson(lam)
    return out

if __name__ == "__main__":
    import time
    import scipy.stats
    x = lognorm_parallel ( 1, 1., 1., 0. )
    n = 10_000_000
    t0 = time.time()
    x = lognorm_parallel ( n, s=1., scale=1., loc=0. )
    t1 = time.time()
    print ( np.mean(x), np.std(x), t1-t0 )
    x = scipy.stats.lognorm.rvs ( s=[1.]*n, scale=[1.]*n, loc=[0.]*n )
    t2 = time.time()
    print ( np.mean(x), np.std(x), t2-t1 )
