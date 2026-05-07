#!/usr/bin/env python3

""" simple numba methods to make random number generation ultra fast
round model. """

__all__ = [ "fast_rvs" ]

import numpy as np
from numba import njit, prange, set_num_threads

# set_num_threads(4)
parallel = True

# NORMAL
@njit(parallel=parallel)
def normal_parallel(n, mu, sigma):
    out = np.empty(n)
    for i in prange(n):
        out[i] = mu + sigma * np.random.normal()
    return out

# LOGNORMAL
@njit(parallel=parallel)
def lognorm_parallel(n, s, scale=1.0, loc=0.0):
    out = np.empty(n)
    for i in prange(n):
        z = np.random.normal()
        out[i] = loc + scale * np.exp(s * z)
    return out

# POISSON
@njit(parallel=parallel)
def poisson_parallel(n, mu ):
    out = np.empty(n, dtype=np.int64)
    for i in prange(n):
        lam = mu[i]
        if lam < 0.:
            lam = 0.0
        out[i] = np.random.poisson( lam )
    return out

def fast_rvs(dist, n, **p):
    if False: # n < 100000:
        import scipy.stats
        if dist == "normal":
            return scipy.stats.norm.rvs ( [p.get("loc",0.0)]*n,
                [p.get("scale",1.0)]*n )
        if dist == "lognorm":
            return scipy.stats.lognorm.rvs ( s=[p.get("s",0.0)]*n,
                scale=[p.get("scale",1.0)]*n,loc=[p.get("loc",0.0)]*n )
        if dist == "poisson":
            return scipy.stats.poisson.rvs ( [p.get("mu",0.0)]*n )
            
    if dist == "normal":
            return normal_parallel(n, p.get("loc",0.0), p.get("scale",1.0))
    if dist == "lognorm":
        return lognorm_parallel(n, p["s"], p.get("scale",1.0), p.get("loc",0.0))
    if dist == "poisson":
        arg = p["mu"]
        if isinstance(arg,float) or isinstance(arg,np.float64):
            arg = np.array( [ arg ]*n )
        return poisson_parallel(n, arg )

if __name__ == "__main__":
    import time
    import scipy.stats
    x = fast_rvs ( "poisson", 1, mu=[1.] )
    n = 100_000_000
    # n = 100000
    t0 = time.time()
    x = fast_rvs ( "poisson", n, mu=np.array([1.]*n) )
    t1 = time.time()
    print ( "numba", x, np.mean(x), np.std(x), t1-t0 )
    """
    x = scipy.stats.lognorm.rvs ( s=[1.]*n, scale=[1.]*n, loc=[0.]*n )
    t2 = time.time()
    print ( "scipy", np.mean(x), np.std(x), t2-t1 )
    """
