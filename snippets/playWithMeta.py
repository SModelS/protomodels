#!/usr/bin/env python3

from scipy import stats
import random
import pickle
import scipy.stats
import numpy
import numpy as np

def drawNuisance ( mu = 0., sigma = 1., lognormal=True, size=1000 ):
    """ draw from the nuisance model.
    :param mu: the expecation value of the distribution
    :param sigma: the standard deviation of the distribution
    :returns: a single fake observation
    """
    ## if lognormal is not selected, or mu is very small
    ## (as the lognormal cannot really produce mu=0)
    if not lognormal or abs(mu)<(sigma/4.):
        ret = stats.norm.rvs ( mu, sigma, size=size )
        ret = ret[ret>0.]
    else:
        loc = mu**2 / numpy.sqrt ( mu**2 + sigma**2 )
        stderr = numpy.sqrt ( numpy.log ( 1 + sigma**2 / mu**2 ) )
        ret = stats.lognorm.rvs ( s=stderr, scale=loc, size=size )
    fakeobs = scipy.stats.poisson.rvs ( ret )
    return fakeobs

def computePWithToys ( obs, bg, bgerr, sigN, lognormal, n ):
    """ compute p value, for now we assume Gaussanity """
    fakes = []
    bigger = 0
    central = bg
    signalmodel = False
    if signalmodel and sigN != None:
        central = bg + sigN
    if lognormal and central > ( bgerr / 4. ):
        loc = central**2 / np.sqrt ( central**2 + bgerr**2 )
        stderr = np.sqrt ( np.log ( 1 + bgerr**2 / central**2 ) )
        lmbda = scipy.stats.lognorm.rvs ( s=[stderr]*n, scale=[loc]*n )
    else: ## Gauss
        lmbda = scipy.stats.norm.rvs ( loc=[central]*n, scale=[bgerr]*n )
        lmbda = lmbda[lmbda>0.]
    fakeobs = scipy.stats.poisson.rvs ( lmbda )
    return sum(fakeobs>obs) / len(fakeobs)

import matplotlib.pyplot as plt

def read():
    f=open("ps.pcl","rb")
    ps = pickle.load(f)
    f.close()
    return ps
    
def write():
    ps = []
    logn = True# False
    for i in range(200):
        if i % 10 == 0:
            print ( ".", end="", flush=True )
        mu = random.uniform ( 1., 30. )
        sigma = -1
        while sigma < 0.:
            sigma = random.gauss( mu, mu )
        ss = drawNuisance ( mu, sigma, logn, 100 )
        for x in ss:
            p = computePWithToys ( x, mu, sigma, 0, logn, 10000 )
            D = { "p": p, "x": x, "mu": mu, "sigma": sigma }
            ps += [ D ]
    f=open("ps.pcl","wb")
    pickle.dump(ps,f)
    f.close()
    return ps

def main():
    ps = write()
    # ps = read()
    pvs =  [ x["p"] for x in ps ] 
    mn,st = np.mean(pvs),np.std(pvs)
    print ( "ps=",mn,"+-",st )
    plt.hist ( pvs )
    plt.title ( "p-values: %.2f+-%.2f" % ( mn,st)  )
    plt.xlabel ( "p-values" )
    plt.savefig ( "ps.png" )

main()
