#!/usr/bin/env python3

import scipy.stats
import matplotlib.pyplot as plt
import numpy as np
import subprocess

def m ( x : float , n : int = 10 ):
    # ret = n*scipy.stats.norm.cdf(x)**(n-1)*scipy.stats.norm.pdf(x)
    ret = n*scipy.stats.t.cdf(x,df=3,loc=-1.18,scale=4.73/5.)**(n-1)*scipy.stats.t.pdf(x,df=3,loc=-1.18,scale=4.73/5.)
    return ret

def getExpected ( xs, ys ):
    S = 0.
    w = 0.
    for x,y in zip ( xs, ys ):
        S+= x*y
        w+= y
    return S / w

def getVariance ( xs, ys, expected ):
    V = 0.
    for x,y in zip ( xs, ys ):
        V+= y*(x-expected)**2
    return V#  / (len(xs)-1)

def plot ( ):
    # xs = np.arange(-3,6,.03)
    xs = np.arange(-7,8,.1)
    ns = [ 1, 2, 5, 20 ]
    dicts = { n: [] for n in ns }
    for x in xs:
        for n in ns:
            dicts[n].append ( m (x, n ) )
    for n in ns:
        e = getExpected ( xs, dicts[n] )
        v = np.sqrt ( getVariance ( xs, dicts[n], e ) )
        plt.plot(xs,dicts[n],label = f"max({n} SN): {e:.2f}+-{v:.2f}" )
    plt.legend()
    plt.savefig("penalty.png")
    o = subprocess.getoutput ( "timg penalty.png" )
    print ( o )
        

if __name__ == "__main__":
    plot ()
