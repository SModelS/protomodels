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
    ns = range(1,200)
    plot = [ 1, 2, 5, 20 ]
    dicts = { n: [] for n in ns }
    for x in xs:
        for n in ns:
            dicts[n].append ( m (x, n ) )
    for n in ns:
        e = getExpected ( xs, dicts[n] )
        v = np.sqrt ( getVariance ( xs, dicts[n], e ) )
        print ( f"n {n} e {e:.2f}+-{v:.2f}" )
    for n in plot:
        plt.plot(xs,dicts[n],label = f"max({n} SR): {e:.2f}+-{v:.2f}" )
    plt.legend()
    plt.savefig("penalty.png")
    o = subprocess.getoutput ( "timg penalty.png" )
    print ( o )
        
def plotTrend ():
    with open ( "penalties", "rt" ) as f:
        lines = f.readlines()
        f.close()
    xs, ys = [] , []
    for line in lines:
        tokens = line.split("," )
        xs.append ( int(tokens[0]) )
        ys.append ( float(tokens[1]) )
    ys = [ y - ys[0] for y in ys ] # subtract first!
    #print ( "xs", xs[:10] )
    #print ( "ys", ys[:10] )
    #sigmas = [ np.sqrt ( y ) + 1e-6 for y in ys ]
    plt.scatter ( xs, ys, label = "data (approx. Jamie's Fig 8.7)" )
    from scipy.optimize import curve_fit
    popt, pcov = curve_fit(lambda t, a: a*np.log(t), xs, ys ) # , sigma=sigmas, absolute_sigma=True)
    #print ( "popt", popt, pcov )
    a = popt[0]
    # b = popt[1]
    #print ( "a", a )
    print ( f"{a} * log ( x )" )
    preds = [ a * np.log(t) for t in xs ]
    simplelog = [ np.log(t) for t in xs ]
    plt.plot ( xs, preds, c="orange", label=f"{a:.3f}*ln(n)" )
    plt.plot ( xs, simplelog, c="red", label="ln(n)" )
    plt.legend()
    plt.title ( "first rough sketch for an <n>-SR penalty term" )
    plt.xlabel ( "n" )
    plt.ylabel ( "expected weight" )
    plt.savefig ( "trend.png" )
    o = subprocess.getoutput ( "timg trend.png" )
    print ( o )
    # import sys, IPython; IPython.embed( colors = "neutral" ); sys.exit()

if __name__ == "__main__":
    plotTrend ()
    # plot ()
