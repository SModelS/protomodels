#!/usr/bin/env python3

import glob, argparse, copy, os
os.environ["DISPLAY"]=""
import matplotlib.pyplot as plt
import scipy.stats
import numpy as np
import matplotlib

matplotlib.use('agg')

def read( which="fake", datadir="./" ):
    pattern = { "fake": "fake*dict", "real": "real*.dict", "realf": "realf*dict",
                "signal": "signal?.dict", "signalf": "signal*f.dict" }
    if not which in pattern:
        pattern[which]=which+"*.dict"
    files = glob.glob( datadir + pattern[which] )
    Ks=[]
    for f in files:
        h=open(f,"rt")
        lines=h.read()
        h.close()
        D=eval(lines)
        Ks.append(D[0]["K"])
    return Ks

def plot( opts: dict, outputfile, datadir ):
    """ plot the money plot.
    :param opts: dictionary detailing what to plot, e.g { "signals": True,
          "fastlim": True, "real": True }
    :param outputfile: the filename of outputfile, eg Kvalues.png
    :param datadir: directory of the data dict files, eg ./
    """
    Ks=read ( opts["fakeprefix"], datadir )
    Kreal = read ( "real", datadir )
    Ksig = read ( opts["signalprefix"], datadir )
    Ksigf = read ( "realf", datadir )
    allK = copy.deepcopy ( Ks )
    fmin, fmax, npoints = .3, 1.2, 100
    if opts["real"]:
        allK += Kreal
    if opts["signal"]:
        allK += Ksig
        fmax = 1.1
    if opts["fastlim"]:
        allK += Ksigf
        fmax = 1.01
    kde = scipy.stats.gaussian_kde ( Ks )
    minKs = min ( allK )
    maxKs = max ( allK )
    # print ( "maxK", maxKs, "fmax", fmax, "product", fmax*maxKs )
    arange = np.arange ( fmin*minKs, fmax*maxKs, (maxKs-minKs)/npoints )
    values = kde.evaluate ( arange )
    plt.plot ( arange, values, c="tab:orange", label="KDE of $K_\mathrm{fake}$" )
    ys = kde.evaluate ( Ks )
    ys = [ x + .001 for x in ys ]
    if opts["fakes"]:
        plt.plot ( Ks, ys, "ro", label="$K_\mathrm{fake}$" )
        print ( "K(bg)=%.3f, [%.3f,%.3f] %d entries" % \
                ( np.mean(Ks), min(Ks), max(Ks), len(Ks) ) )
    if opts["signal"] and len(Ksig)>0:
        ysig = kde.evaluate( Ksig )
        ysig = [ x - .001 for x in ysig ]
        # marker="c*"
        # marker="ro"
        marker_style = dict(color='tab:red', linestyle='', marker='o',
                      markersize=8, fillstyle="none" )
        plt.plot ( Ksig, ysig, label="$K_\mathrm{signal}$", **marker_style )
        print ( "K(signal)=%.3f, [%.3f,%.3f] %d entries" % \
                ( np.mean(Ksig), min(Ksig), max(Ksig), len(Ksig) ) )
    if opts["fastlim"]:
        ysigf = kde.evaluate( Ksigf )
        ysigf = [ x - .002 for x in ysigf ]
        plt.plot ( Ksigf, ysigf, "m*", ms=8, label="K$_\mathrm{signal}^\mathrm{f=0.8}$" )
        print ( "K(realf)=%.3f, [%.3f,%.3f] %d entries" % \
                ( np.mean(Ksigf), min(Ksigf), max(Ksigf), len(Ksigf) ) )
    if opts["real"]:
        yreal = kde.evaluate( Kreal )
        yreal = [ x - .001 for x in yreal ]
        plt.plot ( Kreal, yreal, "g*", ms=8, label="$K_\mathrm{obs}$" )
        Krealmean = np.mean(Kreal)
        print ( "K(real)=%.3f, [%.3f,%.3f] %d entries" % \
                ( Krealmean, min(Kreal), max(Kreal), len(Kreal) ) )
        yrealmean = kde.evaluate ( Krealmean )[0]
        p = kde.integrate_box_1d ( Krealmean, float("inf") )
        pmin = kde.integrate_box_1d ( min(Kreal), float("inf") )
        pmax = kde.integrate_box_1d ( max(Kreal), float("inf") )
        print ( "p(real)=%.3f, [%.3f,%.3f]" % ( p, pmin, pmax ) )
        fromMean = np.arange ( Krealmean, fmax*maxKs+1e-5, (fmax*maxKs-Krealmean)/npoints)
        yFromMean = kde.evaluate ( fromMean )
        plt.plot ( [ Krealmean, Krealmean ], [ yrealmean, 0. ], c="g", label="$\\bar{\\mathrm{K}}_\mathrm{obs}$" )
    # fromMean = [ Krealmean ] + fromMean + [ 1.1*maxKs ]
    # yFromMean = [ 0] + yFromMean + [0]
    # plt.plot ( fromMean, yFromMean, linewidth=.3, c="tab:orange", label="p", zorder=5 )
        plt.fill_between ( fromMean, yFromMean, 0, linewidth=.3, label="$p$",
                           facecolor="tab:green", alpha=.5, zorder=-1 )
        plt.title ( "Determination of $p(\\mathrm{global}) \\approx %.2f$" % p  )
    else:
        plt.title ( "Determination of the Density of $K_\\mathrm{fake}$" )

    plt.ylabel ( "$\\rho(K)$" )
    plt.xlabel ( "$K$" )
    plt.legend ()
    print ( f"[plotKs] saving to {outputfile}." )
    plt.savefig ( outputfile )

if __name__ == "__main__":
    argparser = argparse.ArgumentParser( description="plot the money plots" )
    argparser.add_argument ( '-s', '--signals', help="add the fake signals",
                             action="store_true" )
    argparser.add_argument ( '-b', '--fakes',
                             help="add points for the fake backgrounds",
                             action="store_true" )
    argparser.add_argument ( '-f', '--fastlim', help="add the fastlim real runs",
                             action="store_true" )
    argparser.add_argument ( '-r', '--real', help="add the real Ks",
                             action="store_true" )
    argparser.add_argument ( '-D', '--datadir',
                             help="specify the directory of the dict files [./]",
                             type=str, default="./" )
    argparser.add_argument ( '-o', '--outputfile', help="specify the outputfile [Kvalues.png]",
                             type=str, default="Kvalues.png" )
    argparser.add_argument ( '--fakeprefix', help="specify the prefix for the fakes [fake]",
                             type=str, default="fake" )
    argparser.add_argument ( '--signalprefix', 
                             help="specify the prefix for the signal [signal]",
                             type=str, default="signal" )
    args = argparser.parse_args()
    opts = { "signal": args.signals, "fastlim": args.fastlim, "real": args.real,
             "fakes": args.fakes, "fakeprefix": args.fakeprefix,
             "signalprefix": args.signalprefix }
    plot( opts, args.outputfile, args.datadir )
