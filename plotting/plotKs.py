#!/usr/bin/env python3

import glob, argparse, copy
import matplotlib.pyplot as plt
import scipy.stats
import numpy as np
import matplotlib

matplotlib.use('agg')

def read( which="fake" ):
    pattern = { "fake": "fake*dict", "real": "real?.dict", "realf": "realf*dict",
                "signal": "signal*.dict", "signalf": "signal*f.dict" }
    files = glob.glob( pattern[which] )
    Ks=[]
    for f in files:
        h=open(f,"rt")
        lines=h.read()
        h.close()
        D=eval(lines)
        Ks.append(D[0]["K"])
    return Ks

def plot( opts: dict, outputfile ):
    """ plot the money plot.
    :param opts: dictionary detailing what to plot, e.g { "signals": True, 
          "fastlim": True, "real": True }
    :param outputfile: the filename of outputfile, eg Kvalues.png 
    """
    Ks=read( "fake" )
    Kreal = read ( "real" )
    Ksig = read ( "signal" )
    Ksigf = read ( "realf" )
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
    print ( "maxK", maxKs, fmax, fmax*maxKs )
    arange = np.arange ( fmin*minKs, fmax*maxKs, (maxKs-minKs)/npoints )
    values = kde.evaluate ( arange )
    # label = "KDE of $K_\mathrm{fake}$"
    label = "density, SM hypothesis"
    plt.plot ( arange, values, c="tab:orange", label=label )
    # plt.plot ( arange, values, c="tab:orange" )
    ys = kde.evaluate ( Ks )
    ys = [ x + .001 for x in ys ]
    if opts["fakes"]:
        # label = "$K_\mathrm{fake}$"
        label = "SM hypothesis runs"
        plt.plot ( Ks, ys, "ro", label=label )
        print ( "K(bg)=%.3f, [%.3f,%.3f]" % ( np.mean(Ks), min(Ks), max(Ks) ) )
    if opts["signal"] and len(Ksig)>0:
        ysig = kde.evaluate( Ksig )
        ysig = [ x - .001 for x in ysig ]
        # marker="c*"
        # marker="ro"
        marker_style = dict(color='tab:red', linestyle='', marker='o',
                      markersize=8, fillstyle="none" )
        # label = "$K_\mathrm{signal}$"
        label = "BSM hypothesis"
        plt.plot ( Ksig, ysig, label=label, **marker_style )
        print ( "K(signal)=%.3f, [%.3f,%.3f]" % ( np.mean(Ksig), min(Ksig), max(Ksig) ) )
    if opts["fastlim"]:
        ysigf = kde.evaluate( Ksigf )
        ysigf = [ x - .002 for x in ysigf ]
        label = "K$_\mathrm{signal}^\mathrm{f=0.8}$"
        plt.plot ( Ksigf, ysigf, "m*", ms=8, label=label )
        print ( "K(realf)=%.3f, [%.3f,%.3f]" % ( np.mean(Ksigf), min(Ksigf), max(Ksigf) ) )
    if opts["real"]:
        yreal = kde.evaluate( Kreal )
        yreal = [ x - .001 for x in yreal ]
        # label = "$K_\mathrm{obs}$"
        label = "observed"
        plt.plot ( Kreal, yreal, "g*", ms=8, label=label )
        Krealmean = np.mean(Kreal)
        print ( "K(real)=%.3f, [%.3f,%.3f]" % ( Krealmean, min(Kreal), max(Kreal) ) )
        yrealmean = kde.evaluate ( Krealmean )[0]
        p = kde.integrate_box_1d ( Krealmean, float("inf") )
        pmin = kde.integrate_box_1d ( min(Kreal), float("inf") )
        pmax = kde.integrate_box_1d ( max(Kreal), float("inf") )
        print ( "p(real)=%.3f, [%.3f,%.3f]" % ( p, pmin, pmax ) )
        fromMean = np.arange ( Krealmean, fmax*maxKs+1e-5, (fmax*maxKs-Krealmean)/npoints)
        yFromMean = kde.evaluate ( fromMean )
        label = "$\\bar{\\mathrm{K}}_\mathrm{obs}$"
        plt.plot ( [ Krealmean, Krealmean ], [ yrealmean, 0. ], c="g" )
        # plt.plot ( [ Krealmean, Krealmean ], [ yrealmean, 0. ], c="g", label=label )
    # fromMean = [ Krealmean ] + fromMean + [ 1.1*maxKs ]
    # yFromMean = [ 0] + yFromMean + [0]
    # plt.plot ( fromMean, yFromMean, linewidth=.3, c="tab:orange", label="p", zorder=5 )
        plt.fill_between ( fromMean, yFromMean, 0, linewidth=.3, label="$p$", 
                           facecolor="tab:green", alpha=.5, zorder=-1 )
        plt.title ( "Determination of $p(\\mathrm{global}) \\approx %.2f$" % p, fontsize=16  )
    else:
        plt.title ( "Determination of the Density of $K_\\mathrm{fake}$", fontsize=16 )
        
    plt.ylabel ( "$\\rho(K)$", fontsize=15 )
    plt.xticks ( fontsize=13 )
    plt.yticks ( fontsize=12 )
    plt.xlabel ( "$K$", fontsize=15 )
    plt.legend ( fontsize=14 )
    plt.tight_layout()
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
    argparser.add_argument ( '-o', '--outputfile', help="specify the outputfile",
                             type=str, default="Kvalues.png" )
    args = argparser.parse_args()
    opts = { "signal": args.signals, "fastlim": args.fastlim, "real": args.real,
             "fakes": args.fakes }
    plot( opts, args.outputfile )
