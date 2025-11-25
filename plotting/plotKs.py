#!/usr/bin/env python3

import glob, argparse, copy, os
os.environ["DISPLAY"]=""
import matplotlib.pyplot as plt
import scipy.stats
import numpy as np
import matplotlib
from typing import List

matplotlib.use('agg')

def read( which : str ="fake", datadir : str ="./" ) -> List:
    """ read the K values from hiscores dictionary files.
    the type of run is determine from the filename, fake*dict or bg*dict is
    SM-only synthetic, real*dict corresponds to the actual observations, etc

    :returns: a list of K values
    """
    pattern = { "fake": [ "fake*dict", "bg*dict" ], "real": [ "real*.dict" ], 
        "realf": [ "realf*dict" ], "signal": [ "signal?.dict" ], 
        "signalf": [ "signal*f.dict" ] }
    if not which in pattern:
        pattern[which]=f"{which}*.dict"
    files = []
    for wh in pattern[which]:
        files += list ( glob.glob( f"{datadir}/{wh}" ) )
    Ks=[]
    for f in files:
        h=open(f,"rt")
        lines=h.read()
        h.close()
        try:
            D=eval(lines)
            Ks.append(D[0]["K"])
        except (ValueError,SyntaxError) as e:
            print ( f"[plotKs] error when reading {f}: {e}" )
    return Ks


def plot( opts: dict ):
    """ plot the money plot.
    :param opts: dictionary detailing what to plot, e.g { "signals": True,
          "fastlim": True, "real": True }
    :iparam outputfile: the filename of outputfile, eg Kvalues.png
    :iparam datadir: directory of the data dict files, eg ./
    """
    outputfile = opts["outputfile"]
    datadir = opts["datadir"]
    Ks=read ( opts["fakeprefix"], datadir )
    Kreal = read ( "real", datadir )
    Ksig = read ( opts["signalprefix"], datadir )
    Ksigf = read ( "realf", datadir )
    allK = copy.deepcopy ( Ks )
    fmin, fmax, npoints = .3, 1.2, 100
    if opts["real"]:
        allK += Kreal
    if opts["signals"]:
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
    plt.plot ( arange, values, c="tab:orange", label=r"KDE of $K_\mathrm{fake}$" )
    ys = kde.evaluate ( Ks )
    ys = [ x + .001 for x in ys ]
    if opts["fakes"]:
        plt.plot ( Ks, ys, "ro", label=r"$K_\mathrm{fake}$" )
        print ( f"K(bg)={np.mean(Ks):.3f}, [{min(Ks):.3f},{max(Ks):.3f}] {len(Ks)} entries" )
    if opts["signals"] and len(Ksig)>0:
        ysig = kde.evaluate( Ksig )
        ysig = [ x - .001 for x in ysig ]
        # marker="c*"
        # marker="ro"
        marker_style = dict(color='tab:red', linestyle='', marker='o',
                      markersize=8, fillstyle="none" )
        plt.plot ( Ksig, ysig, label=r"$K_\mathrm{signal}$", **marker_style )
        print ( f"K(signal)={np.mean(Ksig):.3f}, [{min(Ksig):.3f},{max(Ksig):.3f}] {len(Ksig)} entries" )
    if opts["fastlim"]:
        ysigf = kde.evaluate( Ksigf )
        ysigf = [ x - .002 for x in ysigf ]
        plt.plot ( Ksigf, ysigf, "m*", ms=8, label=r"K$_\mathrm{signal}^\mathrm{f=0.8}$" )
        print ( f"K(realf)={np.mean(Ksigf):.3f}, [{min(Ksigf):.3f},{max(Ksigf):.3f}] {len(Ksigf)} entries" )
    if opts["real"]:
        yreal = kde.evaluate( Kreal )
        yreal = [ x - .001 for x in yreal ]
        plt.plot ( Kreal, yreal, "g*", ms=8, label=r"$K_\mathrm{obs}$" )
        Krealmean = np.mean(Kreal)
        print ( f"K(real)={Krealmean:.3f}, [{min(Kreal):.3f},{max(Kreal):.3f}] {len(Kreal)} entries" )
        yrealmean = kde.evaluate ( Krealmean )[0]
        p = kde.integrate_box_1d ( Krealmean, float("inf") )
        pmin = kde.integrate_box_1d ( min(Kreal), float("inf") )
        pmax = kde.integrate_box_1d ( max(Kreal), float("inf") )
        print ( f"p(real)={p:.3f}, [{pmin:.3f},{pmax:.3f}]" )
        fromMean = np.arange ( Krealmean, fmax*maxKs+1e-5, (fmax*maxKs-Krealmean)/npoints)
        yFromMean = kde.evaluate ( fromMean )
        plt.plot ( [ Krealmean, Krealmean ], [ yrealmean, 0. ], c="g", label=r"$\bar{\mathrm{K}}_\mathrm{obs}$" )
    # fromMean = [ Krealmean ] + fromMean + [ 1.1*maxKs ]
    # yFromMean = [ 0] + yFromMean + [0]
    # plt.plot ( fromMean, yFromMean, linewidth=.3, c="tab:orange", label="p", zorder=5 )
        plt.fill_between ( fromMean, yFromMean, 0, linewidth=.3, label="$p$",
                           facecolor="tab:green", alpha=.5, zorder=-1 )
        plt.title ( r"Determination of $p(\mathrm{global}) \approx %.2f$" % p  )
    else:
        plt.title ( r"Determination of the Density of $K_\mathrm{fake}$" )

    plt.ylabel ( r"$\rho(K)$" )
    plt.xlabel ( "$K$" )
    plt.legend ()
    print ( f"[plotKs] saving to {outputfile}." )
    from smodels_utils.helper.various import pngMetaInfo
    metadata = pngMetaInfo()
    plt.savefig ( outputfile, metadata = metadata )

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
    """
    opts = { "signal": args.signals, "fastlim": args.fastlim, "real": args.real,
             "fakes": args.fakes, "fakeprefix": args.fakeprefix,
             "signalprefix": args.signalprefix }
    """
    opts = vars ( args )
    from helpers.various import viewImage
    plot( opts )
    viewImage ( opts["outputfile"] )
