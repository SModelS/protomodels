#!/usr/bin/env python3

""" plot the evolution of Z for every data id as dm increases.

    the input is a dictionary as produced by logZs.py, i.e.
    { dm: { "SR_2018_1": { "nll_sm": .., "nll_bsm": .., "Z": .. }, ... }, ... }
"""

import os, sys, argparse
import numpy as np
import matplotlib.pyplot as plt


def readDict ( filename : str ) -> dict:
    """ read in the (python) dictionary at filename """
    if not os.path.exists ( filename ):
        print ( f"[plotZs] {filename} does not exist." )
        sys.exit()
    with open ( filename, "rt" ) as f:
        return eval ( f.read() )


def collectDataIds ( data : dict ) -> list:
    """ get the sorted list of all data ids that appear anywhere in data """
    ret = set()
    for dms in data.values():
        ret.update ( dms.keys() )
    return sorted ( ret )


def plot ( data : dict, outputfile : str = "Zs.png",
           dataIds : list = None, title : str = None ):
    """ plot Z vs dm, one curve per data id """
    dms = sorted ( data.keys() )
    if dataIds == None:
        dataIds = collectDataIds ( data )

    fig, ax = plt.subplots ( figsize = ( 13, 8 ) )
    cmap = plt.get_cmap ( "turbo" )
    nc = max ( len(dataIds)-1, 1 )
    for i, dId in enumerate ( dataIds ):
        xs = [ dm for dm in dms if dId in data[dm] ]
        ys = [ data[dm][dId]["Z"] for dm in dms if dId in data[dm] ]
        if len(xs)==0:
            continue
        color = cmap ( i/nc )
        ax.plot ( xs, ys, marker = "o", markersize = 3,
                  color = color, label = dId, alpha = 0.8 )
        for x0, x1, y0, y1 in zip ( xs[:-1], xs[1:], ys[:-1], ys[1:] ):
            if abs ( y1 - y0 ) > 1.:
                ax.plot ( [ x0, x1 ], [ y0, y1 ], color = color,
                          linewidth = 5, alpha = 0.8,
                          solid_capstyle = "round", zorder = 3 )
        ax.text ( xs[-1], ys[-1], f"  {dId}", color = color,
                  fontsize = 6, va = "center", ha = "left" )

    ax.axhline ( 0., color = "black", linewidth = 1, linestyle = "--" )
    ax.set_xlabel ( r"$\Delta m$ [GeV]" )
    ax.set_ylabel ( r"$Z$" )
    ax.set_xlim ( left = dms[0], right = dms[0] + ( dms[-1] - dms[0] ) * 1.4 )
    if title == None:
        title = "evolution of Z with dm"
    ax.set_title ( title )
    if len(dataIds) <= 20:
        ax.legend ( fontsize = 7, ncol = 2 )
    fig.tight_layout()
    fig.savefig ( outputfile, dpi = 200 )
    print ( f"[plotZs] plotting {len(dataIds)} data ids -> {outputfile}" )
    return fig, ax


if __name__ == "__main__":
    argparser = argparse.ArgumentParser(
        description = 'plot the evolution of Z with dm' )
    argparser.add_argument ( '-i', '--inputfile',
        help = 'input file [Zs.dict]', type = str, default = "Zs.dict" )
    argparser.add_argument ( '-o', '--outputfile',
        help = 'output file [Zs.png]', type = str, default = "Zs.png" )
    argparser.add_argument ( '-D', '--dataIds', nargs = '*',
        help = 'restrict to these data ids [all]', default = None )
    argparser.add_argument ( '-T', '--title',
        help = 'a title', type = str, default = None )
    args = argparser.parse_args()
    data = readDict ( args.inputfile )
    plot ( data, args.outputfile, args.dataIds, args.title )
