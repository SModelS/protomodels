#!/usr/bin/env python3

from typing import List
from icecream import ic
import numpy as np
import os

def extractPValues( analyses : List, directory : os.PathLike, verbose ):
    from ptools.helpers import readDictionaryFile
    import glob
    files = glob.glob ( f"{directory}/*.dict" )
    pvalues = []
    nuniverses = 0
    print ( f"[plotPAcrossMultiverse] found {len(files)} files in '{directory}/'" )
    for f in files:
        hasEntry = False
        D = readDictionaryFile ( f )["data"]
        for k,v in D.items():
            for ana in anas:
                if ana == k:
                    if "new_p" in v:
                        expectedBG = v["expectedBG"]
                        bgError = v["bgError"]
                        if verbose:
                            print ( f"{k}:expectedBG={expectedBG}+-{bgError} newObs={v['newObs']} p={v['new_p']}" )
                        if True: # expectedBG + 2*bgError > 5.:
                            pvalues.append ( v["new_p"] )
                        hasEntry = True
        nuniverses += 1
    # print ( pvalues )
    return { "pvalues": pvalues, "nuniverses": nuniverses }

def plotPValues( info, anas, outfile ):
    pvalues = info["pvalues"]
    nuniverses = info["nuniverses"]
    from matplotlib import pyplot as plt
    fig, ax = plt.subplots()
    plt.hist ( pvalues, bins=np.arange(0.0,1.01,0.1) )
    plt.xlabel ( "p-values" )
    plt.ylabel ( "# SRs" )
    title = ", ".join ( anas )
    # title += f" [{info['nuniverses']} universes]" 
    plt.title ( title )
    plt.text ( -.1, -.1, f"{nuniverses} universes", transform=ax.transAxes )
    plt.savefig ( outfile )
    import shutil
    if shutil.which ("timg") is not None:
        import subprocess
        o = subprocess.getoutput ( f"timg {outfile}" )
        print ( o )
    # import sys, IPython; IPython.embed( colors = "neutral" ); sys.exit()

def runPlotting( anas : List, directory : os.PathLike, outfile : os.PathLike, verbose : bool ):
    info = extractPValues( anas, directory, verbose )
    plotPValues ( info, anas, outfile )

if __name__ == "__main__":
    import argparse
    argparser = argparse.ArgumentParser( description='simple script to plot p-values across universes' )
    argparser.add_argument ( '-a', '--analyses',
            help='analyses to plot p-values for, comma separated ["ATLAS-SUSY-2019-09"]',
            type=str, default="ATLAS-SUSY-2019-09" )
    argparser.add_argument ( '-d', '--directory',
            help='directory to look for dict files ["dicts"]',
            type=str, default="dicts" )
    argparser.add_argument ( '-o', '--outfile',
            help='output file name ["pvalues.png"]',
            type=str, default="pvalues.png" )
    argparser.add_argument ( '-v', '--verbose',
            help='verbose', action="store_true" )
    args = argparser.parse_args()
    anas = args.analyses.split(",")
    # anas = [ "CMS-EXO-20-004" ]
    runPlotting( anas, args.directory, args.outfile, args.verbose )
