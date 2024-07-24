#!/usr/bin/env python3

from typing import List, Dict
from icecream import ic
import numpy as np
import os

def addPValue ( pvalues : List, v : Dict, verbose : bool , min_expected : float ):
    """ add a pvalue to the list of pvalues """
    if "new_p" in v:
        expectedBG = v["expectedBG"]
        bgError = v["bgError"]
        if verbose:
            s3rd=""
            if "thirdMoment" in v:
                s3rd = ";"+v["thirdMoment"]
            print ( f"{k}:expectedBG={expectedBG}+-{bgError}{s3rd} newObs={v['newObs']} p={v['new_p']}" )
        if expectedBG > min_expected:
            pvalues.append ( v["new_p"] )

def extractPValues( analyses : List, directory : os.PathLike, verbose,
                    min_expected ):
    from ptools.helpers import readDictionaryFile
    import glob
    files = glob.glob ( f"{directory}/*.dict" )
    pvalues = []
    nuniverses = 0
    print ( f"[plotPAcrossMultiverse] found {len(files)} files in '{directory}/'" )
    for f in files:
        D = readDictionaryFile ( f )["data"]
        for k,v in D.items():
            for ana in anas:
                if ":" in ana:
                    if ana == k:
                        addPValue ( pvalues, v, verbose, min_expected )
                elif ana in k:
                    addPValue ( pvalues, v, verbose, min_expected )
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

def runPlotting( anas : List, directory : os.PathLike, outfile : os.PathLike, 
                 verbose : bool, min_expected : float ):
    info = extractPValues( anas, directory, verbose, min_expected )
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
    argparser.add_argument ( '-m', '--min_expected',
            help='minimum expected value to add [0.0]',
            type=float, default=0.0 )
    argparser.add_argument ( '-v', '--verbose',
            help='verbose', action="store_true" )
    args = argparser.parse_args()
    anas = args.analyses.split(",")
    # anas = [ "CMS-EXO-20-004" ]
    runPlotting( anas, args.directory, args.outfile, args.verbose, 
                 args.min_expected )
