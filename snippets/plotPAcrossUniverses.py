#!/usr/bin/env python3

from typing import List
from icecream import ic

def extractPValues( analyses : List ):
    from ptools.helpers import readDictionaryFile
    import glob
    files = glob.glob ( "dicts/*.dict" )
    pvalues = []
    nuniverses = 0
    for f in files:
        hasEntry = False
        D = readDictionaryFile ( f )["data"]
        for k,v in D.items():
            for ana in anas:
                if ana in k:
                    if "new_p" in v:
                        expectedBG = v["expectedBG"]
                        bgError = v["bgError"]
                        if True: # expectedBG + 2*bgError > 5.:
                        # ic ( f"{k}: p={v['new_p']}" )
                            pvalues.append ( v["new_p"] )
                        hasEntry = True
        nuniverses += 1
    # print ( pvalues )
    return { "pvalues": pvalues, "nuniverses": nuniverses }

def plotPValues( info, anas ):
    pvalues = info["pvalues"]
    from matplotlib import pyplot as plt
    plt.hist ( pvalues )
    plt.xlabel ( "p-values" )
    title = ", ".join ( anas )
    title += f" [{info['nuniverses']} universes]" 
    plt.title ( title )
    plt.savefig ( "pvalues.png" )
    import shutil
    if shutil.which ("timg") is not None:
        import subprocess
        o = subprocess.getoutput ( "timg pvalues.png" )
        print ( o )
    # import sys, IPython; IPython.embed( colors = "neutral" ); sys.exit()

def runPlotting( anas : List ):
    info = extractPValues( anas )
    plotPValues ( info, anas )

if __name__ == "__main__":
    import argparse
    argparser = argparse.ArgumentParser( description='simple script to plot p-values across universes' )
    argparser.add_argument ( '-a', '--analyses',
            help='analyses to plot p-values for, comma separated ["ATLAS-SUSY-2019-09"]',
            type=str, default="ATLAS-SUSY-2019-09" )
    args = argparser.parse_args()
    anas = args.analyses.split(",")
    # anas = [ "CMS-EXO-20-004" ]
    runPlotting( anas )
