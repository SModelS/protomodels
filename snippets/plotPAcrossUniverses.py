#!/usr/bin/env python3

from typing import List

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
    import subprocess
    o = subprocess.getoutput ( "timg pvalues.png" )
    print ( o )
    # import sys, IPython; IPython.embed( colors = "neutral" ); sys.exit()

def runPlotting( anas : List ):
    info = extractPValues( anas )
    plotPValues ( info, anas )

if __name__ == "__main__":
    # anas = [ "ATLAS-SUSY-2019-09" ]
    anas = [ "CMS-SUS-21-002" ]
    runPlotting( anas )
