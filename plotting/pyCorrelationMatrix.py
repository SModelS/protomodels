#!/usr/bin/env python3

from __future__ import print_function
import sys, os, time, math
sys.path.insert(0,"../")
from smodels.experiment.databaseObj import Database
from smodels.tools.smodelsLogging import setLogLevel
from smodels.tools.physicsUnits import TeV
from smodels.tools.colors import colors
from smodels_utils.helper.various import hasLLHD
from tester import analysisCombiner
import IPython

def sortBySqrts ( results, sqrts ):
    ret = []
    for res in results:
        if abs (res.globalInfo.sqrts.asNumber(TeV) - sqrts ) < 0.1:
            ret.append ( res )
    return ret

def noFastlim ( results ):
    ret = []
    for res in results:
        if hasattr ( res.globalInfo, "contact" ) and "fastlim" in res.globalInfo.contact:
            continue
        ret.append ( res )
    return ret

def sortOutDupes ( results ):
    """ If an analysis id appears more than once in the list,
    keep only the one with likelihoods. """
    isIn = set() ## mark result as being "in"
    ids = set() ## the analysis ids that are already in isIn
    ids_withoutLLHD = set()
    ret_withoutLLHD = {}
    for res in results:
        force_override = False
        ID = res.globalInfo.id
        if "-agg" in ID:
            ID = ID.replace("-agg","")
            force_override = True
        if ID in ids and not force_override: ## already in
            continue
        hasllhd = hasLLHD ( res )
        if hasllhd and not ID in ids:
            ## not in and should be in: add!
            ids.add ( ID )
            isIn.add ( res.globalInfo.path )
            continue
        if not ID in ids and not hasllhd:
            ## not in but shouldnt be in: add to waiting list
            ids_withoutLLHD.add ( ID )
            ret_withoutLLHD[ID]= res.globalInfo.path
    for i in ids_withoutLLHD:
        if not i in ids: ## nothing with llhd is in, so add!
            isIn.add ( ret_withoutLLHD[i] )
            ids.add ( i )
    ## now sort them like in the original container!
    ret = []
    for res in results:
        if res.globalInfo.path in isIn:
            ret.append ( res )
    return ret

def draw( args : dict ):
    """
    draw the correlation matrix
    :param args: dictionary of args
           triangular: if True, then only plot the upper triangle of this
                         symmetrical matrix
           miscol: color to use when likelihood is missing
           diagcol: color to use for diagonal
           experiment: draw only for specific experiment ("CMS", "ATLAS", "all" )
           database: path to database
           sqrts: draw only for specific sqrts ( "8", "13", "all" )
           drawtimestamp: if true, put a timestamp on plot
           outputfile: file name of output file (matrix.png)
           nofastlim: if True, discard fastlim results
    """
    cols = [ "red", "white", "green", args["miscol"], args["diagcol"] ]

    sqrtses = [ 8, 13 ]
    if args["sqrts"] not in [ "all" ]:
        sqrtses = [ int(args["sqrts"]) ]

    colors.on = True
    setLogLevel ( "debug" )

    # dir = "/home/walten/git/smodels-database/"
    dbdir = args["database"]
    d=Database( dbdir, discard_zeroes = True )
    print(d)
    analysisIds = [ "all" ]
    exps = [ "CMS", "ATLAS" ]
    if args["experiment"] in [ "CMS", "ATLAS" ]:
        analysisIds = [ args["experiment"]+"*" ]
        exps = [ args["experiment"] ]
    results = d.getExpResults( analysisIDs = analysisIds )
    if args["nofastlim"]:
        results = noFastlim ( results )
    results = sortOutDupes ( results )
    if args["sqrts"] in [ "8", "13" ]:
        results = sortBySqrts ( results, int(args["sqrts"]) )

    #results.sort()
    nres = len ( results )

    from matplotlib import pyplot as plt
    import matplotlib
    matplotlib.rc('xtick', labelsize=14, labelcolor = "gray" )
    matplotlib.rc('ytick', labelsize=14, labelcolor = "gray" )

    bins= { "CMS": { 8: [999,0], 13:[999,0] },
            "ATLAS": { 8: [999,0], 13: [999,0] } }

    n = len(results )
    import numpy as np
    h = np.array([[0.]*n]*n)
    labels = []
    for x,e in enumerate(results):
        label = e.globalInfo.id
        hasLikelihood = hasLLHD ( e )
        ana = analysisCombiner.getExperimentName ( e.globalInfo )
        #if not hasLikelihood:
        #    print ( "no likelihood: %s" % label )
        sqrts = int(e.globalInfo.sqrts.asNumber(TeV))
        color = "cyan"
        ymax=0
        if ana == "ATLAS":
            color = "blue"
        if sqrts > 10.:
            color = "dark"+color
        if x < bins[ana][sqrts][0]:
            bins[ana][sqrts][0]=x
        if x > bins[ana][sqrts][1]:
            bins[ana][sqrts][1]=x
            ymax=x
        color = "gray"
        if len(exps)==1 and len(sqrtses)==1:
            label = label.replace("CMS-","").replace("ATLAS-","").replace("-agg","")
        labels.append ( label )
        # label = "#color[%d]{%s}" % (color, label )
        for y,f in enumerate(results):
            if args["triangular"] and y>x:
                continue
            isUn = analysisCombiner.canCombine ( e.globalInfo, f.globalInfo, 
                    args["strategy"] )
            # isUn = e.isUncorrelatedWith ( f )
            if isUn:
                h[n-x-1][y]= 1.
            else:
                h[n-x-1][y]= -1
            if not hasLikelihood or not hasLLHD ( f ): ## has no llhd? cannot be combined
                h[n-x-1][y] = 2.
            if y==x:
                h[n-x-1][y] = 3.

    c = [ "red", "b", "b", "b", "limegreen", "b", "white", "grey" ]
    v = np.arange(0.,1.00001,1. / (len(c)-1) )
    l = list(zip(v,c))
    from  matplotlib.colors import LinearSegmentedColormap
    cmap=LinearSegmentedColormap.from_list('rg',l, N=len(c) )
    plt.matshow ( h, aspect = "equal", origin = "lower", cmap = cmap,
                  vmin = -1, vmax = 3. )
    plt.grid ( visible = False )
    plt.xticks ( rotation=90 )
    fig = plt.gcf()
    fig.set_size_inches(30, 30)
    ax = plt.gca()
    ax.set_xticks ( range(len(labels)) )
    ax.set_xticklabels( labels )
    ax.set_yticks ( range(len(labels)) )
    labels.reverse()
    ax.set_yticklabels( labels ) ## need to invert
    bins, xbins, lines = {}, {}, []
    if len(exps)==1 and len(sqrtses)==1:
        plt.text ( .45, .95, "%s, %d TeV" % ( exps[0], sqrtses[0] ),
                   transform = fig.transFigure )
    if args["drawtimestamp"]:
        plt.text ( .01, .01, "plot produced %s from database v%s" % \
                   ( time.strftime("%h %d %Y" ), d.databaseVersion ), 
                   c="grey", transform = fig.transFigure, fontsize=24 )
    # ROOT.c1.Print("matrix_%s.pdf" % strategy )
    outputfile = args["outputfile"]
    if "@M" in outputfile:
        modifiers = ""
        if len(exps)==1:
            modifiers += exps[0]
        if len(sqrtses)==1:
            modifiers += str(sqrtses[0])
        outputfile = outputfile.replace("@M",modifiers)
    print ( "Plotting to %s" % outputfile )
    plt.savefig ( outputfile, dpi=300 )
    return outputfile

def show ( outputfile ):
    import subprocess
    cmd = f"timg {outputfile}"
    o = subprocess.getoutput ( cmd )
    print ( o )

if __name__ == "__main__":
    import argparse
    argparser = argparse.ArgumentParser(description="correlation/combination matrix plotter")
    argparser.add_argument ( '-S', '--strategy', nargs='?',
            help='combination strategy [aggressive]', type=str, default='aggressive' )
    argparser.add_argument ( '-d', '--database', nargs='?',
            help='path to database [../../smodels-database]',
            type=str, default='../../smodels-database' )
    argparser.add_argument ( '-e', '--experiment', nargs='?',
            help='plot only specific experiment CMS,ATLAS,all [all]',
            type=str, default='all' )
    argparser.add_argument ( '-s', '--sqrts', nargs='?',
            help='plot only specific sqrts 8,13,all [all]',
            type=str, default='all' )
    argparser.add_argument ( '-o', '--outputfile', nargs='?',
            help='outputfile (@M gets replaced by [experiment][sqrts]) [matrix@M.png]',
            type=str, default='matrix@M.png' )
    argparser.add_argument ( '-t', '--triangular',
            help='plot as lower triangle matrix?',
            action="store_true" )
    argparser.add_argument ( '-n', '--nofastlim',
            help='discard fastlim results',
            action="store_true" )
    argparser.add_argument ( '-N', '--notimestamp',
            help='dont put a timestamp on it',
            action="store_true" )
    args=argparser.parse_args()
    args.drawtimestamp = not args.notimestamp
    args.miscol = "gold" ## missing likelihood color, golden
    args.miscol = "white" ## missing likelihood color, white
    args.diagcol = "black"
    args.diagcol = "grey"
    outputfile = draw( vars ( args ) )
    show ( outputfile )