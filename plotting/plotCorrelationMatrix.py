#!/usr/bin/env python3

__all__ = [ "draw", "show" ]

import sys, os, time, math
sys.path.insert(0,"../")
from smodels.experiment.databaseObj import Database
from smodels.tools.smodelsLogging import setLogLevel
from smodels.tools.physicsUnits import TeV
from smodels.tools.colors import colors
from smodels_utils.helper.various import hasLLHD
from tester import analysisCombiner
import IPython
from typing import Union, Dict
from os import PathLike
import subprocess

def getCombinationsMatrix ( path : Union[None,Dict,PathLike] ):
    """ get the combinations matrix. If path is matrix dictionary itself, return it.
        If path is None, retrieve matrix from tester.combinationsmatrix.getMatrix.
    """
       
    if type ( path ) == type ( None ):
        from tester.combinationsmatrix import getMatrix
        return getMatrix()
    if type ( path ) == dict:
        return path
    import importlib
    spec = importlib.util.spec_from_file_location( "getMatrix", path )
    imp = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(imp)
    return imp.getMatrix()

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

def checkForPartialCombinability ( e1, e2 ):
    """ check if a and b are partially combinable """
    ads = e1.datasets
    bds = e2.datasets
    for a in ads:
        for b in bds:
            if a.isCombinableWith ( b ):
                # print ( "a and b!", e1.globalInfo.id, a, e2.globalInfo.id, b )
                return True
    return False

def draw( args : dict ):
    """
    draw the correlation matrix
    :param args: dictionary of args
           triangular: if True, then only plot the upper triangle of this
                         symmetrical matrix
           experiment: draw only for specific experiment ("CMS", "ATLAS", "all" )
           database: path to database
           sqrts: draw only for specific sqrts ( "8", "13", "all" )
           drawtimestamp: if true, put a timestamp on plot
           outputfile: file name of output file (matrix.png)
           nofastlim: if True, discard fastlim results
    """
    matrix = getCombinationsMatrix ( args["combinationsmatrix"] )
    sqrtses = [ 8, 13 ]
    if args["sqrts"] not in [ "all" ]:
        sqrtses = [ int(args["sqrts"]) ]

    colors.on = True
    setLogLevel ( "debug" )

    # dir = "/home/walten/git/smodels-database/"
    dbdir = args["database"]
    d=Database( dbdir, discard_zeroes=False, combinationsmatrix = matrix )
    print(d)
    analysisIds = [ "all" ]
    if "analyses" in args and args["analyses"]!=None:
        analysisIds = args["analyses"].split(",")
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
    matplotlib.use('agg')
    labelsize = 14
    # x- and y- tickpads are to adjust the position of the analysis id labels
    xtickpad, ytickpad = -55, -55 
    if nres < 60:
        xtickpad = 0
        ytickpad = 0
        labelsize = 26
    if nres < 5:
        xtickpad = -580
        ytickpad = -580
        labelsize = 80
    matplotlib.rc('xtick', labelsize=labelsize, labelcolor = "gray" )
    matplotlib.rc('ytick', labelsize=labelsize, labelcolor = "gray" )

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
        ymax=0
        if x < bins[ana][sqrts][0]:
            bins[ana][sqrts][0]=x
        if x > bins[ana][sqrts][1]:
            bins[ana][sqrts][1]=x
            ymax=x
        label = label.replace("-agg","")
        if len(exps)==1 and len(sqrtses)==1:
            label = label.replace("CMS-","").replace("ATLAS-","").replace("-agg","")
        labels.append ( label )
        for y,f in enumerate(results):
            if args["triangular"] and y<x:
                h[x][n-y-1]= float("nan")
                continue
            isComb = e.isCombinableWith ( f )
            partial = False
            if not isComb:
                partial = checkForPartialCombinability ( e, f )
            #if partial:
            #    print ( f"{label}+{f.globalInfo.id}: partially combinable" )
            # sys.exit()
            #isUn = analysisCombiner.canCombine ( e.globalInfo, f.globalInfo, 
            #        args["strategy"] )
            # isUn = e.isUncorrelatedWith ( f )
            v = 0.
            if isComb:
                v = 1.
            else:
                v = 2.
                if partial:
                    v = 3.
            if not hasLikelihood or not hasLLHD ( f ): ## has no llhd? cannot be combined
                v = 4.
            if y==x:
                v = 5.
            h[x][n-y-1]= v
            # h[n-x-1][y]= v

    c = [ "b", "limegreen", "red", "orange", "white", "grey" ]
    # c[3]="darkgreen"
    v = np.arange(0.,1.00001,1. / (len(c)-1) )
    l = list(zip(v,c))
    # print ( "l", l )
    from  matplotlib.colors import LinearSegmentedColormap
    cmap=LinearSegmentedColormap.from_list('rg',l, N=len(c) )
    plt.matshow ( h, aspect = "equal", origin = "lower", cmap = cmap,
                  vmin = 0, vmax = 5. )
    plt.grid ( visible = False )
    # plt.xticks ( rotation=90, horizontalalignment="center" )
    fig = plt.gcf()
    fig.set_size_inches(30, 30)
    ax = plt.gca()
    ax.xaxis.set_ticks_position("bottom")
    ax.tick_params(axis='x', pad=xtickpad )
    ax.tick_params(axis='y', pad=ytickpad )
    plt.setp(ax.get_xticklabels(), rotation=90,
         ha="center", rotation_mode="default")
    ax.set_xticks ( range(len(labels)) )
    labels.reverse()
    ax.set_xticklabels( labels )
    ax.set_yticks ( range(len(labels)) )
    labels.reverse()
    ax.set_yticklabels( labels ) ## need to invert
    if len(exps)==1 and len(sqrtses)==1:
        plt.text ( .45, .95, "%s, %d TeV" % ( exps[0], sqrtses[0] ),
                   fontsize = 3 * labelsize, transform = fig.transFigure )
    ct = 0
    for ana in exps:
        for sqrts in sqrtses:
            name= "%s%d" % ( ana, sqrts )
            xcoord = .5 * ( bins[ana][sqrts][0] + bins[ana][sqrts][1] )
            ycoord = n- .5 * ( bins[ana][sqrts][0] + bins[ana][sqrts][1] ) -3
            if len(sqrtses)>1 or len(exps)>1:
                plt.text(-5,xcoord-3,"%s\n%d TeV" % ( ana, sqrts ),
                         fontsize=44, c="black", rotation=90,
                         horizontalalignment="center" )
                plt.text(ycoord,-8,"%s\n%d TeV" % ( ana, sqrts ) ,
                         fontsize=44, c="black", horizontalalignment="center" )
            yt = bins[ana][sqrts][1] +1
            extrudes = 3 # how far does the line extrude into tick labels?
            xmax = n
            if args["triangular"]:
                xmax = n-yt
            lc = "black"
            alpha = 1
            if nres<5:
                lc = "white"
                alpha = 0.
            ymax = n
            if args["triangular"]:
                ymax = yt
            for s in [ "bottom", "top", "left", "right" ]:
                ax.spines[s].set_visible(False)
            if ct>0:
                plt.plot ( [ -extrudes, xmax ], [ yt-.5, yt-.5 ], c=lc, alpha=alpha )
                plt.plot ( [ n-yt-.5, n-yt-.5], [ymax, -extrudes ], c=lc, alpha=alpha )
            ct += 1
    if args["drawtimestamp"]:
        plt.text ( .01, .01, "plot produced %s from database v%s" % \
                   ( time.strftime("%h %d %Y" ), d.databaseVersion ), 
                   c="grey", transform = fig.transFigure, fontsize=24 )
    outputfile = args["outputfile"]
    if "@M" in outputfile:
        modifiers = ""
        if len(exps)==1:
            modifiers += exps[0]
        if len(sqrtses)==1:
            modifiers += str(sqrtses[0])
        outputfile = outputfile.replace("@M",modifiers)
    print ( f"Plotting to {outputfile}" )
    dpi = 200
    if nres < 5:
        ## fewer than 5 results? make it very small!
        dpi = 15
    plt.tight_layout( )
    plt.savefig ( outputfile, dpi=dpi )
    if "trim" in args and args["trim"]:
        cmd = f"convert {outputfile} -trim trimmed.{outputfile}"
        subprocess.getoutput ( cmd )
        cmd = f"mv trimmed.{outputfile} {outputfile}"
        subprocess.getoutput ( cmd )
    return outputfile

def show ( outputfile ):
    cmd = f"timg {outputfile}"
    o = subprocess.getoutput ( cmd )
    print ( o )

def plotHandCrafted():
    """ modify this to produce your special version of this plot """
    sys.exit()

if __name__ == "__main__":
    # plotHandCrafted()
    import argparse
    argparser = argparse.ArgumentParser(description="correlation/combination matrix plotter")
    argparser.add_argument ( '-S', '--strategy', nargs='?',
            help='combination strategy [aggressive]', type=str, default='aggressive' )
    # dbpath = "../../smodels-database"
    dbpath = "official"
    argparser.add_argument ( '-d', '--database', nargs='?',
            help=f'path to database [{dbpath}]', type=str, default=dbpath )
    argparser.add_argument ( '-c', '--combinationsmatrix', nargs='?',
            help='path to combinationsmatrix file (will call getMatrix() within that file). If none, get it from protomodels.tester.combinationsmatrix.getMatrix() [None]',
            type=str, default=None )
    argparser.add_argument ( '-e', '--experiment', nargs='?',
            help='plot only specific experiment CMS,ATLAS,all [all]',
            type=str, default='all' )
    argparser.add_argument ( '-s', '--sqrts', nargs='?',
            help='plot only specific sqrts 8,13,all [all]',
            type=str, default='all' )
    argparser.add_argument ( '-o', '--outputfile', nargs='?',
            help='outputfile (@M gets replaced by [experiment][sqrts]) [matrix@M.png]',
            type=str, default='matrix@M.png' )
    argparser.add_argument ( '-a', '--analyses', 
            help='select for comma separated list of analyses [None]',
            type=str, default=None )
    argparser.add_argument ( '-t', '--triangular',
            help='plot as lower triangle matrix?',
            action="store_true" )
    argparser.add_argument ( '-T', '--trim',
            help='trim the figure in the end',
            action="store_true" )
    argparser.add_argument ( '-n', '--nofastlim',
            help='discard fastlim results',
            action="store_true" )
    argparser.add_argument ( '-N', '--notimestamp',
            help='dont put a timestamp on it',
            action="store_true" )
    args=argparser.parse_args()
    args.drawtimestamp = not args.notimestamp
    outputfile = draw( vars ( args ) )
    show ( outputfile )
