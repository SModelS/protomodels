#!/usr/bin/env python3

__all__ = [ "draw", "show" ]

import sys, os, time, math
sys.path.insert(0,"../")
from smodels.experiment.databaseObj import Database
from smodels.base.smodelsLogging import setLogLevel
from smodels.base.physicsUnits import TeV
# from smodels.tools.colors import colors
from smodels.base.smodelsLogging import colors
from smodels_utils.helper.various import hasLLHD
from tester import analysisCombiner
import IPython
from typing import Union, Dict, List
from os import PathLike
import subprocess
from tester.combinationsmatrix import getYamlMatrix

def getCombinationsMatrix ( path : Union[None,Dict,PathLike] ):
    """ get the combinations matrix. If path is matrix dictionary itself, return it.
        If path is None, retrieve matrix from tester.combinationsmatrix.getYamlMatrix.
    """

    if type ( path ) == type ( None ):
        return getYamlMatrix()
    if type ( path ) == dict:
        return path
    import importlib
    spec = importlib.util.spec_from_file_location( "getYamlMatrix", path )
    imp = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(imp)
    return imp.getYamlMatrix()

def sortBySqrts ( results : List, sqrts : TeV ) -> List:
    ret = []
    for res in results:
        if abs (res.globalInfo.sqrts.asNumber(TeV) - sqrts ) < 0.1:
            ret.append ( res )
    return ret

def noFastlim ( results : List ) -> List:
    """ remove fastlim results, FIXME should recycle method from smodels-utils """
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

def checkForPartialCombinability ( e1, e2 ) -> bool:
    """ check if a and b are partially combinable """
    ads = e1.datasets
    bds = e2.datasets
    for a in ads:
        for b in bds:
            if a.isCombinableWith ( b ):
                # print ( "a and b!", e1.globalInfo.id, a, e2.globalInfo.id, b )
                return True
    return False

def filterResults ( results : list, excludes : list, renames : dict ) -> list:
    """ given the list of excludes, filter results """
    import fnmatch
    ret = []
    for result in results:
        anaId = result.globalInfo.id
        if anaId in renames.keys():
            result.globalInfo.id = renames[anaId]
        isExcluded = False
        for exclude in excludes:
            if fnmatch.fnmatch ( anaId, exclude ):
                if exclude == anaId:
                    print ( f"[plotBAM] dropping {anaId}" )
                else:
                    print ( f"[plotBAM] dropping {anaId}: matches {exclude}" )
                isExcluded = True
                break
        if not isExcluded:
            ret.append ( result )
    return ret

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
    combinationsmatrix, status = getCombinationsMatrix ( args["combinationsmatrix"] )
    if not combinationsmatrix or status != 0:
        sys.exit("Combination matrix not loaded correctly.")

    sqrtses = [ 8, 13 ]
    if args["sqrts"] not in [ "all" ]:
        sqrtses = [ int(args["sqrts"]) ]

    colors.on = True
    setLogLevel ( "debug" )

    # dir = "/home/walten/git/smodels-database/"
    dbdir = args["database"]
    d=Database( dbdir, combinationsmatrix = combinationsmatrix )
    print(d)
    analysisIds = [ "all" ]
    if "analyses" in args and args["analyses"]!=None:
        analysisIds = args["analyses"].split(",")
    exps = [ "CMS", "ATLAS" ]
    if args["experiment"] in [ "CMS", "ATLAS" ]:
        analysisIds = [ args["experiment"]+"*" ]
        exps = [ args["experiment"] ]
    dataTypes = [ "all" ]
    if args["effmaps_only"]:
        dataTypes = [ "efficiencyMap" ]
    results = d.getExpResults( analysisIDs = analysisIds, dataTypes = dataTypes )
    if args["nofastlim"]:
        results = noFastlim ( results )
    results = sortOutDupes ( results )
    if args["sqrts"] in [ "8", "13" ]:
        results = sortBySqrts ( results, int(args["sqrts"]) )

    excludes = args["exclude"].split(",")
    renames = args["rename"]
    if renames not in [ None, "" ]:
        renames = eval(renames)
    else:
        renames = {}
    results = filterResults ( results, excludes, renames )

    #results.sort()
    nres = len ( results )

    from matplotlib import pyplot as plt
    import matplotlib
    matplotlib.use('agg')
    labelsize = 14
    # x- and y- tickpads are to adjust the position of the analysis id labels
    xtickpad, ytickpad = -55, -55
    xoff_title, yoff_title = .35, .95
    ts_off_x = .6 ## timestamp offset x
    ts_off_y = 0. ## timestamp offset y
    if nres < 60:
        xtickpad = 0
        ytickpad = 0
        labelsize = 26
    if nres < 30:
        xtickpad = 0
        ytickpad = 0
        labelsize = 40
        xoff_title, yoff_title = .47, .9
        ts_off_x = -.01
    if nres < 18:
        ts_off_x = -.05
        ts_off_y = -.0
        xtickpad = 0
        ytickpad = 0
        labelsize = 60
    if nres < 10:
        ts_off_x = -.2
        ts_off_y = -.1
        xoff_title, yoff_title = .35, .95
        xtickpad = 0 # -580
        ytickpad = 0#  -580
        labelsize = 100
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
        #    print ( f"no likelihood: {label}" )
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
    c = [ "b", "tab:green", "tab:red", "tab:orange", "white", "tab:grey" ]
    # c[3]="darkgreen"
    v = np.arange(0.,1.00001,1. / (len(c)-1) )
    l = list(zip(v,c))
    # print ( "l", l )
    from  matplotlib.colors import LinearSegmentedColormap
    cmap=LinearSegmentedColormap.from_list('rg',l, N=len(c) )
    plt.matshow ( h, aspect = "equal", origin = "lower", cmap = cmap,
                  vmin = 0, vmax = 5. )
    drawGrid = True
    if drawGrid:
        # This is very hack-ish
        #xticks_ = plt.gca().get_xticks()[1:-1]
        # xticks = [x - 0.5 for x in plt.gca().get_xticks()][1:-1]
        #print ( "xticks", xticks_ )
        #xticks = list ( range ( int(min(xticks_)), int(max(xticks_)+2 )) )
        xticks = [ x+.5 for x in range(len(results)) ]
        plt.gca().set_xticks( xticks, minor='true')
        plt.gca().set_yticks( xticks, minor='true')
        # plt.grid ( visible = True )
        plt.grid(which='minor')
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
        title = args["title"]
        title = title.replace("@exp@",exps[0]).replace("@sqrts@",str(sqrtses[0]))
        plt.text ( xoff_title, yoff_title, title,
                   fontsize = 2 * labelsize, transform = fig.transFigure,
                   horizontalalignment = "center" )
    ct = 0
    for ana in exps:
        for sqrts in sqrtses:
            name= f"{ana}{sqrts}"
            xcoord = .5 * ( bins[ana][sqrts][0] + bins[ana][sqrts][1] )
            ycoord = n- .5 * ( bins[ana][sqrts][0] + bins[ana][sqrts][1] ) -3
            if len(sqrtses)>1 or len(exps)>1:
                plt.text(-5,xcoord-3,f"{ana}\n{sqrts} TeV", fontsize=44,
                         c="black", rotation=90, horizontalalignment="center" )
                plt.text(ycoord,-8, f"{ana}\n{sqrts} TeV",
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
        t = time.strftime("%h %d %Y" )
        dbver = d.databaseVersion
        plt.text ( ts_off_x, ts_off_y, f"plot produced {t}\nfrom database v{dbver}",
                   va="bottom", c="lightgrey", transform = fig.transFigure, 
                   fontsize = .5 * labelsize )
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
    # plt.tight_layout( )
    plt.savefig ( outputfile, dpi=dpi, bbox_inches="tight" )
    if "trim" in args and args["trim"]:
        cmd = f"convert {outputfile} -trim trimmed.{outputfile}"
        subprocess.getoutput ( cmd )
        cmd = f"mv trimmed.{outputfile} {outputfile}"
        subprocess.getoutput ( cmd )
    return outputfile

def show ( outputfile ):
    import shutil
    if shutil.which ( "timg" ) != None:
        cmd = f"timg {outputfile}"
        # print ( cmd )
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
            help='path to combinationsmatrix file (will call getYamlMatrix() within that file). If none, get it from protomodels.tester.combinationsmatrix.getYamlMatrix() [None]',
            type=str, default=None )
    argparser.add_argument ( '-e', '--experiment', nargs='?',
            help='plot only specific experiment CMS,ATLAS,all [all]',
            type=str, default='all' )
    argparser.add_argument ( '-s', '--sqrts', nargs='?',
            help='plot only specific sqrts 8,13,all [all]',
            type=str, default='all' )
    argparser.add_argument ( '--exclude',
            help='exclude this comma-separated list of analysis, wildcards allowed [none]',
            type=str, default='' )
    argparser.add_argument ( '--rename',
            help="dictionary (given as string) of analyses to rename, e.g.: { 'ATLAS-SUSY-2018-22-multibin': 'ATLAS-SUSY-2018-22' } [none]",
            type=str, default=None )
    argparser.add_argument ( '-o', '--outputfile', nargs='?',
            help='outputfile (@M gets replaced by [experiment][sqrts]) [matrix@M.png]',
            type=str, default='matrix@M.png' )
    argparser.add_argument ( '-a', '--analyses',
            help='select for comma separated list of analyses [None]',
            type=str, default=None )
    argparser.add_argument ( '--title',
            help='specify the title [@exp@, @sqrts@ TeV]',
            type=str, default="@exp@, @sqrts@ TeV" )
    argparser.add_argument ( '-t', '--triangular',
            help='plot as lower triangle matrix?',
            action="store_true" )
    argparser.add_argument ( '--effmaps_only',
            help='plot only for efficiency map results',
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
    argparser.add_argument ( '--show', help='show plot', action="store_true" )
    args=argparser.parse_args()
    args.drawtimestamp = not args.notimestamp
    outputfile = draw( vars ( args ) )
    if args.show:
        show ( outputfile )
