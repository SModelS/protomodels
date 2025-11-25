#!/usr/bin/env python3

""" the posteriori plots, see e.g. posterior_Xt.png """

from typing import Union
import numpy as np
from math import exp
import os, glob

def getAllModels( directory : os.PathLike = "../data/fake_stops1",
       minF : float = .8, maxfiles : Union[None,int] = None,
       burnin : int = 0 ) -> tuple[list[dict],float]:
    """ get all the model dictionaries
    :param directory: directory to search for
    :param minF: ignore all points with K < minF * minK
    :param maxfiles: if not none, the cap on the number of files
    :param burnin: if not none, then throw away this number of burnin steps
    """
    dictfpattern = f"{directory}/dictfiles/pmodel_*.dict"
    dictfiles = glob.glob ( dictfpattern )
    print ( f"[plotPosterior] found {len(dictfiles)} files in {dictfpattern}" )
    from natsort import natsorted
    dictfiles = natsorted ( dictfiles )
    all_models = []
    if maxfiles != None:
        dictfiles = dictfiles[:maxfiles]
    for dictfile in dictfiles[:]:
        dfname = os.path.basename ( dictfile )
        dfname = dfname.replace("pmodel_","").replace(".dict","")
        with open ( dictfile, "rt" ) as f:
            models = eval ( f.read() )
            if burnin != None:
                models = models[burnin:]
            for model in models:
                model["dictfile"]=dfname
                if model["Accepted"]==0:
                    all_models.append ( model )
    minK = min ( x["K"] for x in all_models )
    print ( f"[plotPosterior] minK is {minK:.2f}" )
    filtered = [d for d in all_models if d["K"] is not None and d["K"] > (1 + minF ) * minK ]
    print ( f"[plotPosterior] filtered from {len(all_models)} to {len(filtered)}" )
    return filtered, minK

def getAllPModels( directory : os.PathLike = "../data/fake_stops1" ) -> \
                                              list[dict]:
    """ get all the model dictionaries """
    dictfpattern = f"{directory}/Pmodels/pmodel*.dict"
    dictfiles = glob.glob ( dictfpattern )
    all_models = []
    for dictfile in dictfiles:
        with open ( dictfile, "rt" ) as f:
            model = eval ( f.read() )
            #    if model["Accepted"]==0:
            all_models.append ( model )
    return all_models

def splitByDictFile ( models : list ):
    ret = {}
    for model in models:
        df = model["dictfile"]
        if not df in ret:
            ret[df]=[]
        ret[df].append ( model )
    return ret

def getCoordsFromModel ( model : dict, coords : dict ) -> dict:
    """ get the right coordinates from a single model """
    xpid, ypid = coords["x"], coords["y"]
    if xpid in model["masses"] and ypid in model["masses"]:
        return { "x": model["masses"][xpid], "y": model["masses"][ypid], 
                 "K": model["K"], "legit": True }
    return { "x": float("nan"), "y": float("nan"), "K": float("nan"),
             "legit": False }

def getCoordinates ( models : Union[list,dict], minK : float, minF : float,
       coords : dict ) -> dict:
    if type(models) == dict:
        mcoords = getCoordsFromModel ( models, coords )
        return { "x": [ mcoords["x"] ], "y": [ mcoords["y"] ], 
                 "K": [ mcoords["K"] ], "w": 1 }
    print ( f"[plotPosterior] we have {len(models)} models" )
    points = {}
    for model in models:
        mcoords = getCoordsFromModel ( model, coords )
        if mcoords["legit"]==True:
            xv, yv =  mcoords[ "x" ], mcoords[ "y" ]
            K = mcoords["K"] - minK
            hashCode = 1e6*xv+yv
            if not hashCode in points:
                points[hashCode] = { "x": xv, "y": yv, "w": 0, "K": K }
            points[hashCode]["w"]+=1
    xvalues, yvalues, Kvalues, weights = [], [], [], []
    for hashC,point in points.items():
        xv, yv =  point[ "x" ], point[ "y" ]
        xvalues.append ( xv )
        yvalues.append ( yv )
        weights.append ( point["w"] )
        Kvalues.append ( point["K"] )
    return { "x": xvalues, "y": yvalues, "w": weights, "K": Kvalues }

def getTruthModel( directory : os.PathLike ):
    truthfile = f"{directory}/truth.dict"
    if not os.path.exists ( truthfile ):
        return None
    with open ( truthfile, "rt" ) as f:
        model = eval ( f.read() )
        return model

# Function to darken colors
def darken(color, amount=0.6):
    """Return a darker shade of a given matplotlib color."""
    import matplotlib.colors as mcolors
    c = mcolors.to_rgb(color)
    return tuple(max(0, min(1, i * amount)) for i in c)

def getPlottingCoords ( args, truth ):
    """ determine what gets plotted on what axis """
    coords={ "x": 1000023, "y": 1000022, "type_x": "mass", "type_y": "mass" }
    if truth is not None and 1000006 in truth["masses"]:
        coords[ "x" ] = 1000006
    if "xcoordinate" in args and args["xcoordinate"] is not None:
        xc = args["xcoordinate"]
        if xc.startswith ( "M" ):
            coords["x"]=int(xc[1:])
            coords["type_x"]="mass"
        if xc.startswith ( "S" ):
            xc = xc.lower().replace("ssms","").replace("ssm","").replace("s","")
            coords["x"]=eval(xc)
            coords["type_x"]="ssm"
    if "ycoordinate" in args and args["ycoordinate"] is not None:
        yc = args["ycoordinate"]
        if yc.startswith ( "M" ):
            coords["y"]=int(yc[1:])
            coords["type_y"]="mass"
        if yc.startswith ( "S" ):
            yc = yc.lower().replace("ssms","").replace("ssm","").replace("s","")
            coords["y"]=eval(yc)
            coords["type_y"]="ssm"
    return coords

def plotPosterior( args : dict ):
    """ closure plot

    args (dict):
        path (os.PathLike) - the path to all the data
        minF (float) - the minimum fraction of K to plot a point
    """
    dirname = os.path.abspath ( os.path.dirname(__file__) + "/../" )
    args["path"]=args["path"].replace("__file__", dirname )
    from matplotlib import pyplot as plt
    # models = getAllPModels( path )
    print ( f"[plotPosterior] obtaining data from {args['path']}" )
    models, minK = getAllModels( args["path"], args["minF"], args["maxfiles"],
            args["burnin"] )
    truth = getTruthModel( args["path"] )
    coords = getPlottingCoords ( args, truth )
    sinjection = "ewkino"
    if coords["x"] == 1000006:
        sinjection = "stop"
    if truth is not None:
        true_coords = getCoordinates ( truth, minK, 0., coords )
    # splitm = splitByDictFile ( models )
    # colors = plt.cm.viridis(np.linspace(0.3, 0.9, len(splitm)))
    dcoords = getCoordinates ( models, minK, args["minF"], coords )
    x,y,w = dcoords["x"], dcoords["y"], dcoords["w"]
    plt.scatter ( x, y, s = np.sqrt(w),
            alpha=0.5, color = "green",
            edgecolors = "darkgreen" )
    from scipy.stats import gaussian_kde
    # x, y are your 1D arrays of data points
    # x, y = ...

    # stack data for KDE
    data = np.vstack([x, y])
    kde = gaussian_kde(data,weights=w)

    # create a grid for evaluating the KDE
    xmin, xmax = min(x), max(x)
    ymin, ymax = min(y), max(y)
    xx, yy = np.meshgrid(
        np.linspace(xmin, xmax, 300),
        np.linspace(ymin, ymax, 300)
    )

    # evaluate KDE on grid
    zz = kde(np.vstack([xx.ravel(), yy.ravel()])).reshape(xx.shape)

    # compute contour levels for 67%, 95%, 100%
    # sort the density values from high to low
    z_sorted = np.sort(zz.ravel())[::-1]
    cumulative = np.cumsum(z_sorted)
    cumulative /= cumulative[-1]

    def find_level(threshold):
        """Returns the density value such that the enclosed probability is 'threshold'."""
        return z_sorted[np.searchsorted(cumulative, threshold)]

    # clevels= [ 0.99, 0.95, 0.67 ]
    # clevels= [ 0.99, 0.86, 0.39 ]
    # .5 sigma 0.1175030974
    # 1 sigma 0.3934693403
    # 2 sigma 0.8646647168
    # 3 sigma 0.9888910035
    # 4 sigma 0.9996645374
    clevels= [ 0.9888910035, 0.8646647168, 0.3934693403 ]
    clevels= [ 0.95, 0.8646647168, 0.3934693403 ]

    levels = [find_level(p) for p in clevels ]

    # plot
    # plt.figure(figsize=(7,6))
    # plt.scatter(x, y, s=5, alpha=0.3)  # original points

    cs = plt.contour(
        xx, yy, zz,
        levels=levels,
        colors=['0.4', '0.2', '0.0'],
        linewidths=2
    )

    labels = plt.clabel(cs, inline=True, fontsize=12,
        fmt={lvl: f"{p*100:.0f}%" for lvl, p in zip(levels, clevels )})

    for txt in labels:
        txt.set_fontweight('bold')

    ## change the xrange
    # plt.xlim ( 300, 2300 )

    # plt.contour ( X, Y, Z )
    if truth is not None:
        plt.scatter ( x_true, y_true, s=280, marker="+", color="white",
            linewidths=4 )
        plt.scatter ( x_true, y_true, s=140, marker="+", color="red",
            label="truth" )
    loc = "best"
    loc = "lower right"
    # plt.legend( loc = loc )
    from ptools.sparticleNames import SParticleNames
    namer = SParticleNames()
    filename = args["outfile"]
    filename = filename.replace( "@@X@@", namer.asciiName(coords['x']) )
    filename = filename.replace( "@@Y@@", namer.asciiName(coords['y']) )
    plt.xlabel ( f"mass, ${namer.texName(coords['x'])}$ [GeV]" )
    plt.ylabel ( f"mass, ${namer.texName(coords['y'])}$ [GeV]" )
    plt.title ( f"a posteriori distribution" )
    print ( f"[plotPosterior] saving to {filename}" )
    plt.savefig ( filename )
    from smodels_utils.plotting.mpkitty import timg
    timg ( filename )
    if args["interact"]:
        import sys, IPython; IPython.embed( colors = "neutral" ); sys.exit()

if __name__ == "__main__":
    import argparse
    argparser = argparse.ArgumentParser(description="plots closure tests")
    argparser.add_argument ( '-p', '--path',
            help='path [../data/rundir4/]', type=str,
            default='__file__/data/rundir4/' )
    argparser.add_argument ( '-i', '--interact',
            help='interactive shell', action="store_true" )
    argparser.add_argument ( '-m', '--minF',
            help='minF [0.1]',
            type=float, default=0.1 )
    argparser.add_argument ( '-x', '--xcoordinate',
            help='what to plot on the x axis, e.g. "M1000006". None is automatic. [None]',
            type=str, default=None )
    argparser.add_argument ( '-y', '--ycoordinate',
            help='what to plot on the y axis, e.g. "SSM10000061000006". None is automatic [None]',
            type=str, default=None )
    argparser.add_argument ( '-o', '--outfile',
            help='Name of output file, replacing @@X@@ and @@Y@@ [posterior_@@X@@_@@Y@@.png]',
            type=str, default="posterior_@@X@@_@@Y@@.png" )
    argparser.add_argument ( '--maxfiles', type=int,
            help='maximum numbers of files [None]', default=None )
    argparser.add_argument ( '--burnin', type=int,
            help='throw away first n burnin steps [100]', default=100 )
    args=vars ( argparser.parse_args() )
    #path = "../data/fake_ewk1/"
    #path = "../data/fake_ewkoff1/"
    #minF = .1
    # path = "../data/fake_stops1/"
    # args = { "path": path, "minF": minF, "interact": interact }
    plotPosterior( args )
