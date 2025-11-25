#!/usr/bin/env python3

""" the closure plots, see e.g. closure_Xt.png """

from typing import Union
import numpy as np
from math import exp
import os, glob

def getAllModels( directory : os.PathLike = "../data/fake_stops1",
       minF : float = .8, maxfiles : Union[None,int] = None,
       burnin : int = 0 ) -> tuple[list[dict],float]:
    """ get all the model dictionaries
    :param directory: directory to search for
    :param minF: ignore all points with K < minF * maxK
    :param maxfiles: if not none, the cap on the number of files
    :param burnin: if not none, then throw away this number of burnin steps
    """
    dictfpattern = f"{directory}/dictfiles/pmodel_*.dict"
    dictfiles = glob.glob ( dictfpattern )
    print ( f"[plotClosure] found {len(dictfiles)} files in {dictfpattern}" )
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
    maxK = max ( x["K"] for x in all_models )
    filtered = [d for d in all_models if d["K"] > minF * maxK ]
    return filtered, maxK

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

def getCoordinates ( models : Union[list,dict], maxK : float, minF : float,
       coords : dict ):
    xvalues, yvalues, Kvalues = [], [], []
    xpid, ypid = coords["x"], coords["y"]
    if type(models) == dict:
        if xpid in models["masses"] and ypid in models["masses"]:
            return models["masses"][xpid], models["masses"][ypid], models["K"]
        return float("nan"), float("nan"), float("nan")
    print ( f"[plotClosure] we have {len(models)} models" )
    for model in models:
        if xpid in model["masses"] and ypid in model["masses"]:
            xvalues.append ( model["masses"][ xpid ] )
            yvalues.append ( model["masses"][ ypid ] )
            K = model["K"]
            if K == None:
                K = float("nan")
            else:
                pass
                # K = K**2 / (maxK**2)
                K = max ( 1, 10 * ( K - minF * maxK ) )
                # K /= 100.
                # K = exp(K) / exp(60.264) * 80
            Kvalues.append ( K )
    return xvalues, yvalues, Kvalues

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

def plotClosure( args : dict ):
    """ closure plot

    args (dict):
        path (os.PathLike) - the path to all the data
        minF (float) - the minimum fraction of K to plot a point
    """
    from matplotlib import pyplot as plt
    # models = getAllPModels( path )
    print ( f"[plotClosure] obtaining data from {args['path']}" )
    models, maxK = getAllModels( args["path"], args["minF"], args["maxfiles"],
            args["burnin"] )
    truth = getTruthModel( args["path"] )
    coords={ "x": 1000023, "y": 1000022, "type_x": "mass", "type_y": "mass" }
    sinjection = "ewkino"
    if truth is not None and 1000006 in truth["masses"]:
        coords[ "x" ] = 1000006
        sinjection = "stop"
    if truth is not None:
        x_true, y_true, K_true = getCoordinates ( truth, maxK, 0., coords )
    splitm = splitByDictFile ( models )
    colors = plt.cm.viridis(np.linspace(0.3, 0.9, len(splitm)))
    for i,(df,models) in enumerate(splitm.items()):
        x, y, K = getCoordinates ( models, maxK, args["minF"], coords )
        plt.scatter ( x, y, s = K, alpha=0.5, color = colors[i],
                      label=f"walker #{df}", edgecolors = darken(colors[i]) )
    if truth is not None:
        plt.scatter ( x_true, y_true, s=280, marker="+", color="white",
            linewidths=4 )
        plt.scatter ( x_true, y_true, s=140, marker="+", color="red",
            label="truth" )
    loc = "best"
    loc = "lower right"
    plt.legend( loc = loc )
    from ptools.sparticleNames import SParticleNames
    namer = SParticleNames()
    filename = f"closure_{namer.asciiName(coords['x'])}.png"
    plt.xlabel ( f"mass, ${namer.texName(coords['x'])}$ [GeV]" )
    plt.ylabel ( f"mass, ${namer.texName(coords['y'])}$ [GeV]" )
    plt.title ( f"closure test, {sinjection} injection" )
    from smodels_utils.helper.various import pngMetaInfo
    metadata = pngMetaInfo()
    from installation import version as protomodels_version
    metadata["protomodels"] = protomodels_version()
    plt.savefig ( filename, metadata= metadata )
    from smodels_utils.plotting.mpkitty import timg
    timg ( filename )
    if args["interact"]:
        import sys, IPython; IPython.embed( colors = "neutral" ); sys.exit()

if __name__ == "__main__":
    import argparse
    argparser = argparse.ArgumentParser(description="plots closure tests")
    argparser.add_argument ( '-p', '--path',
            help='path [../data/fake_stops1/]', type=str,
            default='../data/fake_stops1/' )
    argparser.add_argument ( '-i', '--interact',
            help='interactive shell', action="store_true" )
    argparser.add_argument ( '-m', '--minF',
            help='minF [0.7]',
            type=float, default=0.7 )
    argparser.add_argument ( '--maxfiles', type=int,
            help='maximum numbers of files [None]', default=None )
    argparser.add_argument ( '--burnin', type=int,
            help='throw away first n burnin steps [None]', default=0 )
    args=vars ( argparser.parse_args() )
    #path = "../data/fake_ewk1/"
    #path = "../data/fake_ewkoff1/"
    #minF = .1
    # path = "../data/fake_stops1/"
    # args = { "path": path, "minF": minF, "interact": interact }
    plotClosure( args )
