#!/usr/bin/env python3

from typing import Union
import numpy as np
from math import exp
import os, glob

def getAllModels( directory : os.PathLike = "../data/fake_stops1",
       minF : float = .8 ) -> tuple[list[dict],float]:
    """ get all the model dictionaries 
    :param directory: directory to search for
    :param minF: ignore all points with K < minF * maxK
    """
    dictfpattern = f"{directory}/dictfiles/pmodel_*.dict"
    dictfiles = glob.glob ( dictfpattern )
    all_models = []
    for dictfile in dictfiles[:10]:
        dfname = os.path.basename ( dictfile ) 
        dfname = dfname.replace("pmodel_","").replace(".dict","")
        with open ( dictfile, "rt" ) as f:
            models = eval ( f.read() )
            for model in models[500:]:
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
    print ( f"[plotClosure] obtaining data from {path}" )
    models, maxK = getAllModels( args["path"], args["minF"] )
    truth = getTruthModel( args["path"] )
    coords={ "x": 1000023, "y": 1000022, "type_x": "mass", "type_y": "mass" }
    sinjection = "ewkino"
    if 1000006 in truth["masses"]:
        coords[ "x" ] = 1000006
        sinjection = "stop"
    x_true, y_true, K_true = getCoordinates ( truth, maxK, 0., coords )
    splitm = splitByDictFile ( models )
    colors = plt.cm.viridis(np.linspace(0.3, 0.9, len(splitm)))
    for i,(df,models) in enumerate(splitm.items()):
        x, y, K = getCoordinates ( models, maxK, minF, coords )
        plt.scatter ( x, y, s = K, alpha=0.5, color = colors[i],
                      label=df, edgecolors = darken(colors[i]) )
    plt.scatter ( x_true, y_true, s=280, marker="+", color="white",linewidths=4 )
    plt.scatter ( x_true, y_true, s=140, marker="+", color="red", 
                  label="truth" )
    plt.legend()
    from ptools.sparticleNames import SParticleNames
    namer = SParticleNames()
    filename = f"closure_{namer.asciiName(coords['x'])}.png"
    plt.xlabel ( f"mass, ${namer.texName(coords['x'])}$ [GeV]" )
    plt.ylabel ( f"mass, ${namer.texName(coords['y'])}$ [GeV]" )
    plt.title ( f"closure test, {sinjection} injection" )
    plt.savefig ( filename )
    from smodels_utils.plotting.mpkitty import timg
    timg ( filename )
    if args["interact"]:
        import sys, IPython; IPython.embed( colors = "neutral" ); sys.exit()

if __name__ == "__main__":
    path = "../data/fake_ewk1/"
    path = "../data/fake_ewkoff1/"
    minF = .1
    # path = "../data/fake_stops1/"
    interact = False
    args = { "path": path, "minF": minF, "interact": interact }
    plotClosure( args )
