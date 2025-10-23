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

def getCoordinates ( models : Union[list,dict], maxK : float, minF : float ):
    plot={ "x": 1000006, "y": 1000022, "type_x": "mass", "type_y": "mass" }
    plot={ "x": 1000023, "y": 1000022, "type_x": "mass", "type_y": "mass" }
    xvalues, yvalues, Kvalues = [], [], []
    xpid, ypid = plot["x"], plot["y"]
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

def plotClosure():
    from matplotlib import pyplot as plt
    path = "../data/fake_ewk1/"
    # path = "../data/fake_stops1/"
    # models = getAllPModels( path )
    minF = .7
    models, maxK = getAllModels( path, minF )
    truth = getTruthModel( path )
    x_true, y_true, K_true = getCoordinates ( truth, maxK, 0. )
    splitm = splitByDictFile ( models )
    colors = plt.cm.viridis(np.linspace(0.3, 0.9, len(splitm)))
    for i,(df,models) in enumerate(splitm.items()):
        x, y, K = getCoordinates ( models, maxK, minF )
        plt.scatter ( x, y, s = K, alpha=0.5, c = colors[i],
                      label=df, edgecolors = darken(colors[i]) )
    plt.scatter ( x_true, y_true, s=140, marker="+", color="black", 
                  label="truth" )
    plt.legend()
    filename = "closure.png"
    plt.xlabel ( "mass, Xt [GeV]" )
    plt.ylabel ( "mass, X1Z [GeV]" )
    plt.title ( "closure test, stop injection" )
    plt.savefig ( filename )
    from smodels_utils.plotting.mpkitty import timg
    timg ( filename )

if __name__ == "__main__":
    plotClosure()
