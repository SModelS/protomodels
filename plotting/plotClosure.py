#!/usr/bin/env python3

from typing import Union
from math import exp

def getAllModels() -> list[dict]:
    """ get all the model dictionaries """
    dictfpattern = f"../data/dictfiles/pmodel_*.dict"
    import glob
    dictfiles = glob.glob ( dictfpattern )
    all_models = []
    for dictfile in dictfiles:
        with open ( dictfile, "rt" ) as f:
            models = eval ( f.read() )
            for model in models[500:]:
                if model["Accepted"]==0:
                    all_models.append ( model )
    return all_models

def getAllPModels() -> list[dict]:
    """ get all the model dictionaries """
    dictfpattern = f"../data/Pmodels/pmodel*.dict"
    import glob
    dictfiles = glob.glob ( dictfpattern )
    all_models = []
    for dictfile in dictfiles:
        with open ( dictfile, "rt" ) as f:
            model = eval ( f.read() )
            #    if model["Accepted"]==0:
            all_models.append ( model )
    return all_models


def getCoordinates ( models : Union[list,dict] ):
    plot={ "x": 1000006, "y": 1000022, "type_x": "mass", "type_y": "mass" }
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
                K = K**2 / 60.264
                # K /= 100.
                # K = exp(K) / exp(60.264) * 80
            Kvalues.append ( K )
    return xvalues, yvalues, Kvalues

def getTruthModel():
    truthfile = f"../data/truth.dict"
    with open ( truthfile, "rt" ) as f:
        model = eval ( f.read() )
        return model

def plotClosure():
    models = getAllPModels()
    truth = getTruthModel()
    x, y, K = getCoordinates ( models )
    x_true, y_true, K_true = getCoordinates ( truth )
    from matplotlib import pyplot as plt
    plt.scatter ( x, y, s = K )
    plt.scatter ( x_true, y_true, s=80, marker="+", color="red" )
    filename = "closure.png"
    plt.xlabel ( "mass, [Xt]" )
    plt.ylabel ( "mass, [X1Z]" )
    plt.title ( "closure test, stop injection " )
    plt.savefig ( filename )
    from smodels_utils.plotting.mpkitty import timg
    timg ( filename )

if __name__ == "__main__":
    plotClosure()
