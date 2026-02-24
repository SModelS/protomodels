#!/usr/bin/env python3

""" script that creates the data for the dilepton analysis """

def write ( values : dict ):
    """ write out the dictionary, add meta info """
    import time
    meta = { "created": time.asctime() }
    meta["comment"]="keys in dictionary are lepton BRs"
    values["meta"]=meta
    from ptools.helpers import py_dumps 
    d = py_dumps ( values )
    with open ( "dileptons.dict", "wt" ) as f:
        f.write ( d+"\n" )

def writePoint ( lr, values ):
    fname = f"points/pt{lr}.dict"
    with open ( fname, "wt" ) as f:
        f.write ( f"{values}\n" )
        f.close()

def create():
    import os
    if not os.path.exists ( "points" ):
        os.mkdir ( "points" )
    from base.runEnviron import RunEnviron
    environ = RunEnviron()
    walkerid = "dilepton"
    from ptools.hiscoreTools import fetchHiscoresObj
    infile = "hiscores_global.dict"
    hi = fetchHiscoresObj ( infile, None, environ = environ, walkerid = walkerid )
    protomodel = hi.hiscores[0]
    from builder.manipulator import Manipulator
    from tester.predictor import Predictor
    from tester.critic import Critic
    ma = Manipulator ( protomodel, environ )
    pr = Predictor( walkerid, environ = environ ) # instantiate for convenience
    cr = Critic ( walkerid, environ = environ )
    # print ( f"decays {protomodel.decays}" )
    import numpy as np
    values = {}
    import tqdm
    points = np.arange ( .0, .15,.001)
    points = np.arange ( .0, .15,.03)
    for lept_ratio in tqdm.tqdm ( points ):
        lr = float(lept_ratio)
        protomodel.decays[1000023][(1000022,11,11)]=lr
        if lr == 0.:
            protomodel.decays[1000023].pop ( (1000022,11,11) )
        nu_ratio = 1. - 3. * lr
        protomodel.decays[1000023][(1000022,12,12)]=nu_ratio
        pr.predict ( ma, keep_predictions = True, force_computation_K = True )
        v = { "K": ma.M.K, "TL": ma.M.TL }
        robsmax, rexpmax = 0., 0.
        allowed, n_sensitive, n_excluding = cr.ul_critic ( protomodel, 
                pr.predictions, keep_predictions = True )
        for d in cr.predictions:
            if d["robs"]>robsmax:
                robsmax = d["robs"]
            if d["rexp"]>rexpmax:
                rexpmax = d["rexp"]
            v[ d["anaid"]+":"+d["dataid"] ] = d
        v["robsmax"]=robsmax
        v["rexpmax"]=rexpmax
        v["allowed"]=allowed
        v["n_sensitive"]=n_sensitive
        v["n_excluding"]=n_excluding
        writePoint ( lr, v )
        values[lr]=v
    write ( values )

if __name__ == "__main__":
    create()
