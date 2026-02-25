#!/usr/bin/env python3

""" script that creates the data for the dilepton analysis """
    
import os

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

def loadPoint ( lr ):
    fname = f"points/pt{lr}.dict"
    if not os.path.exists ( fname ):
        return None
    with open ( fname, "rt" ) as f:
        d = eval (f.read() )
        f.close()
        return d

def create():
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
    from smodels.base.physicsUnits import fb, GeV, TeV
    points = np.arange ( .0, .15,.001)
    # points = np.arange ( .0, .15,.03)
    for lept_ratio in tqdm.tqdm ( points ):
        lr = float(lept_ratio)
        v = loadPoint ( lr )
        if v is not None:
            values[lr]=v
            continue
        protomodel.decays[1000023][(1000022,11,11)]=lr
        if lr == 0.:
            protomodel.decays[1000023].pop ( (1000022,11,11) )
        nu_ratio = 1. - 3. * lr
        protomodel.decays[1000023][(1000022,12,12)]=nu_ratio
        pr.predict ( ma, keep_predictions = True, force_computation_K = True )
        v = { "K": ma.M.K, "TL": ma.M.TL }
        robsmax, rexpmax = 0., 0.
        ## FIXME call ul_critic just like predict_critic would do!!!
        #allowed, n_sensitive, n_excluding = cr.ul_critic ( protomodel,
        #        pr.predictions, keep_predictions = True )
        # Run SModelS to get for UL-type predictions, and best SR preditcions if no UL-type result.
        slhafile = protomodel.createSLHAFile()
        sigmacut = 0.02*fb
        mingap = 10*GeV
        mingapISR = 1*GeV
        UL_preds, bestSR_preds = None, None
        rSM = cr.runSModelS( slhafile, combineSRs=False, ULpreds=True, sigmacut=sigmacut, mingap=mingap, mingapISR=mingapISR)
        if rSM not in ( None, [] ):
            UL_preds, bestSR_preds = rSM

        # Use best SR preds only if no UL-type result.
        predictions = cr.merge_preds(UL_preds,bestSR_preds)
        allowed, n_sensitive, n_excluding = cr.ul_critic(protomodel, predictions, keep_predictions=True )

        for d in cr.predictions:
            if d["robs"]>robsmax:
                robsmax = d["robs"]
            if d["rexp"]>rexpmax:
                rexpmax = d["rexp"]
            v[ f"{d['anaid']}:{d['dataid']}" ] = d
        v["robsmax"]=robsmax
        v["rexpmax"]=rexpmax
        v["allowed"]=allowed
        v["n_sensitive"]=n_sensitive
        v["n_excluding"]=n_excluding

        ## llhd-based critic
        predictions = cr.runSModelS( slhafile, combineSRs=True, ULpreds=False, sigmacut=sigmacut, mingap=mingap, mingapISR=mingapISR )
        allowed_by_llhd_critic, mostSensiComb, robsComb, rexpComb = cr.llhd_critic(predictions, cut=0.1, keep_predictions = False )
        v["llhd_allowed"]=allowed_by_llhd_critic
        v["robsComb"] = robsComb
        v["rexpComb"] = rexpComb
        writePoint ( lr, v )
        values[lr]=v
    write ( values )

if __name__ == "__main__":
    create()
