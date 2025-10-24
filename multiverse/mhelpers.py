#!/usr/bin/env python3

""" module to contain helper methods for the multiverse
"""

import os
from typing import Union

def createMyFile ( signal_model : str = "signal_model.dict",
        dbpath : Union[None,os.PathLike] = "./signal.pcl", 
        outfile : str = "truth.dict",
        interactive : bool = False ):
    if signal_model == "":
        print ( f"[predictForTruth] no signal model defined." )
        return
    if not os.path.exists ( signal_model ):
        print ( f"[predictForTruth] {signal_model} does not exist." )
        return
    from builder.manipulator import Manipulator
    # from builder.protomodel import ProtoModel
    from tester.predictor import Predictor
    from tester.critic import Critic
    from base.runEnviron import RunEnviron
    from ptools.helpers import formatObject
    from smodels.experiment.databaseObj import Database
    environ = RunEnviron ( )
    if dbpath == None:
        dbpath = environ.dbpath
    print ( f"[predictForTruth] instantiate database {dbpath}" )
    db = Database ( dbpath )
    # protomodel = Protomodel()
    walkerid = "truth"
    ma = Manipulator( signal_model, walkerid = "truth", environ = environ )
    ma.M.dbpath = dbpath
    print ( f"[predictForTruth] starting" )
    predictor = Predictor( "truth", environ )
    critic = Critic ( "truth", environ ) 
    cr, response = critic.predict_critic ( ma.M )
    print ( f"[predictForTruth] critic: {cr}, {response}" )
    predictor.predict ( ma, keep_predictions = True, force_computation_K = True )
    print ( f"[predictForTruth] predict K={formatObject(ma.M.K)} TL={formatObject(ma.M.TL)}" )
    K,TL = ma.M.K, ma.M.TL
    ma = Manipulator( signal_model, walkerid = "truth", environ = environ )
    ma.M.dbpath = dbpath
    ma.M.K, ma.M.TL = K, TL
    ma.writeDictFile( outfile )
    print ( f"[predictForTruth] wrote truth into {outfile}" )
    if interactive:
        import sys, IPython; IPython.embed( colors = "neutral" ); sys.exit()

