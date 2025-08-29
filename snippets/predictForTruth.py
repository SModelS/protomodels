#!/usr/bin/env python3

""" simple code snippet that runs the true model, computes predictions, 
    test statistic, etc
"""

from builder.manipulator import Manipulator
# from builder.protomodel import ProtoModel
from tester.predictor import Predictor
from tester.critic import Critic
from smodels.experiment.databaseObj import Database

def predictForTruth():
    signal_model = "./signal_model.dict"
    dbpath = "./signal.pcl"
    print ( f"[predictForTruth] instantiate database {dbpath}" )
    db = Database ( dbpath )
    # protomodel = Protomodel()
    ma = Manipulator( signal_model )
    ma.M.dbpath = dbpath
    print ( f"[predictForTruth] starting" )
    predictor = Predictor( 0, db, do_srcombine=True )
    critic = Critic ( "critic", db, do_srcombine = True ) 
    cr, response = critic.predict_critic ( ma.M )
    print ( f"[predictForTruth] critic: {cr}, {response}" )
    predictor.predict ( ma, keep_predictions = True, force_computation_K = True )
    print ( f"[predictForTruth] predict K={ma.M.K} TL={ma.M.TL}" )
    import sys, IPython; IPython.embed( colors = "neutral" ); sys.exit()

if __name__ == "__main__":
    predictForTruth()
