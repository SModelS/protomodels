#!/usr/bin/env python3

""" module to contain helper methods for the multiverse
"""

def createMyFile ( signal_model : str = "signal_model.dict",
        dbpath : str = "./signal.pcl", outfile : str = "my.truth",
        interactive : bool = False ):
    from builder.manipulator import Manipulator
    # from builder.protomodel import ProtoModel
    from tester.predictor import Predictor
    from tester.critic import Critic
    from smodels.experiment.databaseObj import Database
    print ( f"[predictForTruth] instantiate database {dbpath}" )
    db = Database ( dbpath )
    # protomodel = Protomodel()
    walkerid = "truth"
    ma = Manipulator( signal_model, walkerid = "truth" )
    ma.M.dbpath = dbpath
    print ( f"[predictForTruth] starting" )
    predictor = Predictor( "truth", db, do_srcombine=True )
    critic = Critic ( "truth", db, do_srcombine = True ) 
    cr, response = critic.predict_critic ( ma.M )
    print ( f"[predictForTruth] critic: {cr}, {response}" )
    predictor.predict ( ma, keep_predictions = True, force_computation_K = True )
    print ( f"[predictForTruth] predict K={ma.M.K} TL={ma.M.TL}" )
    ma = Manipulator( signal_model, walkerid = "truth" )
    ma.M.dbpath = dbpath
    ma.writeDictFile( outfile )
    print ( f"[predictForTruth] wrote truth into {outfile}" )
    if interactive:
        import sys, IPython; IPython.embed( colors = "neutral" ); sys.exit()

