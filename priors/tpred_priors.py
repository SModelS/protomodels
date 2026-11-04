#!/usr/bin/env python3
    
""" Collect here priors that are based on the theory predictions
"""

def tpredsBasedPrior ( protomodel ) -> tuple[float,str]:
    """ a prior that penalizes based on the theorypredictions """
    from tester.predictor import Predictor
    pr = Predictor ( "myprior", protomodel.environ )
    pr.obtainPredictions ( protomodel, keep_predictions = True )
    n_predictions = len ( pr.predictions )
    return 0, "not yet in use"
