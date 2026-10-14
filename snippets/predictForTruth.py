#!/usr/bin/env python3

""" simple code snippet that runs the true model, computes predictions, 
    test statistic, etc
"""

from multiverse.mhelpers import createMyTruthFile

if __name__ == "__main__":
    createMyTruthFile ( signal_model = "signal_model.dict",
        outfile = "truth.dict", interactive = False )
