#!/usr/bin/env python3

""" simple code snippet that runs the true model, computes predictions, 
    test statistic, etc
"""

from multiverse.mhelpers import createMyFile

if __name__ == "__main__":
    createMyFile ( signal_model = "signal_model.dict",
        outfile = "truth.dict", interactive = False )
