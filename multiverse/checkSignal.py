#!/usr/bin/env python3

"""
.. module:: checkSignal
   :synopsis: simple snippet to check if signal was injected correctly

.. moduleauthor:: Wolfgang Waltenberger <wolfgang.waltenberger@gmail.com>

"""

import sys, os
sys.path.insert(0,"../")

def checkSignal():
    from smodels.experiment.databaseObj import Database
    # dbpath = "official"
    orig = Database( f"{os.environ['HOME']}/git/smodels-database"  )
    import sys, IPython; IPython.embed( colors = "neutral" )
    print ( orig )

if __name__ == "__main__":
    checkSignal()
