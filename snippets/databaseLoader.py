#!/usr/bin/env python3

"""
.. module:: databaseLoader
   :synopsis: When running the complete test suite, we need to
              load the database only once

.. moduleauthor:: Wolfgang Waltenberger <wolfgang.waltenberger@gmail.com>

"""

import sys
sys.path.insert(0,"../")
from smodels.experiment.databaseObj import Database

dbpath = "../../smodels-database"
# dbpath = "./signal.pcl"
# dbpath = "unittest"

database = Database( dbpath)

if __name__ == "__main__":
    print ( "database", database )
    import IPython; IPython.embed( colors = "neutral" )
