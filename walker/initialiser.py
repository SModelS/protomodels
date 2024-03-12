#!/usr/bin/env python3

""" A class that encapsulates the notion of starting from a sensible protomodel
"""

__all__ = [ "Initialiser" ]

from builder.loggerbase import LoggerBase
from smodels.experiment.expResultObj import ExpResult
from typing import List, Set

class Initialiser ( LoggerBase ):
    """ class to come up with a sensible first guess of a protomodel,
    from data. """

    def __init__ ( self, listOfExpRes : List[ExpResult] ):
        """ constructor.

        :param listOfExpRes: a complete of experimental results.
        we will sift through them to build our candidate model.
        """
        self.listOfExpRes = listOfExpRes
        self.getSetOfTxNames()
        print ( txnames )

    def getListOfAnalyses ( self, txname : str ) -> List:
        """ given a txname, give us a list of all analyses that have a result
        for it.

        :param txname: a txname, e.g. "T1"
        :returns: list of analyses
        """

    def getSetOfTxNames( self )
        """ determine a set of all txnames in the database """
        ret=set()
        for er in self.listOfExpRes:
            for ds in er.datasets:
                for txn in ds.txnameList:
                    ret.add ( txn.txName )
        self.txnames = ret


if __name__ == "__main__":
    from smodels.experiment.databaseObj import Database
    db = Database( "official" )
    ini = Initialiser( db.getExpResults() )
