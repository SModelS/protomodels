#!/usr/bin/env python3

""" A class that encapsulates the notion of starting from a sensible protomodel
"""

__all__ = [ "Initialiser" ]

import os, glob, time
import numpy as np
from builder.loggerbase import LoggerBase
from smodels.experiment.expResultObj import ExpResult
from typing import List, Set

class Initialiser ( LoggerBase ):
    """ class to come up with a sensible first guess of a protomodel,
    from data. """
    __slots__ = [ "dbpath", "db", "listOfExpRes", "llhdRatios", "txnames",
                  "picklefilename" ]

    def __init__ ( self, dbpath : str ):
        """ constructor.

        :param dbpath: path to database.
        """
        super ( Initialiser, self ).__init__ ( 0 )
        dbpath = os.path.expanduser ( dbpath )
        self.picklefilename = "initialiser.pcl"
        if os.path.exists ( self.picklefilename ):
            self.loadFromPickleFile()
        else:
            self.pprint ( f"initialising with {dbpath}" )
            self.dbpath = dbpath
            self.db = Database( dbpath )
            self.listOfExpRes = self.db.getExpResults()
            self.getDictOfTxNames()
            self.llhdRatios = {}
            self.getLlhdRatios()
            self.saveToPickleFile()
        self.interact()

    def getLlhdRatios ( self ):
        """ for all topo/analyses pairs, get the llhd ratios for all available
        mass points, and store in a huge dictionary """
        self.pprint ( f"computing llhd ratios" )
        for txname,analyses in self.txnames.items():
            # go to the respective validation folder, parse the validation dict 
            # files.
            for analysis in analyses:
                self.getLlhdRatiosFor ( txname, analysis )

    def interact ( self ):
        """ open an interactive shell """
        import IPython
        IPython.embed( colors = "neutral" )

    def saveToPickleFile ( self ):
        """ save all the information to a big pickle file """
        self.pprint ( f"saving to {self.picklefilename}" )
        d = { "llhdratios": self.llhdRatios, "dbpath": self.dbpath,
              "time": time.asctime(), "dbver": self.db.databaseVersion, }
        import pickle
        with open ( self.picklefilename, "wb" ) as f:
            pickle.dump ( d, f )
            f.close()

    def loadFromPickleFile ( self ):
        """ fetch all the info from a pickle file """
        import pickle
        with open ( "initialiser.pcl", "rb" ) as f:
            d = pickle.load ( f )
            self.llhdRatios = d["llhdratios"]
            self.dbpath = d["dbpath"]
            f.close()

    def getLlhdRatiosFor ( self, txname : str, analysis : str ):
        """ get the llhd ratios for the given txname/analysis pair """
        from smodels_utils.helper.various import getSqrts, getCollaboration
        sqrts = getSqrts ( analysis )
        collab = getCollaboration ( analysis )
        base = f"{self.dbpath}/{sqrts}TeV/{collab}/"
        path = f"{base}/{analysis}-eff/validation/{txname}*.py"
        valfiles = glob.glob ( path )
        points = {}
        for valfile in valfiles:
            with open(valfile) as f:
                exec(f.read(), globals())
                ## now we have validationData!!
                for point in validationData:
                    if not "axes" in point:
                        continue
                    if not "llhd" in point or point["llhd"] == 0.:
                        continue
                    if not "l_SM" in point:
                        # FIXME this we can compute post-mortem?
                        continue
                    ratio = point["l_SM"] / point["llhd"]
                    ## FIXME for more complicated cases this needs to change.
                    pa = point["axes"]
                    if "z" in pa:
                        axistuple = ( pa["x"], pa["y"], pa["z"] )
                    if "y" in pa:
                        axistuple = (pa["x"], pa["y"] )
                    else:
                        axistuple = (pa["x"],)
                    points[ axistuple ] =  -2 * np.log ( ratio )
        if len ( points ) == 0:
            return
        if not txname in self.llhdRatios:
            self.llhdRatios[txname]={}
        self.llhdRatios[txname][analysis]=points

    def getDictOfTxNames( self ):
        """ determine a dictionary of the database content with the txnames as keys, 
        and analyses names as values
        """
        self.debug ( f"create dictionary of txnames """ )
        ret={}
        for er in self.listOfExpRes:
            for ds in er.datasets:
                for txn in ds.txnameList:
                    if not txn.txName in ret:
                        ret[txn.txName]=set()
                    ret[txn.txName].add ( er.globalInfo.id )
        self.txnames = ret


if __name__ == "__main__":
    from smodels.experiment.databaseObj import Database
    dbpath = "~/git/smodels-database/"
    ini = Initialiser( dbpath )
