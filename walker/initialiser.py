#!/usr/bin/env python3

""" A class that encapsulates the notion of starting from a sensible protomodel
"""

__all__ = [ "Initialiser" ]

import os, glob, time, scipy
import numpy as np
from builder.loggerbase import LoggerBase
from ptools.helpers import computeZFromP
from smodels.experiment.expResultObj import ExpResult
from typing import List, Set, Dict
from colorama import Fore as ansi
import matplotlib.pyplot as plt

# this block is simply to compute lSM!
from smodels_utils.helper.slhaManipulator import extractSLHAFileFromTarball
from smodels.matching.theoryPrediction import theoryPredictionsFor
from smodels.share.models.SMparticles import SMList
from smodels.particlesLoader import load
from smodels.base.model import Model
from smodels.base.physicsUnits import fb, GeV
from smodels.decomposition import decomposer

class Initialiser ( LoggerBase ):
    """ class to come up with a sensible first guess of a protomodel,
    from data. """

    def __init__ ( self, dbpath : str, ignore_pickle : bool = False ):
        """ constructor.

        :param ignore_pickle: if true, then ignore possible pickle files,
        recompute from scratch
        :param dbpath: path to database.
        """
        super ( Initialiser, self ).__init__ ( 0 )
        dbpath = os.path.expanduser ( dbpath )
        self.picklefilename = "initialiser.pcl"
        self.plotfname = "initialiser.png"
        self.significanceMap = {}
        self.analysesForTxName = {} # keep track of which analyses made it into
        # which txname
        if os.path.exists ( self.picklefilename ) and not ignore_pickle:
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
        self.produceSignificanceMaps()
        self.interact()

    def produceSignificanceMaps ( self ):
        """ produce the significance maps for all txnames,
        so we can later find maxima within these maps. """
        for txname,anaresults in self.llhdRatios.items():
            self.produceSignificanceMapForTxName ( txname, anaresults )

    def produceSignificanceMapForTxName ( self, txname : str, anaresults : Dict ):
        """ produce a significance map for a single txname

        :param txname: the txname, e.g. T1
        :param anaresults: dictionary with analsis name as key, and
        a dictionary with the point-wise results as values
        """
        smap = {}
        self.analysesForTxName[txname] = set()
        for ananame, points in anaresults.items():
            for point, Z in points.items():
                if not point in smap:
                    smap[point]={}
                smap[point][ananame]=Z
                self.analysesForTxName[txname].add ( ananame )
        self.significanceMap[txname] = smap


    def getLlhdRatios ( self ):
        """ for all topo/analyses pairs, get the llhd ratios for all available
        mass points, and store in a huge dictionary """
        self.pprint ( f"computing llhd ratios" )
        for txname,analyses in self.txnames.items():
            # go to the respective validation folder, parse the validation dict
            # files.
            for analysis in analyses:
                if analysis in [ "CMS-SUS-16-050", "ATLAS-SUSY-2018-22" ]:
                    continue
                self.getLlhdRatiosFor ( txname, analysis )

    def interact ( self ):
        """ open an interactive shell """
        import sys, os
        from importlib import reload
        import IPython
        import matplotlib.pyplot as plt
        print ( f"starting interfactive shell:" )
        print ( f"data members self.: llhdRatios, significanceMap, analysesForTxName" )
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

    def computeLSMFor ( self, analysis : str, txname : str, valfile : str,
                        point : Dict ) -> float:
        """ compute the SM likelihood for this point. it wasnt in the dict file.
        """
        slhafilename = point["slhafile"]
        print ( f"@@A computing for {analysis},{txname},{valfile},{slhafile}" )
        # print ( f"@@A point: {point}" )
        slhafile = extractSLHAFileFromTarball ( slhafilename )
        BSMList = load()
        model = Model(BSMparticles=BSMList, SMparticles=SMList)
        # ignorePQN = ['eCharge','colordim','spin']
        ignorePQN = [ ]
        model.updateParticles(inputFile=slhafile,
                              ignorePromptQNumbers = ignorePQN )
        sigmacut = 0.001*fb
        mingap = 5.*GeV
        topDict = decomposer.decompose(model, sigmacut,
                  massCompress=True, invisibleCompress=True, minmassgap=mingap)
        # self.db.selectExpResults ( )
        # analyses = [ analysis, analysis+"-agg" ]
        analyses = [ analysis ]
        # datasets = [ "all" ]
        dataset = point["dataset"]
        if dataset == "(combined)": ## combined results
            dataset = "all"
            self.db.selectExpResults ( analysisIDs = analyses,
                                       datasetIDs = [ dataset ],
                                       dataTypes = [ "efficiencyMap" ] )
            preds = theoryPredictionsFor( self.db, topDict, useBestDataset = True,
                                          combinedResults=True )
        else: ## sets of SRs
            self.db.selectExpResults ( analysisIDs = analyses,
                                       datasetIDs = [ dataset ],
                                       dataTypes = [ "efficiencyMap" ] )
            preds = theoryPredictionsFor( self.db, topDict, useBestDataset = True,
                                          combinedResults=False )
        if len(preds)==1:
            os.unlink ( slhafile )
            ret = preds[0].lsm()
            print ( f"@@C returning {ret:.2g} for {analysis}:{txname}:{dataset}" )
            return ret
        #preds = theoryPredictionsFor( self.db, topDict, useBestDataset = True,
        #                              combinedResults=False )
        print ( f"@@A something is off, got {len(preds)} results:" )
        import sys, IPython; IPython.embed( colors = "neutral" ); sys.exit()

        os.unlink ( slhafile )
        return float("nan")

    def getLlhdRatiosFor ( self, txname : str, analysis : str ):
        """ get the llhd ratios for the given txname/analysis pair """
        from smodels_utils.helper.various import getSqrts, getCollaboration
        sqrts = getSqrts ( analysis )
        collab = getCollaboration ( analysis )
        base = f"{self.dbpath}/{sqrts}TeV/{collab}/"
        # first we try the combined val files
        path = f"{base}/{analysis}-eff/validation/{txname}*_combined.py"
        valfiles = glob.glob ( path )
        if len(valfiles)==0: ## seems like there are no combined ones
            # go for the individual ones
            path = f"{base}/{analysis}-eff/validation/{txname}*.py"
            valfiles = glob.glob ( path )
            return # FIXME for now
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
                        lSM = self.computeLSMFor ( analysis, txname, valfile, point )
                        # FIXME this we can compute post-mortem?
                        point["l_SM"] = lSM
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

    def significanceFromT ( self, T : float ) -> float:
        """ compute the significance from T """
        p = scipy.stats.chi2.cdf ( T, df = 1 )
        Z = computeZFromP ( p )
        return Z

    def plot ( self, txname : str ):
        """ plot for txname """
        plt.clf()
        anas = self.analysesForTxName[txname]
        smap = self.significanceMap[txname]
        ## first determine the grid
        x, y = set(), set()
        for point, anaresults in smap.items():
            x.add ( point[0] )
            y.add ( point[1] )
        x, y = list ( x ), list (y )
        x.sort(); y.sort()
        X,Y = np.meshgrid(x,y) # we have a mesh grid, now fill Z
        Z =[]
        ## store also the position of the maximum
        maxx,maxy,maxvalue = float("nan"),float("nan"),-1e9
        for yi in range(len(y)): # index for y
            row = []
            for xi in range(len(x)): # index for x
                masses = (x[xi],y[yi])
                value = float("nan")
                if masses in smap:
                    if True: # len(smap[masses])==len(anas):
                        value = sum ( smap[masses].values() )
                        value = self.significanceFromT ( value  )
                if np.isinf ( value ):
                    value = float("nan")
                if value > maxvalue:
                    maxvalue = value
                    maxx, maxy = masses
                row.append ( value )
            Z.append ( row )
        Z = np.array ( Z )
        # vmin, vmax = -200, 200
        vmin, vmax = np.nanmin(Z), np.nanmax(Z)
        self.pprint ( f"plotting values between {vmin:.2g} and {vmax:.2g}" )
        plt.pcolormesh(X,Y,Z,shading="nearest", vmin=vmin, vmax=vmax)
        plt.scatter ( [ maxx ], [ maxy ], s=60, c="black", marker="*" )
        plt.colorbar()
        plt.title ( f"T, {txname}" )
        plt.xlabel ( f"m(x)" )
        plt.ylabel ( f"m(y)" )
        plt.savefig ( self.plotfname )
        self.show()

    def show ( self ) -> bool:
        """ show plot.

        :returns: true, if successful
        """
        import shutil
        if shutil.which ( "timg" ) is None:
            return False
        cmd = f"timg {self.plotfname}"
        import subprocess
        o = subprocess.getoutput ( cmd )
        print ( o )
        return True


if __name__ == "__main__":
    from smodels.experiment.databaseObj import Database
    dbpath = "~/git/smodels-database/"
    ini = Initialiser( dbpath, ignore_pickle = False )
