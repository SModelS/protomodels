#!/usr/bin/env python3

""" A class that encapsulates the notion of starting from a sensible protomodel
"""

__all__ = [ "Initialiser" ]

import os, glob, sys
import numpy as np
import random
from builder.loggerbase import LoggerBase
# from ptools.helpers import computeZFromP
from typing import List, Set, Dict, Tuple
# from colorama import Fore as ansi

class Initialiser ( LoggerBase ):
    """ class to come up with a sensible first guess of a protomodel,
    from data. """

    def __init__ ( self, dictfile : str ):
        """ constructor.

        :param dictfile: path to the database dict file we will base this on
        """
        super ( Initialiser, self ).__init__ ( 0 )
        dictfile = os.path.expanduser ( dictfile )
        self.dictfile = dictfile
        from ptools.expResModifier import readDictFile
        d = readDictFile ( dictfile )
        self.meta = d["meta"]
        self.data = d["data"]
        self.Zmax = 1. # disregard all results below this
        self.computePDict()
        self.setMassRanges()
        self.getTxParamsFromTemplates()

    def setMassRanges ( self ):
        """ set the mass ranges to draw from. for now set by hand. 
        """
        self.massRanges = { 1000022: [ 0, 500 ] } # N1
        self.massRanges[1000023] = [50, 800 ] # N2
        self.massRanges[1000024] = [50, 800 ] # C1

    def getTxParamsFor ( self, filename : str ):
        """ get pids, decays for slha template <filename>
        :param filename: e.g. ..../T1.template
        """
        txname = filename.replace(".template","")
        pr = txname.rfind("/")
        txname = txname[pr+1:]
        f = open ( filename, "rt" )
        lines = f.readlines()
        f.close()
        ret = {}
        if not txname in self.pidsForTxnames:
            self.pidsForTxnames[txname]={}
        if not txname in self.decaysForTxnames:
            self.decaysForTxnames[txname]={}
        tempf = "/dev/shm/temp.slha"
        tmpfile = open ( tempf, "wt" )
        pids = set()
        for line in lines:
            p1 = line.find("#")
            if p1 > 0:
                line = line[:p1]
            for x in [ 0, 1, 2 ]:
                if f"M{x}" in line or f"m{x}" in line:
                    tokens = line.split()
                    pids.add ( int(tokens[0]) )
                    line = line.replace ( f"M{x}", "100" )
            tmpfile.write ( line )
        tmpfile.close()
        self.pidsForTxnames[txname][x]=pids
        import pyslha
        r = pyslha.readSLHAFile(tempf)
        for pid in pids:
            if pid == 1000022:
                continue
            decays = r.decays[pid].decays
            if not pid in self.decaysForTxnames[txname] and len(decays)>0:
                self.decaysForTxnames[txname][pid]=set()
            for decay in decays:
                ids = tuple ( decay.ids )
                self.decaysForTxnames[txname][pid].add(ids)
        os.unlink ( tempf )
        if True:
            with open ( "pids.cache", "wt" ) as f:
                f.write ( f"{self.pidsForTxnames}\n" )
                f.write ( f"{self.decaysForTxnames}\n" )
                f.close()

    def getTxParamsFromTemplates ( self ):
        """ get particle ids from template files in 
        smodels-utils/slha/templates/ """
        pathname = "../../smodels-utils/slha/templates/"
        self.pidsForTxnames = {}
        self.decaysForTxnames = {}
        files = glob.glob ( f"{pathname}/T*.template" )
        for f in files:
            self.getTxParamsFor ( f )

    def computePDict ( self ):
        """ compute the probabilities with which we choose a result """
        prels = {}
        ptot = 0.
        for anaAndSRName,stats in self.data.items():
            if not "Z" in stats:
                continue
            Z = stats["Z"] # make sure we have unique Zs
            if Z < self.Zmax: # we dont look at underfluctuations, or small Zs
                continue
            ## FIXME for now we shoose by exp(Z), maybe
            ## we do sth better motivated
            prel = np.exp ( Z )
            while prel in prels:
                prel+=1e-10
            value = stats
            value["id"]=anaAndSRName
            prels[ prel ] = value
            ptot += prel
        self.probs = dict( [ (k/ptot,v) for k,v in prels.items() ] )
        probkeys = list ( self.probs.keys() )
        probkeys.sort (reverse = True )
        self.probkeys = probkeys

    def randomlyChooseOneResult ( self ):
        """ randomly choose one result from self.probs """
        txn = "TRV1"
        while txn in [ "TRV1", "TRS1" ]: 
            # dont yet know how to handle these
            choice = np.random.choice(list(self.probs.values()), 
                    1, p=list(self.probs.keys()) )
            result = choice[0]
            txns = result["txns"].split(",")
            ## choose a random txname
            txn  = random.choice ( txns )
        return txn

    def getRandomSubmodelForTxname ( self, txname : str ):
        """ given a txname, create a random submodel. """
        if not txname in self.pidsForTxnames:
            self.pprint ( "we dont seem to have pids for {txname}" )
            return None
        pids = self.pidsForTxnames[txname]
        masses = {}
        for pid in pids:
            if not pid in self.massRanges:
                self.error ( f"we dont have mass ranges for {pid}" )
                sys.exit()
            mass = random.uniform ( *self.massRanges[pid] )
            self.pprint ( f"setting mass of {pid} to {mass}" )
            masses[pid]=mass

    def createRandomSubmodel ( self ) -> Dict:
        """ create a random submodel for one txname.
        we will merge later.
        """
        ## choose a random txname
        txn = self.randomlyChooseOneResult()
        self.getRandomSubmodelForTxname ( txn )

    def propose ( self ):
        """ propose a random initial model. """
        # choose a random txn
        self.createRandomSubmodel()

    def interact ( self ):
        """ interactive shell, for debugging and development """
        import IPython
        IPython.embed( colors = "neutral" )

if __name__ == "__main__":
    dictfile = "../300.dict"
    ini = Initialiser( dictfile )
    ini.interact()
