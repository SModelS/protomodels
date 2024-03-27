#!/usr/bin/env python3

""" A class that encapsulates the notion of starting from a sensible protomodel
"""

__all__ = [ "Initialiser" ]

import os, glob
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
        self.getParticleIdsFromTemplates()

    def getParticleIdsForTxname ( self, filename : str ):
        txname = filename.replace(".template","")
        pr = txname.rfind("/")
        txname = txname[pr+1:]
        f = open ( filename, "rt" )
        lines = f.readlines()
        f.close()
        ret = {}
        if not txname in self.pidsForTxnames:
            self.pidsForTxnames[txname]={}
        for line in lines:
            for x in [ 0, 1, 2 ]:
                if f"M{x}" in line:
                    tokens = line.split()
                    self.pidsForTxnames[txname][x]=int(tokens[0])
        if True:
            with open ( "pids.cache", "wt" ) as f:
                f.write ( self.pidsForTxnames+"\n" )
                f.close()

    def getParticleIdsFromTemplates ( self ):
        """ get particle ids from template files in 
        smodels-utils/slha/templates/ """
        pathname = "../../smodels-utils/slha/templates/"
        self.pidsForTxnames = {}
        files = glob.glob ( f"{pathname}/T*.template" )
        for f in files:
            self.getParticleIdsForTxname ( f )

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
        self.pprint ( f"we need to find random masses and decays for {txname}:{pids}" )

    def propose ( self ):
        """ propose a random initial model. """
        pass

    def interact ( self ):
        """ interactive shell, for debugging and development """
        import IPython
        IPython.embed( colors = "neutral" )

if __name__ == "__main__":
    dictfile = "../300.dict"
    ini = Initialiser( dictfile )
    ini.interact()
