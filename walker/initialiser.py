#!/usr/bin/env python3

""" A class that encapsulates the notion of starting from a sensible protomodel
"""

__all__ = [ "Initialiser" ]

import os
import numpy as np
from builder.loggerbase import LoggerBase
# from ptools.helpers import computeZFromP
from typing import List, Set, Dict, Tuple
from colorama import Fore as ansi

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
        choice = np.random.choice(list(self.probs.values()), 
                1, p=list(self.probs.keys()) )
        return choice[0]

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
