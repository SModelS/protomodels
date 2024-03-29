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
from ptools.sparticleNames import SParticleNames
from builder.protomodel import ProtoModel

namer = SParticleNames ( susy = False )

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
        self.massRanges[1000005] = [50, 1200 ] # ~b
        self.massRanges[1000006] = [100, 1400 ] # ~t
        self.massRanges[1000021] = [1000, 3000 ] # ~t
        squarkrange = [ 200, 1800 ]
        for i in [ 1000001, 2000001 ]:
            self.massRanges[i] = squarkrange # ~b

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
        pids = {}
        for line in lines:
            p1 = line.find("#")
            if p1 > -1:
                line = line[:p1]
            if line.endswith("\n"):
                line = line[:-1]
            for x in [ 0, 1, 2 ]:
                if f"M{x}" in line or f"m{x}" in line:
                    tokens = line.split()
                    if not x in pids:
                        pids[x]=set()
                    pids[x].add ( int(tokens[0]) )
                    line = line.replace ( f"M{x}", "100" )
                    line = line.replace ( f"m{x}", "100" )
            tmpfile.write ( line+ "\n" )
        tmpfile.close()
        flatpids = set()
        for k,v in pids.items():
            for i in v:
                flatpids.add ( i )
        for k,v in pids.items():
            self.pidsForTxnames[txname][k]=v
        import pyslha
        r = pyslha.readSLHAFile(tempf)
        for pid in flatpids:
            if pid == ProtoModel.LSP:
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
            self.pprint ( f"choosing random txn from {result['id']}: {txn}" )
        return txn

    def getRandomMassesForTxname ( self, txname : str ) -> Dict:
        """ sample random mass values for the given txname """
        pidsdict = self.pidsForTxnames[txname]
        masses = {}
        pid = ProtoModel.LSP
        lspmass = random.uniform ( *self.massRanges[pid] )
        self.pprint ( f"setting mass of {namer.asciiName(pid)} to {lspmass:.1f}" )
        masses[pid]=lspmass

        for position,pids in pidsdict.items():
            for pid in pids:
                if pid == ProtoModel.LSP:
                    continue
                if not pid in self.massRanges:
                    self.error ( f"we dont have mass ranges for pid={pid}({namer.asciiName(pid)})" )
                    sys.exit()
                mass = -1.
                while mass < lspmass:
                    mass = random.uniform ( *self.massRanges[pid] )
                ## for C1 and N2: with a certain change we set them to the same 
                ## value
                if pid in [ 1000023, 1000024 ] and random.uniform(0,1)<.5:
                    self.pprint ( f"setting mass of {namer.asciiName(1000023)} and {namer.asciiName(1000024)} to {mass:.1f}" )
                    masses[1000023]=mass
                    masses[1000024]=mass
                else:
                    self.pprint ( f"setting mass of {namer.asciiName(pid)} to {mass:.1f}" )
                    masses[pid]=mass
        return masses

    def getDecaysForTxname ( self, txname : str ) -> Dict:
        """ get some random decays starting points """
        if not txname in self.decaysForTxnames:
            self.error ( f"we dont have any decays??" )
            sys.exit()
        decays = { ProtoModel.LSP: {} }
        tmp = self.decaysForTxnames[txname]
        for mother,daughters in tmp.items():
            if not mother in decays:
                decays[mother]={}
            for daughterpids in daughters:
                keys = []
                for daughterpid in daughterpids:
                    if daughterpid > 0:
                        keys.append ( daughterpid )
                decays[mother][tuple(keys)]=1.0
        return decays

    def getRandomSubmodelForTxname ( self, txname : str ) -> Dict:
        """ given a txname, create a random submodel. """
        if not txname in self.pidsForTxnames:
            self.pprint ( "we dont seem to have pids for {txname}" )
            return None
        masses = self.getRandomMassesForTxname ( txname )
        decays = self.getDecaysForTxname ( txname )
        ssms = {}
        model = { "masses": masses, "decays": decays, "ssms": ssms }
        return model

    def createRandomSubmodel ( self ) -> Dict:
        """ create a random submodel for one txname.
        we will merge later.
        """
        ## choose a random txname
        txn = self.randomlyChooseOneResult()
        self.pprint ( f"creating random submodel for {txn}" )
        submodel = self.getRandomSubmodelForTxname ( txn )
        return submodel

    def propose ( self ):
        """ propose a random initial model. """
        # choose a random txn
        submodel = self.createRandomSubmodel()
        return submodel

    def mergeSubmodels ( self, submodels : List ) -> Dict:
        """ merge several submodels to a model.
        """
        from ptools.hiscoreTools import mergeTwoModels
        return mergeTwoModels ( *submodels )
        
    def interact ( self ):
        """ interactive shell, for debugging and development """
        import IPython
        IPython.embed( colors = "neutral" )

if __name__ == "__main__":
    dictfile = "../300.dict"
    ini = Initialiser( dictfile )
    ini.interact()
