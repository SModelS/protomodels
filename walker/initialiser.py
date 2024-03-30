#!/usr/bin/env python3

""" A class that encapsulates the notion of starting from a sensible protomodel
"""

__all__ = [ "Initialiser" ]

import os, glob, sys
import numpy as np
import scipy
import random
import pyslha
from builder.loggerbase import LoggerBase
from typing import List, Set, Dict, Tuple, Union
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
        self.tempslha = "/dev/shm/temp.slha"
        self.cachefile = "pids.cache"
        self.meta = d["meta"]
        self.data = d["data"]
        self.Zmax = 1. # disregard all results below this
        self.computePDict()
        self.setMassRanges()
        re = self.readInitialData()
        if not re:
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
        for i in [ 1000001, 2000001, 1000002, 2000002,
                   1000003, 2000003, 1000004, 2000004 ]:
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
        if not txname in self.ssmsForTxnames:
            self.ssmsForTxnames[txname]={}
        tmpfile = open ( self.tempslha, "wt" )
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
        r = pyslha.readSLHAFile(self.tempslha)
        for pid in flatpids:
            if pid == ProtoModel.LSP:
                continue
            decays = r.decays[pid].decays
            if not pid in self.decaysForTxnames[txname] and len(decays)>0:
                self.decaysForTxnames[txname][pid]=set()
            for decay in decays:
                ids = tuple ( decay.ids )
                self.decaysForTxnames[txname][pid].add(ids)
        os.unlink ( self.tempslha )
        self.getDefaultSSMs ( filename )

        if True:
            with open ( self.cachefile, "wt" ) as f:
                f.write ( f"{self.pidsForTxnames}\n" )
                f.write ( f"{self.decaysForTxnames}\n" )
                f.write ( f"{self.ssmsForTxnames}\n" )
                f.close()

    def readInitialData ( self ) -> bool:
        """ read in all the data (pids,decays,ssms) from the slha files. 
        :returns: False, if no cache file found.
        """
        if not os.path.exists ( self.cachefile ):
            self.error ( f"did not find {self.cachefile}" )
            return False
        self.pprint ( f"reading in all initial data from {self.cachefile}" )
        with open ( self.cachefile, "rt" ) as f:
            lines = f.readlines()
            f.close()
            self.pidsForTxnames = eval(lines[0])
            self.decaysForTxnames = eval(lines[1])
            self.ssmsForTxnames = eval(lines[2])
        return True

    def getDefaultSSMs ( self, templatename : str ):
        """ get default ssms for templatename
        :param templatename: e.g. ../../smodels-utils/slha/templa
        """
        txname = templatename.replace(".template","")
        pr = txname.rfind("/")
        txname = txname[pr+1:]
        tarball = templatename.replace(".template",".tar.gz").replace("templates/","")
        if not os.path.exists ( tarball ):
            self.debug ( f"cannot find {tarball}, cannot get default productions." )
            return
        self.pprint ( f"get first file in {tarball}" )
        import tarfile
        tar = tarfile.open ( tarball, "r:gz" )
        files = tar.members
        fobj = tar.extractfile ( files[0].name )
        txt = fobj.read() 
        with open ( self.tempslha, "wt" ) as f:
            f.write ( txt.decode("ascii") )
            f.close()
        tar.close()
        r = pyslha.readSLHAFile(self.tempslha)
        xsecs = r.xsections
        ssmpids = set()
        for k,v in xsecs.items():
            pids = tuple ( filter(lambda x: x not in [ 2212 ], k) )
            ssmpids.add ( pids )
        self.ssmsForTxnames[txname]=ssmpids

    def getTxParamsFromTemplates ( self ):
        """ get particle ids from template files in 
        smodels-utils/slha/templates/ """
        pathname = "../../smodels-utils/slha/templates/"
        self.pidsForTxnames = {}
        self.decaysForTxnames = {}
        self.ssmsForTxnames = {}
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
            # prel = np.exp ( Z )
            prel = 1. / ( 1. - scipy.stats.norm.cdf ( Z ) )
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

    def tiePids ( self, pid : int, pids : List[int] ) -> Union[None,int]:
        """ determine if we tie this pid to another pid, meaning
        we set the mass of another pid to the value of this pid.
        :param pid: check for this pid
        :param pids: these are all the pids that are there. if the 
        alternative pid is not in pids, dont tie

        :returns: None if we dont tie pids, pid of other particle if yes
        """
        ties = [ ( 1000023, 1000024 ) ]
        for tie in ties:
            if pid in ties:
                pos = ties.index ( pid )# pos is our guy
                apos = 1 - pos # thats the alternative pid
                if not ties[apos] in pids:
                    return None
                if random.uniform ( 0, 1 ) < .5:
                    # tie them only in half the cases
                    return ties[apos]
                return None
        return None

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
                masses[pid]=mass
                self.pprint ( f"setting mass of {namer.asciiName(pid)} to {mass:.1f}" )
                apid = self.tiePids ( pid, pids )
                if apid != None:
                    self.pprint ( f"setting mass of {namer.asciiName(apid)} to {mass:.1f}" )
                    masses[apid]=mass
        return masses

    def getDecaysForTxname ( self, txname : str ) -> Dict:
        """ get some random decays starting points """
        if not txname in self.decaysForTxnames:
            self.error ( f"we dont have any decays??" )
            sys.exit()
        decays = { ProtoModel.LSP: {} }
        tmp = self.decaysForTxnames[txname]
        for mother,daughters in tmp.items():
            ndaughters = len(daughters)
            if not mother in decays:
                decays[mother]={}
            for daughterpids in daughters:
                keys = []
                for daughterpid in daughterpids:
                    if daughterpid > 0:
                        keys.append ( daughterpid )
                decays[mother][tuple(keys)]=1./ndaughters
        self.pprint ( f"decays for {txname}: {decays}" )
        return decays

    def getSSMsForTxname ( self, txname : str ) -> Dict:
        """ get typical ssms for the given txname. """
        ssms = {}
        self.debug ( f"getSSMsForTxname for {txname}" )
        if not txname in self.ssmsForTxnames:
            return ssms
        for pids in self.ssmsForTxnames[txname]:
            ssms[pids]=1.
        return ssms

    def getRandomSubmodelForTxname ( self, txname : str ) -> Dict:
        """ given a txname, create a random submodel. """
        if not txname in self.pidsForTxnames:
            self.pprint ( "we dont seem to have pids for {txname}" )
            return None
        masses = self.getRandomMassesForTxname ( txname )
        decays = self.getDecaysForTxname ( txname )
        ssms = self.getSSMsForTxname ( txname )
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
        submodels = []
        nmodels = random.choice ( [1,2,3] )
        self.pprint ( f"proposed model will consist of {nmodels} submodels." )
        for i in range(nmodels):
            submodels.append ( self.createRandomSubmodel() )
        self.submodels = submodels
        from ptools.hiscoreTools import mergeNModels
        model = mergeNModels ( submodels )
        return model

    def interact ( self ):
        """ interactive shell, for debugging and development """
        import IPython
        IPython.embed( colors = "neutral" )

if __name__ == "__main__":
    dictfile = "../300.dict"
    ini = Initialiser( dictfile )
    ini.interact()
