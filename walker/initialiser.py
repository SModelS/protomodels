#!/usr/bin/env python3

""" A class that encapsulates the notion of starting from a sensible protomodel
"""

__all__ = [ "Initialiser" ]

import os, glob, sys, copy
import numpy as np
import scipy
import random
import pyslha
from base.loggerbase import LoggerBase
from typing import List, Set, Dict, Tuple, Union
from ptools.sparticleNames import SParticleNames
from builder.protomodel import ProtoModel
from builder.manipulator import Manipulator

LSP = ProtoModel.LSP

namer = SParticleNames ( susy = False )

def mergeTwoModels ( model1 : str, model2: str ) -> Union[None,Dict]:
    """ merge two models, add all particles from both models.
    If a particle appears in both models, its mass will be the average
    of the two, etc.

    :returns: merged model
    """
    if not os.path.exists ( model1 ):
        print ( f"[initialiser] {model1} does not exist" )
        return None
    if not os.path.exists ( model2 ):
        print ( f"[initialiser] {model2} does not exist" )
        return None
    f=open ( model1, "rt" )
    txt=f.read()
    f.close()
    dict1 = eval(txt)
    f=open ( model2, "rt" )
    txt=f.read()
    f.close()
    dict2 = eval(txt)
    import copy
    ret = copy.deepcopy(dict1)
    lsp1 = dict1["masses"][LSP]
    lsp2 = dict2["masses"][LSP]
    for pid,m in dict2["masses"].items():
        if not pid in ret["masses"]:
            # ok, we just copy <pid> from dict2
            ret["masses"][pid]=m-lsp2+lsp1
            ret["decays"][pid]=dict2["decays"][pid]
            for pidpair,ssms in dict2["ssmultipliers"].items():
                if pid in pidpair or -pid in pidpair:
                    ret["ssmultipliers"][pidpair]=ssms
        else:
            ret["masses"][pid]=.5*m + .5*ret["masses"][pid]
    for drop in [ "K", "TL", "xsecs[fb]", "walkerid", "step" ]:
        if drop in ret:
            ret.pop( drop )
    ret["timestamp"] = time.asctime()
    return ret

def mergeNModels ( models : List[Dict], add_timestamp : bool = True ) \
        -> Union[None,Dict]:
    """ merge two models, add all particles from both models.
    If a particle appears in both models, its mass will be the average
    of the two, etc.

    :returns: merged model
    """
    import time
    if len(models)== 0:
        print ( f"[initialiser] provided empty list of models" )
        return None
    if len(models)==1: # trivial merge
        if add_timestamp:
            models[0]["timestamp"] = time.asctime()
        return models[0]
    ret = { "masses": {}, "decays": {}, "ssmultipliers": {} }

    def collectPids ( models : List[Dict] ) -> Set:
        """ collect all pids in all models """
        pids = set()
        for m in models:
            for k in m["masses"].keys():
                pids.add ( k )
        return pids

    def computeAverageMassesForLSP ( pid : int, models : List[Dict] ) -> Dict:
        """ for the lsp we really average """
        masses=[]
        for model in models:
            if pid in model["masses"]:
                masses.append ( model["masses"][pid] )
        return float ( np.mean ( masses ) )

    def computeAverageMassesForPid ( pid : int, models : List[Dict] ) -> Dict:
        """ for the other particles we average over the distance to
        the LSP """
        masses=[]
        if pid == LSP:
            return computeAverageMassesForLSP ( pid, models )
        lspmasses = []
        for model in models:
            if pid in model["masses"]:
                masses.append ( model["masses"][pid] - model["masses"][LSP] )
                lspmasses.append ( model["masses"][LSP] )
        avg_delta = float ( np.mean ( masses ) )
        return float ( np.mean ( lspmasses ) ) + avg_delta

    def computeAverageDecaysForPid ( pid : int, models : List[Dict] ) -> Dict:
        decays = {}
        Stot = 0.
        for model in models:
            if pid in model["decays"]:
                mdecays = model["decays"][pid]
                for daughterpids, br in mdecays.items():
                    if not daughterpids in decays:
                        decays[daughterpids] = 0.
                    decays[daughterpids] += br
                    Stot += br
        for daughterpids, br in decays.items():
            decays[daughterpids]/=Stot
        return decays

    def computeAverageSSMs ( models : List[Dict] ) -> Dict:
        ssms = {}
        for model in models:
            mssms = model["ssmultipliers"]
            for k,v in mssms.items():
                ssms[k]=v # just bluntly take them over.
        return ssms

    # first, we collect all pids
    pids = collectPids ( models )
    # next we compute average masses
    masses, decays, ssms = {}, {}, {}
    for pid in pids:
        mass = computeAverageMassesForPid ( pid, models )
        masses[pid]=mass
        piddecays = computeAverageDecaysForPid ( pid, models )
        decays[pid] = piddecays
    ssms = computeAverageSSMs ( models )
    ret["masses"]=masses
    ret["decays"]=decays
    ret["ssmultipliers"]=ssms
    if add_timestamp:
        ret["timestamp"] = time.asctime()
    print ( f"[initialiser] merged {len(models)} models into:" )
    return ret


class Initialiser ( LoggerBase ):
    """ class to come up with a sensible first guess of a protomodel,
    from data. """
    cachefile = "pids.cache"

    def __init__ ( self, walkerid : Union[str,int] = 0,
            dictfile : str = "signal_database.dict",
            allowN1N1Prod : bool = True, verbose : bool = False ):
        """ constructor.

        :param walkerid: the walkerid we run this under
        :param dictfile: path to the database dict file we will base this on.
        dictfile is usally sth like signal_database.dict, *_database.dict, <dbver>.dict.
        :param allowN1N1Prod: do we allow N1 N1 production?
        """
        super ( Initialiser, self ).__init__ ( "ini" )
        self.printLogMessages = verbose
        dictfile = os.path.expanduser ( dictfile )
        self.dictfile = dictfile
        self.allowN1N1Prod = allowN1N1Prod
        from multiverse.expResModifier import readDictFile
        self.log ( f"reading database dict {dictfile}" )
        d = readDictFile ( dictfile )
        self.tempslha = "/dev/shm/temp.slha"
        self.meta = d["meta"]
        self.data = d["data"]
        # self.TLmax = 1. # disregard all results below this
        self.pmax = 0.25 # disregard all results above this
        self.computePDict()
        self.setMassRanges()
        self.ignore_pids = [ 1000002, 1000003, 1000004, 2000001,
                             2000002, 2000003, 2000004, 2000011,
                             2000013, 20000015, 1000014, 1000016, 
                             2000012, 20000014, 2000016 ]
        re = self.readInitialData()
        if not re:
            self.getTxParamsFromTemplates()

    def setMassRanges ( self ):
        """ set the mass ranges to draw from. for now set by hand.
        """
        self.massRanges = { 1000022: [ 50, 500 ] } # N1
        self.massRanges[1000023] = [60, 800 ] # N2
        self.massRanges[1000024] = [60, 800 ] # C1
        self.massRanges[1000005] = [60, 1500 ] # ~b1
        self.massRanges[2000005] = [300, 1500 ] # ~b2
        self.massRanges[1000006] = [100, 1800 ] # ~t
        self.massRanges[1000021] = [500, 3500 ] # ~g
        self.massRanges[1000011] = [100, 2000 ] # ~e
        self.massRanges[1000012] = [100, 2000 ] # ~nu_e
        self.massRanges[1000013] = [100, 2000 ] # ~mu
        self.massRanges[1000015] = [100, 2000 ] # ~tau
        squarkrange = [ 200, 1800 ]
        #squarks = [ 1000001, 2000001, 1000002, 2000002,
        #           1000003, 2000003, 1000004, 2000004 ]
        squarks = [ 1000001 ]

        for i in squarks:
            self.massRanges[i] = squarkrange # ~b

    def getTxParamsFor ( self, filename : str ):
        """ get pids, decays for slha template <filename>
        :param filename: e.g. ..../T1.template
        """
        self.log ( f"getTxParamsFor {filename}" )
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
            tmpfile.write ( f"{line}\n" )
        tmpfile.close()
        flatpids = set()
        for k,v in pids.items():
            for i in v:
                if i not in self.ignore_pids:
                    flatpids.add ( i )
        for k,v in pids.items():
            nv = set()
            for vi in v:
                if not vi in self.ignore_pids:
                    nv.add ( vi )
            self.pidsForTxnames[txname][k]=nv
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
            from base.locker import lock, unlock
            lock ( self.cachefile )
            with open ( self.cachefile, "wt" ) as f:
                f.write ( f"{self.pidsForTxnames}\n" )
                f.write ( f"{self.decaysForTxnames}\n" )
                f.write ( f"{self.ssmsForTxnames}\n" )
                f.close()
            unlock ( self.cachefile )

    def readInitialData ( self ) -> bool:
        """ read in all the data (pids,decays,ssms) from the slha files.
        :returns: False, if no cache file found.
        """
        if not os.path.exists ( self.cachefile ):
            self.error ( f"did not find {self.cachefile}" )
            return False
        self.log ( f"reading in all initial data from {self.cachefile}" )
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
        self.log ( f"get first file in {tarball}" )
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
        try:
            xsecs = r.xsections
            ssmpids = set()
            ignore_pids = self.ignore_pids + [ -x for x in self.ignore_pids ] + [ 2212]
            for k,v in xsecs.items():
                # print ( "k=", k, "  v=", v, "type v", type(v) )
                pids = tuple ( filter(lambda x: x not in ignore_pids, k) )
                # print ( "pids", pids )
                ssmpids.add ( pids )
            self.ssmsForTxnames[txname]=ssmpids
        except Exception as e:
            self.error ( "caught {e}: will skip for now" )

    def getTxParamsFromTemplates ( self ):
        """ get particle ids from template files in
        smodels-utils/slha/templates/ """
        self.log ( f"recreate {self.cachefile}" )
        pathname = "../../smodels-utils/slha/templates/"
        altpathname = "~/git/smodels-utils/slha/templates/"
        altpathname = os.path.expanduser ( altpathname )
        self.pidsForTxnames = {}
        self.decaysForTxnames = {}
        self.ssmsForTxnames = {}
        files = glob.glob ( f"{pathname}/T*.template" )
        files += glob.glob ( f"{altpathname}/T*.template" )
        if len(files)==0:
            self.error ( f"could not find template files!" )
            sys.exit()
        for f in files:
            self.getTxParamsFor ( f )

    def computePDict ( self ):
        """ compute the probabilities with which we choose a result
        """
        prels = {}
        ptot = 0.
        self.pvalues = {}
        for anaAndSRName,stats in self.data.items():
            #if not "TL" in stats:
            #    continue
            #TL = stats["TL"] # make sure we have unique TLs
            #if TL < self.TLmax: # we dont look at underfluctuations, or small TLs
            #    continue
            # Z = stats["new_Z"] # significance of the fake data (incl signal)
            p = stats["new_p"] # p-value of the fake data (incl signal)
            if p > self.pmax:
                continue
            ## FIXME for now we shoose by exp(Z), maybe
            ## we do sth better motivated
            # prel = np.exp ( Z )
            # prel = 1. / ( 1. - scipy.stats.norm.cdf ( np.sqrt(TL) ) )
            prel = 1. / p
            while prel in prels:
                prel+=1e-10
            value = stats
            value["id"]=anaAndSRName
            prels[ prel ] = value
            ptot += prel
            self.pvalues[anaAndSRName]={ "p": p }
        if len(prels)==0:
            self.error ( "computePDict: no results returned" )
        self.probs = dict( [ (k/ptot,v) for k,v in prels.items() ] )
        probkeys = list ( self.probs.keys() )
        probkeys.sort (reverse = True )
        self.probkeys = probkeys

    def randomlyChooseOneResult ( self ) -> Tuple[str,Dict]:
        """ randomly choose one result from self.probs

        :returns: tuple(txname, result-dictionary)
        result-dictionary is the whole result dictionary of the excess we are
        exploiting
        """
        txn = "TRV1"
        while txn in [ "TRV1", "TRS1" ]:
            # dont yet know how to handle these
            choice = np.random.choice(list(self.probs.values()),
                    1, p=list(self.probs.keys()) )
            result = choice[0]
            txns = result["txns"] # .split(",")
            if "," in txns:
                txns = txns.split(",")
            ## choose a random txname
            txn = txns
            ## FIXME this should be smarter:
            ## we should choose the one
            self.log ( f"randomly chosing among the txn {txn}" )
            self.log ( f"FIXME make this smarter, choose by efficiency!" )
            if type(txn) in [ list, tuple ]:
                txn  = str(np.random.choice ( txns ))
        self.log ( f"choosing random txn from {result['id']}: {txn}" )
        return txn, result

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
                if np.random.uniform ( 0, 1 ) < .5:
                    # tie them only in half the cases
                    return ties[apos]
                return None
        return None

    def getRandomMassesForTxname ( self, txname : str ) -> Dict:
        """ sample random mass values for the given txname """
        pidsdict = copy.deepcopy ( self.pidsForTxnames[txname] )
        masses = {}
        # masses[ProtoModel.LSP]=self.lspmass
        pid = ProtoModel.LSP
        lspmass = float(np.random.uniform ( *self.massRanges[pid] ))
        masses[pid]=lspmass
        self.log ( f"setting mass of LSP/{namer.asciiName(pid)} to {lspmass:.1f} -- I chose from [{self.massRanges[pid][0]:.1f},{self.massRanges[pid][1]:.1f}]" )
        #leftsquarks = [ 1000001, 1000002, 1000003, 1000004 ]
        leftsquarks = [ 1000001 ]
        #rightsquarks = [ 2000001, 2000002, 2000003, 2000004 ]
        rightsquarks = [  ]
        squarks = leftsquarks + rightsquarks
        mylightsquark = int(np.random.choice ( leftsquarks ))
        offshell = txname.endswith ( "off" ) or "ISR" in txname
        onshell = False
        if "tt" in txname and not "off" in txname:
            onshell = True
        if "W" in txname and not "off" in txname:
            onshell = True
        if "Z" in txname and not "off" in txname:
            onshell = True
        massgap = 0.
        if "W" in txname:
            massgap = 80.
        if "Z" in txname:
            massgap = 90.
        if "t" in txname:
            massgap = 170.

        for position,pids in pidsdict.items():
            if mylightsquark in pids:
                hasWarned = False
                for rm in squarks:
                    if rm == mylightsquark:
                        continue
                    if not hasWarned:
                        self.log ( f"there are many light quark-partners, will keep only {namer.asciiName(mylightsquark)}" )
                        hasWarned =True
                    if rm in pids:
                        pids.remove(rm)
        for position,pids in pidsdict.items():
            for pid in pids:
                if pid == ProtoModel.LSP:
                    continue
                if not pid in self.massRanges:
                    self.error ( f"we dont have mass ranges for pid={pid}({namer.asciiName(pid)}). will skip." )
                    continue
                    # sys.exit()
                mass = -1.
                massRanges = copy.deepcopy ( self.massRanges[pid] )
                # print ( f"orig massranges {massRanges}" )
                if lspmass > massRanges[0]:
                    massRanges[0] = lspmass + .5 # at least .5 gev gap
                if offshell:
                    massRanges[1] = float ( massRanges[0]+massgap )
                if onshell:
                    massRanges[0] = lspmass + massgap
                # print ( f"massranges {massRanges}" )
                mass = float(np.random.uniform ( *massRanges ))
                ## for C1 and N2: with a certain change we set them to the same
                ## value
                masses[pid]=mass
                self.log ( f"setting mass of {namer.asciiName(pid)} to {mass:.1f} -- I chose from [{self.massRanges[pid][0]:.1f},{self.massRanges[pid][1]:.1f}]" )
                apid = self.tiePids ( pid, pids )
                if apid != None:
                    self.log ( f"setting mass of {namer.asciiName(apid)} to {mass:.1f}" )
                    masses[apid]=mass

        # make sure the LSP is the lowest
        lspmass = masses[ProtoModel.LSP]
        for pid,mass in masses.items():
            if pid == ProtoModel.LSP:
                continue
            if mass < lspmass:
                masses[ProtoModel.LSP]=mass
                masses[pid]=lspmass

        if 1000023  in masses and 1000024 in masses:
            # with some probability we make them equal
            if random.uniform(0.,1.)<.1:
                masses[1000023] = masses[1000024]
            if random.uniform(0.,1.)>.9:
                masses[1000024] = masses[1000023]
        return masses

    def getAllowedParticles ( self, constraints ):
        allowed_particles = set()
        test_particles = { "W(": (24,), "Z(": (23,), "e": (11,),
            "mu": (13,), "tau": (15,), "l": (11,13), "L": (11,13,15),
            "W+(": ( 24,), "W-(": ( 24, ), "t(": (6,), "t+(": (6,),
            "t-(": (6,), "b": ( 5, ), "jet": ( 1,4 ), "q": ( 1,4 ),
            "nu": (12,) }
        for constraint in constraints:
            for particle,pids in test_particles.items():
                if particle in constraint:
                    for pid in pids:
                        allowed_particles.add ( pid )
        return allowed_particles

    def pidsAreGood ( self, pids : set , allowed_particles : set ):
        """ check if pids are what created the excess

        :returns: true if pids are in the list
        """
        # print ( f"@@X01checking pids {pids} against {allowed_particles}" )
        for pid in pids:
            if abs(pid) > 1000000: # throw out all BSM particles
                continue
            if abs(pid) in allowed_particles:
                continue
            return False
        return True

    def getPossibleChannelsFromConstraints ( self, txname : str,
            constraints : tuple[str] ) -> dict:
        """ for a given txname, alongside with the constraints of
        the result, return a dictionary of the decay channels we want to
        see being open, for multiple particles """
        all_channels = self.decaysForTxnames[txname]
        allowed_particles = self.getAllowedParticles ( constraints )

        good_channels = {}
        for mother,decays in all_channels.items():
            good_decays = decays
            if len(decays)>1: # no choice for == 1
                good_decays = set()
                for pids in decays:
                    if self.pidsAreGood ( pids, allowed_particles ):
                        good_decays.add ( pids )
            if len(good_decays)>0:
                good_channels[mother]=good_decays
        self.log ( f"for {txname}, {constraints}" )
        self.log ( f"we started with {all_channels}" )
        self.log ( f"we selected {good_channels}" )
        return good_channels

    def getDecaysForTxname ( self, txname : str, result : Dict ) -> Dict:
        """ get some random decays starting points
        :param result: the whole result dictionary of the excess we are
        """
        if not txname in self.decaysForTxnames:
            self.error ( f"we dont have any decays??" )
            sys.exit()
        decays = { ProtoModel.LSP: {} }

        if not "constraints" in result:
            print ( f"@@22 no constraint:  {result}" )
        constraints = result["constraints"]
        tmp = self.getPossibleChannelsFromConstraints ( txname, constraints )
        for mother,daughters in tmp.items():
            nbr_tot = 0.
            ndaughters = len(daughters)
            if not mother in decays:
                decays[mother]={}
            for daughterpids in daughters:
                #keys = []
                #for daughterpid in daughterpids:
                #    if daughterpid > 0:
                #        keys.append ( daughterpid )
                keys = tuple ( set ( [ abs(dp) for dp in daughterpids ] ) )
                # nbr = 1./ndaughters
                nbr = random.uniform(0.,1.)
                decays[mother][keys]=nbr
                nbr_tot += nbr
            for keys, nbr in decays[mother].items():
                decays[mother][keys]= nbr / nbr_tot
        self.log ( f"decays for {txname}: {decays}" )
        # self.log ( f"constraints were {result['constraints']}" )
        return decays

    def getSSMsForTxname ( self, txname : str ) -> Dict:
        """ get typical ssms for the given txname. """
        ssms = {}
        self.debug ( f"getSSMsForTxname for {txname}" )
        if not txname in self.ssmsForTxnames:
            return ssms
        for pids in self.ssmsForTxnames[txname]:
            if len(pids)<2:
                ## we currently ignore
                continue
            # ssm = 1.
            ssm = np.exp ( scipy.stats.norm.rvs() )
            ssms[pids]= float ( ssm )
        return ssms

    def getRandomSubmodelForTxname ( self, txname : str, result : dict ) -> Dict:
        """ given a txname, create a random submodel.
        :param txname: txname for which to create submodel
        :param result: the whole result dictionary of the excess we are
        exploiting
        """
        if not txname in self.pidsForTxnames:
            self.log ( f"we dont seem to have pids for {txname}" )
            return None
        masses = self.getRandomMassesForTxname ( txname )
        decays = self.getDecaysForTxname ( txname, result )
        ssms = self.getSSMsForTxname ( txname )
        model = { "masses": masses, "decays": decays, "ssmultipliers": ssms }
        return model

    def createRandomSubmodel ( self ) -> Dict:
        """ create a random submodel for one txname.
        we will merge later.

        :returns: model dict
        """
        ## choose a random txname
        txn, result = self.randomlyChooseOneResult()
        self.log ( f"creating random submodel for {txn}" )
        submodel = self.getRandomSubmodelForTxname ( txn, result )
        return submodel

    def propose ( self ):
        """ propose a random initial model. """
        # choose a random txn
        # self.getRandomMassForLSP()
        submodels = []
        nmodels = np.random.choice ( [1,2,3] )
        self.log ( f"proposed model will consist of {nmodels} submodels." )
        for i in range(nmodels):
            submodel = self.createRandomSubmodel()
            if submodel == None:
                self.warn ( f"got none as random submodel" )
                continue
            submodels.append ( submodel )
        self.submodels = submodels
        model = mergeNModels ( submodels )
        if self.allowN1N1Prod and model is not None:
            ssm = float ( np.exp ( scipy.stats.norm.rvs() ) )
            model["ssmultipliers"][(1000022,1000022)]=ssm
        self.log ( f"initialiser.propose proposes {model}" )
        return model

    def create ( self ) -> Manipulator:
        """ create the protomodel """
        dct = self.propose()
        ma = Manipulator ( dct )
        return ma

    def interact ( self, dbpath : os.PathLike  ):
        """ interactive shell, for debugging and development """
        from tester.predictor import Predictor
        pr = Predictor("ini", dbpath, do_srcombine = True )
        import IPython
        IPython.embed( colors = "neutral" )

    def readDBPath ( self, dbpath ):
        if dbpath != None:
            return dbpath
        if not os.path.exists ( "run.dict" ):
            return "official"
        with open ( "run.dict", "rt" ) as f:
            txt = f.read()
            d = eval ( txt )
            dbpath = d["dbpath"]
            self.log ( f"extracted dbpath {dbpath} from run.dict" )
            return dbpath

if __name__ == "__main__":
    import argparse
    argparser = argparse.ArgumentParser(
            description='CLI of initialiser' )
    argparser.add_argument ( '-d', '--dictfile',
            help='input database dict file ["signal_database.dict"]',
            type=str, default="*_database.dict" )
    argparser.add_argument ( '--dbpath',
            help='path to database, if none then read from run.dict [none]',
            type=str, default=None )
    argparser.add_argument ( '-v', '--verbose',
            help='verbose', action="store_true" )
    argparser.add_argument ( '-r', '--recompute_cache',
            help='verbose', action="store_true" )
    ## 310.dict is also a good default, for the actual observations
    args = argparser.parse_args()
    if args.recompute_cache and os.path.exists ( Initialiser.cachefile ):
        print ( f"[Initialiser] removing old cachefile" )
        os.unlink ( Initialiser.cachefile )
    if "*" in args.dictfile:
        import glob
        files = glob.glob ( args.dictfile )
        if len(files)==1:
            args.dictfile = files[0]
        elif len(files)==0:
            print ( f"[initialiser] could not find dict files with {args.dictfile}. specify!" )
            sys.exit(-1)
        else:
            print ( f"[initialiser] potential dict files are {files}. specify!")
            sys.exit(-1)
    allowN1N1Prod = True
    ini = Initialiser( "ini", args.dictfile, allowN1N1Prod = allowN1N1Prod,
            verbose = args.verbose )
    dbpath = ini.readDBPath ( args.dbpath )
    ini.interact( dbpath )
