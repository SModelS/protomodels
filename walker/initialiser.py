#!/usr/bin/env python3

""" A class that encapsulates the notion of starting from a sensible protomodel
"""

__all__ = [ "Initialiser" ]

import os, glob, sys, copy
import numpy as np
import scipy
import random
import pyslha
from typing import Union, Dict, List, Set
from functools import lru_cache

from smodels_utils.helper.terminalcolors import *

from base.loggerbase import LoggerBase
from base.runEnviron import RunEnviron
from base.locker import lock, unlock

from ptools.sparticleNames import SParticleNames
from ptools.helpers import formatObject
from ptools.refxsecComputer import RefXSecComputer

from builder.protomodel import ProtoModel
from builder.manipulator import Manipulator
from tester.predictor import Predictor
from tester.critic import Critic

LSP = ProtoModel.LSP

namer = SParticleNames ( susy = False )

class TopNDict(dict):
    """ a dictionary that only keeps the top 20 entries """
    def __init__ ( self, nmax : int = 20 ):
        super ( TopNDict, self ).__init__()
        self.nmax = nmax
    def __setitem__(self, key, value):
        # insert the key-value pair
        super().__setitem__(key, value)
        # if we now have more than 10 entries, drop the smallest key
        if len(self) > self.nmax:
            smallest_key = min(self.keys())
            del self[smallest_key]
    def update(self, d : dict ):
        for k,v in d.items():
            self[k]=v

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

class Initialiser ( LoggerBase ):
    """ class to come up with a sensible first guess of a protomodel,
    from data. """
    cachefile = "pids.cache"
    effscachefile = "xsecs.cache"

    def __init__ ( self, walkerid : Union[str,int],
            dictfile : os.PathLike, environ : RunEnviron, 
            verbose : bool = False ):
        """ constructor.

        :param walkerid: the walkerid we run this under
        :param environ: a run environment
        :param verbose: true if verbose
        """
        super ( Initialiser, self ).__init__ ( walkerid )
        self.dictfile = dictfile
        self.printLogMessages = verbose
        self.environ = environ
        self.warnings= {}
        self._pm = ProtoModel ( walkerid, environ = self.environ )
        self.dsIdsInProposal = set()
        self.xsecComputer = RefXSecComputer( allowN1N1Prod = self.environ.allowN1N1Prod )
        self.mapTxnames = { "TRS1": None, "TChiWWISRqq": "TChiWWoff",
            "TRV1": None, "TDTM1F": None, "TDTM2F": None, "TRHadGM1": None,
            "THSCPM1": None, "THSCPM2": None, "T1Disp": None }
        from multiverse.expResModifier import readDictFile
        self.log ( f"reading database dict {RED}{dictfile}{RESET}" )
        d = readDictFile ( dictfile )
        self.tempslha = "/dev/shm/temp.slha"
        self.meta = d["meta"]
        self.data = d["data"]
        # self.TLmax = 1. # disregard all results below this
        self.pmax = 0.25 # disregard all results above this
        self.computePDict()
        self.setMassRanges()
        self.ignore_pids = [ 1000003, 2000001,
                             2000002, 2000003, 2000004, 2000011,
                             2000013, 20000015, 1000014, 1000016,
                             2000012, 20000014, 2000016 ]
        force_build  = False
        re = None
        if not force_build:
            re = self.readInitialData()
        if not re:
            self.getTxParamsFromTemplates()
        self.getHighestXSecsFromDatabase( force_build )

    def fmtP ( self, p : float ) -> str:
        return f"{YELLOW}(p={p:.3f}){RESET}"

    def mergeNModels ( self, models : List[Dict], add_timestamp : bool = True )\
            -> Union[None,Dict]:
        """ merge two models, add all particles from both models.
        If a particle appears in both models, its mass will be the average
        of the two, etc.

        :returns: merged model
        """
        import time
        if len(models)== 0:
            self.log ( f"provided empty list of models" )
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

        def computeAverageMassesForLSP ( pid: int, models: List[Dict] ) -> Dict:
            """ for the lsp we really average """
            masses=[]
            for model in models:
                if pid in model["masses"]:
                    masses.append ( model["masses"][pid] )
            return float ( np.mean ( masses ) )

        def computeAverageMassesForPid ( pid : int, models : List[Dict], lspmass = None ) -> Dict:
            """ for the other particles we average over the distance to
            the LSP """
            dm=[]
            if pid == LSP:
                return computeAverageMassesForLSP ( LSP, models )
            lspmasses = []
            for model in models:
                if pid in model["masses"]:
                    if not LSP in model["masses"]:
                        dm.append ( model["masses"][pid] )
                        lspmasses.append ( 0. )
                    else:
                        dm.append ( model["masses"][pid] - model["masses"][LSP] )
                        lspmasses.append ( model["masses"][LSP] )
            avg_delta = float ( np.mean ( dm ) )
            if lspmass == None:
                lspmass = float ( np.mean ( lspmasses ) )
            return lspmass + avg_delta

        def renameParticle ( newpid : int, oldpid : int, model : dict ) -> dict:
            """ rename the particle oldpid to newpid in model 

            :returns: model with oldpid replaced by newpid
            """
            self.log ( f"renaming {oldpid} to {newpid}" )
            if not oldpid in model["masses"]:
                return
            self.log ( f"set m({newpid}) to m({oldpid})={model['masses'][oldpid]}" )
            model["masses"][newpid]=model["masses"][oldpid]
            model["masses"].pop(oldpid)
            if oldpid in model["decays"]:
                model["decays"][newpid] = model["decays"][oldpid]
                model["decays"].pop(oldpid)
            if oldpid in model["decays"]:
                model["decays"][newpid] = model["decays"][oldpid]
                model["decays"].pop(oldpid)
            for mpid, daughters in model["decays"].items():
                newdaughters = {}
                for pdaughters, br in daughters.items():
                    if oldpid in pdaughters:
                        newkeys = tuple( newpid if v == oldpid else v for v in pdaughters )
                        newdaughters[newkeys] = br
                    else:
                        newdaughters[pdaughters]=br
                model["decays"][mpid] = newdaughters
            newssms = {}
            for pids,ssm in model["ssmultipliers"].items():
                if oldpid in pids:
                    newpids = tuple( newpid if v == oldpid else v for v in pids )
                    newssms[newpids]=ssm
                elif -oldpid in pids:
                    newpids = tuple( -newpid if v == -oldpid else v for v in pids )
                    newssms[newpids]=ssm
                else:
                    newssms[pids]=ssm
            model["ssmultipliers"] = newssms
            self.debug ( f'renameParticle returning {model}' )
            return model
                

        def computeAverageDecaysForPid ( pid : int, models : List[Dict] ) -> Dict:
            decays = {}
            Stot = 0.
            nentries = {}
            for model in models:
                if pid in model["decays"]:
                    mdecays = model["decays"][pid]
                    for daughterpids, br in mdecays.items():
                        if not daughterpids in decays:
                            decays[daughterpids] = 0.
                        decays[daughterpids] += br
                        if not daughterpids in nentries:
                            nentries[daughterpids]=0
                        nentries[daughterpids]+=1
                        Stot += br
            for daughterpids, ct in nentries.items():
                if ct>1:
                    decays[daughterpids]/=ct
            #for daughterpids, br in decays.items():
            #    decays[daughterpids]/=Stot
            return decays

        def computeAverageSSMs ( models : List[Dict] ) -> Dict:
            ssms = {}
            for model in models:
                mssms = model["ssmultipliers"]
                for k,v in mssms.items():
                    ssms[k]=v # just bluntly take them over.
            return ssms

        # step #0, rename some particles
        onoffParticles = [ 1000006, 1000023, 1000024 ]
        # for the particles above, we dont merge on with offshell variants
        hasOffshell = set() ## for 1000023, 1000024, 1000006
        for model in models:
            for particle in onoffParticles:
                if particle in model["masses"]:
                    isOn = self.isOnshell ( particle, model["masses"] )
                    if isOn == False:
                        hasOffshell.add ( particle )
        for i,model in enumerate(models):
            for particle in onoffParticles:
                if self.isOnshell ( particle, model["masses"] ) and particle in hasOffshell:
                    newPids = { 1000023: 1000025, 1000024: 1000037,
                                1000006: 2000006 }
                    models[i]=renameParticle ( newPids[particle], particle, model )

        # first, we collect all pids
        pids = collectPids ( models )
        # next we compute average masses
        masses, decays, ssms = {}, {}, {}
        lspmass = computeAverageMassesForLSP ( LSP, models )
        self.log ( f"when merging, new lspmass -> {lspmass:.1f}" )
        masses[LSP] = lspmass
        for pid in pids:
            if pid != LSP:
                mass = computeAverageMassesForPid ( pid, models, lspmass )
                self.log ( f"when merging m({pid}) -> {mass:.1f}" )
                masses[pid]=mass
            piddecays = computeAverageDecaysForPid ( pid, models )
            decays[pid] = piddecays
        ties = [ ( 1000023, 1000024 ) ]
        for pids in ties:
            if pids[0] in masses and pids[1] in masses:
                if self.isOnshell ( pids[0], masses ) == \
                       self.isOnshell ( pids[1], masses ):
                    p = np.random.uniform ( 0, 1 )
                    if p < .1:
                        masses[pids[1]] = masses[pids[0]]
                    if p > .9:
                        masses[pids[0]] = masses[pids[1]]

        ssms = computeAverageSSMs ( models )
        ret["masses"]=masses
        ret["decays"]=decays
        ret["ssmultipliers"]=ssms
        if add_timestamp:
            ret["timestamp"] = time.asctime()
        self.log ( f"we merged {len(models)} models" )
        return ret

    def getHighestXSecsFromDatabase( self, force_build : bool = False ):
        """ we sift through the database, and note the points with the
        highest efficiencies, write into a file 
        :param force_build: if true, then ignore cache
        """
        self.log ( f"effs cache {self.effscachefile}" )
        if not os.path.exists ( self.effscachefile ) and not force_build:
            effscachefile = os.path.abspath ( f"{self.getMyPathName()}/../share/{self.effscachefile}" )
            # self.log ( f"we might have a version at {effscachefile}" )
            if os.path.exists ( effscachefile ):
                #os.symlink ( effscachefile, self.effscachefile )
                lock ( self.effscachefile )
                self.log ( f"copying {self.effscachefile} from {effscachefile}" )
                import shutil
                shutil.copy ( effscachefile, self.effscachefile )
                unlock ( self.effscachefile )

        if os.path.exists ( self.effscachefile ) and not force_build:
            with open ( self.effscachefile, "rt" ) as f:
                txt = f.read()
                self.highestXSecs = eval(txt)
                f.close()
            return

        self.log ( f"now get the highest fiducial xsecs for all results" )
        self.highestXSecs = {}
        from ptools.helpers import computeP
        ers = self.environ.database.getExpResults( dataTypes=["efficiencyMap"] )
        ners = len(ers)
        for i,er in enumerate(ers):
            if True: # i % 10 == 0:
                self.log ( f"computing xsecs for #{i}/{ners}: {er.globalInfo.id}" )
            for ds in er.datasets:
                obsN = ds.dataInfo.observedN
                expBG = ds.dataInfo.expectedBG
                if False: # obsN < expBG: # not interesting
                    continue
                bgErr = ds.dataInfo.bgError
                p = computeP ( obsN, expBG, bgErr )
                if p > 0.3:
                    # not interesting
                    continue
                self.findHighestXSecsFor ( ds )
        lock ( self.effscachefile )
        with open ( self.effscachefile, "wt" ) as f:
            f.write ( f"{self.highestXSecs}\n" )
            f.close()
        unlock ( self.effscachefile )

    def massVecToDict ( self, masses, txname ):
        """ given a masses vector and the txname object,
        translate into a dictionary with pids as keys and
        masses as values
        """
        ret = {}
        txn = txname.txName
        if txn in self.mapTxnames and self.mapTxnames[txn] is not None:
            txn = self.mapTxnames[txn]
        if not txn in self.pidsForTxnames:
            self.error ( f"txname {txn} not in pidsForTxnames" )
            return ret
        pidsD = self.pidsForTxnames[ txn ]
        for idx,pids in pidsD.items():
            if len(pids)==0:
                continue
            # pid = tuple(pids)[0] # FIXME
            for pid in pids:
                ret[pid] = masses[idx]
        return ret

    @lru_cache
    def getXSecDictFor ( self, pids : tuple ):
        sqrts, ewk = 13, "wino"
        massvec = [ 300, 100 ] # doesnt matter for us here
        xsecall,order,comment = self.xsecComputer.getXSecsFor ( pids[0], pids[1],
            sqrts, ewk, massvec )
        if xsecall is None:
            pass
            # self.log ( f"didnt get xsecall for pids={pids}" )
            # return None
        return xsecall,order,comment

    @lru_cache
    def getRefXSecsFor ( self, txname, masses : tuple, pids : tuple ) -> \
            Union[float,None]:
        """ get the reference cross sections for the pt
        :param txname: the txname object
        :param massvec: a list of the masses

        :returns: the cross sections in fb, or None
        """
        if len(pids)==0:
            self.log ( f"pids length?? {pids} {txname} {masses}" )
            return None
            # pids = [ 1000022, 1000022 ]
        pids = list ( pids )
        masses = dict (masses )
        for i,p in enumerate(pids):
            if p == 2000015:
                pids[i]=-1000015
            if p == 2000005:
                pids[i]=-1000005
        if len(pids)==1:
            self.log ( f"pids length?? {pids} {txname} {masses}" )
            pids = [ pids[0], None ]
        if pids == [ 1000021, 1000004, 1000002, 1000001 ]:
            pids = [ 1000001, 1000021 ]
        if pids == [ 1000004, 1000002, 1000001 ]:
            pids = [ -1000001, 1000001 ]
        if len(pids)>2 and 1000022 in pids:
            pids.remove ( 1000022 )
        if pids == [ 1000015, 1000015 ]:
            pids = [ -1000015, 1000015 ]
        if len(pids)==2 and pids[1]<pids[0]:
            pids = [ pids[1], pids[0] ]
        if pids == [ 1000011, 1000013 ]:
            # can use just e+ e-, good enough
            pids = [ -1000011, 1000011 ]
        massvec = []
        for pid in pids: ## thats the order
            if pid in masses:
                massvec.append ( masses[pid] )
            elif abs(pid) in masses:
                massvec.append ( masses[abs(pid)] )
            else:
                self.log ( f"could not find {pid} in {masses}?" )
        if pids == [ 1000002, 1000022 ] or pids == [ 1000001, 1000022 ] or \
                pids == [ 1000004, 1000022 ]:
            # for now replace with C1N2
            self.logThrice ( f"sq-N1 production: {pids}. will replace with C1N2" )
            pids = [ 1000023, 1000024 ]
        # self.log ( f"asking for xsecs for {txname} {massvec} {masses} {pids}" )
        xsecall,order,comment = self.getXSecDictFor ( tuple(pids) )
        if xsecall is None:
            self.logThrice ( f"didnt get xsecall for {txname} pids={pids}" )
            return None
        # self.log ( f"xsecall is {len(xsecall)} order {order} comment {comment}" )
        xsec = self.xsecComputer.interpolate ( massvec, xsecall )
        # self.log ( f"xsecs for {massvec} are {xsec}" )
        if xsec is None:
            xsecmin, xsecmax = min(xsecall.keys()), max(xsecall.keys())
            if massvec[0]< xsecmin:
                xsec = self.xsecComputer.interpolate ( [ xsecmin, xsecmin ], xsecall )
            elif massvec[0]> xsecmax:
                xsec = self.xsecComputer.interpolate ( [ xsecmax, xsecmax], xsecall ) * xsecmax / massvec[0] ## linearly decrease
            else:
                # just randomly try this
                if massvec[1] < 500:
                    massvec[1]+=1e-8
                else:
                    massvec[1]-=1e-8
                xsec = self.xsecComputer.interpolate ( massvec, xsecall )

            if xsec is None:
                self.log ( f"did not get xsec for {txname} {pids} {masses} {massvec}: xsecall were {xsecmin}: {xsecall[xsecmin]} ... {xsecmax}: {xsecall[xsecmax]}" )
        return xsec

    def findHighestXSecsFor ( self, dataset ):
        """ search for highest xsecs in this dataset """
        nmax = 20
        d = TopNDict(nmax=nmax)
        for txname in dataset.txnameList:
            txn = txname.txName
            if txn in self.mapTxnames:
                if self.mapTxnames[txn] is None: # like TRS1, skip 
                    continue
                txn =  self.mapTxnames[txn]
            if not txn in self.pidsForTxnames: # skip this
                self.log ( f"skipping {txn}: not in pidsForTxnames" )
                continue # FIXME we sure?
            pids = list (self.pidsForTxnames[txn][0]) ## the mother pids
            pids.sort ( reverse = True )
            #if len(pids)<2:
            #    self.log ( f"for {txn} we have {pids}" )
            if len(pids)==1:
                if pids[0] in self.xsecComputer.samesignmodes:
                    pids = [ pids[0], pids[0] ]
                else:
                    pids = [ -pids[0], pids[0] ]
            data = txname.txnameData
            d_txn = TopNDict(nmax=int(np.ceil(nmax/len(dataset.txnameList))))
            txn_xsec_tot = 0.
            for pt,eff in zip(data.tri.points,data.y_values):
                # in the mass plane
                massvec = data.inversePCAtransf(pt)
                masses = self.massVecToDict ( massvec, txname )
                refxsec = self.getRefXSecsFor ( 
                        txname, tuple(masses.items()), tuple(pids) )
                if refxsec is None:
                    refxsec = 1e-6 # worst case we fall back to effs,
                    # but multiply with 1e-6 so xsecs win out
                xsec = float ( refxsec * eff )
                while xsec in d: # make sure we dont overwrite
                    xsec += 1e-10
                if type(pt) not in [ list, tuple ]:
                    pt = pt.tolist()
                txn_xsec_tot += xsec
                d_txn[xsec]={ "masses": masses, "txn": txn, "eff": float(eff),
                          "pt": pt, "pids": tuple(pids), "refxsec": float(refxsec),
                          "xsee": float ( refxsec * eff ) }
            for k,v in d_txn.items():
                d[k/txn_xsec_tot]=v
            # d.update ( d_txn ) # didnt rewrite this method for capping
        totxsec = 0. # normalize the keys
        for k,v in d.items():
                totxsec += k
        newd = {}
        for k,v in sorted ( d.items(), reverse = True ):
            newd[k/totxsec] = v
        label = f"{dataset.globalInfo.id}:{dataset.getID()}"
        self.highestXSecs[label]=newd

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

    def decaysAllowedBySLHA ( self, decays : dict ) -> dict:
        """ for txname are decays <ids> allowed for 
        mother <pid>, according to the slha template file? """
        self.log ( f"filter decays allowed by template slha file" )
        ret = {}
        counts = {}
        for mpid, dpids in decays.items():
            if mpid not in [ 1000023, 1000024 ]:
                ret[mpid]=dpids # ugly hack for now
                continue
            temp_tuples = list ( self._pm.decay_tuples[mpid].values() )
            newdpids = set()
            for dpid in dpids:
                adpid = tuple ( map(abs,dpid) )
                ct = temp_tuples.count ( adpid )
                if ct > 0:
                    newdpids.add ( dpid )
                    if not mpid in counts:
                        counts[mpid]={}
                    counts[mpid][adpid]=ct
            ret[mpid]=newdpids
            self.debug ( f"for {mpid} we allow {newdpids}" )
        return ret, counts

    def getTxParamsFor ( self, filename : str ):
        """ get pids, decays for slha template <filename>
        :param filename: e.g. ..../T1.template
        """
        # self.log ( f"getTxParamsFor {filename}" )
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
                ids = decay.ids
                if ids[-1] > 1000000:
                    ids = [ ids[-1] ] + ids[:-1]
                ids = tuple ( ids )
                self.decaysForTxnames[txname][pid].add(ids)
        os.unlink ( self.tempslha )
        self.getDefaultSSMs ( filename )

        if True:
            from base.locker import lock, unlock
            self.fixPidsForTxnames()
            lock ( self.cachefile )
            with open ( self.cachefile, "wt" ) as f:
                f.write ( f"{self.pidsForTxnames}\n" )
                f.write ( f"{self.decaysForTxnames}\n" )
                f.write ( f"{self.ssmsForTxnames}\n" )
                f.close()
            unlock ( self.cachefile )

    def getMyPathName ( self ):
        """ get the pathname of this very file, but symlinks resolved """
        return os.path.dirname ( os.path.realpath(__file__) )

    def readInitialData ( self ) -> bool:
        """ read in all the data (pids,decays,ssms) from the slha files.

        :param force_build: if True, force rebuilding this file
        :returns: False, if no cache file found.
        """
        if not os.path.exists ( self.cachefile ):
            cachefile = os.path.abspath ( f"{self.getMyPathName()}/../share/{self.cachefile}" )
            # self.log ( f"we do not have {self.cachefile}, lets see if we can pull in {cachefile}" )
            if os.path.exists ( cachefile ):
                # os.symlink ( cachefile, self.cachefile )
                lock ( self.cachefile )
                self.log ( f"copying {self.cachefile} from {cachefile}" )
                import shutil
                shutil.copy ( cachefile, self.cachefile )
                unlock ( self.cachefile )
                # self.cachefile = cachefile
            else:
                return False
        self.log ( f"reading in all initial data from {self.cachefile}" )
        with open ( self.cachefile, "rt" ) as f:
            lines = f.readlines()
            f.close()
            self.pidsForTxnames = eval(lines[0])
            self.decaysForTxnames = eval(lines[1])
            self.ssmsForTxnames = eval(lines[2])
        self.fixPidsForTxnames()
        return True

    def fixPidsForTxnames  ( self ):
        # change manually for a few cases
        self.pidsForTxnames["TChiQ"]={0:{1000002,1000022},1:{1000022}}
        self.pidsForTxnames["TScharm"]={0:{1000004,1000022},1:{1000022}}

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
        # self.log ( f"get first file in {tarball}" )
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
        associate_modes = [ ( -1000024, 1000023 ),
                            ( -1000024, 1000022 ),
                            ( -1000024, 1000025 ),
                            ( -1000037, 1000022 ),
                            ( -1000037, 1000023 ) ]
        try:
            xsecs = r.xsections
            ssmpids = set()
            ignore_pids = self.ignore_pids + [ -x for x in self.ignore_pids ] + [ 2212]
            for k,v in xsecs.items():
                # print ( "k=", k, "  v=", v, "type v", type(v) )
                pids = list ( filter(lambda x: x not in ignore_pids, k) )
                pids.sort()
                pids = tuple( pids )
                skip_associate = False
                for am in associate_modes:
                    if pids == am:
                        skip_associate = True
                if skip_associate:
                    continue
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
            txns = set()
            for txn in stats["txns"]:
                if txn in self.mapTxnames and self.mapTxnames[txn] is None:
                    continue
                txns.add ( txn )
            if len(txns)>7 or len(txns)==0:
                ## so many txnames, there isnt much info
                ## this is too vague, or no txns we are using
                continue
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
            if p == 0:
                self.error ( f"for {anaAndSRName} we got p=0" )
                p = 1e-8
            prel = 1. / p
            while prel in prels:
                prel+=1e-6
            value = stats
            value["id"]=anaAndSRName
            prels[ prel ] = value
            ptot += prel
            self.pvalues[anaAndSRName]={ "p": p }
        if len(prels)==0:
            self.error ( "computePDict: no results returned" )
        self.probs = dict ( sorted ( [ (k/ptot,v) for k,v in prels.items() ], 
                            reverse = True ) )
        #probkeys = list ( self.probs.keys() )
        #probkeys.sort (reverse = True )
        #self.probkeys = probkeys

    def randomlyChooseFromDataset ( self, result : dict ) -> dict:
        """ ok, we found a dataset, now we choose a random
        txn and mass point 

        :returns: result dictionary but enriched with txn and mass point
        """

        Id = result["id"]
        idx = list(self.probs.values()).index ( result )
        p = list(self.probs.keys())[idx]
        self.log ( f"of {len(self.probs)} entries we randomly choose #{idx+1}:" )
        if not "new_Z" in result:
            result["new_Z"]=result["orig_Z"]
            result["new_p"]=result["orig_p"]
#            self.pprint ( f"result is {result}" )
        newZ = result["new_Z"]
        self.log ( f"  {GREEN}{Id}{RESET} {self.fmtP(p)} Z={newZ:.2f}" )
        self.log ( f"  `- txns = {', '.join(result['txns'])}" )

        # txns = result["txns"] # .split(",")
        hi_effs = self.highestXSecs[ Id ]
        tot_effs = sum(hi_effs.keys())
        norm_effs = {} # normalized
        #self.log ( f"choosing random txn from {result['id']}: {txn}" )
        for p,pt in hi_effs.items():
            norm_effs[p/tot_effs]=pt
        ctr = 0
        while True:
            choose_pt = np.random.choice(list(norm_effs.values()),
                                         p=list(norm_effs.keys()))
            idx = list(norm_effs.values()).index(choose_pt)
            p = list(norm_effs.keys())[idx]
            txn = choose_pt["txn"]
            result["txn"]=txn
            masses=choose_pt["masses"] # FIXME smear them, and turn into dictionary
            result["masses"]=masses
            sm = ", ".join ( [ f"m({k})={v:.1f}" for k,v in masses.items() ] )
            self.log ( f"for {Id} we randomly pick entry #{idx+1}/{len(norm_effs)}:" )
            self.log ( f"  {GREEN}{txn}: {sm}{RESET} {self.fmtP(p)}" )
            if txn not in self.mapTxnames or self.mapTxnames[txn] is not None:
                break
            ctr += 1
            if ctr > 10:
                self.log ( f"couldnt escape the while-loop at initialiser.py:A" )
                break
            # import sys, IPython; IPython.embed( colors = "neutral" )
        return result

    def randomlyChooseOneResult ( self ) -> Dict:
        """ randomly choose one result from self.probs

        :returns: tuple(txname, result-dictionary)
        result-dictionary is the whole result dictionary of the excess we are
        exploiting

        :returns: result dictionary object
        """
        Id, result = "?", {}
        probs = copy.deepcopy ( self.probs )
        while (Id not in self.highestXSecs) or (Id in self.dsIdsInProposal):
            keys = list(probs.keys())
            key = np.random.choice(keys,p=keys)
            result = probs[key]
            #result = np.random.choice(list(probs.values()),
            #        p=list(probs.keys()) )
            Id = result["id"]
            probs.pop ( key ) ## take it out
            prob_tot = 1. - key
            newprobs = {}
            for k,v in probs.items():
                newprobs[k/prob_tot]=v
            probs = newprobs
        self.dsIdsInProposal.add ( Id )
        result = self.randomlyChooseFromDataset ( result )
        return result
 
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
            if pid in tie:
                otherpid = tie[0] if tie[1] == pid else tie[1]
                if otherpid in pids:
                    return otherpid
        return None

    def isOnshell ( self, pid : int, masses : dict ) -> bool:
        """ determine if a pid like 1000023, 1000024 is onshell

        :returns: true if onshell, false if ofshell, none if
        doesnt apply
        """
        if not pid in masses:
            return None
        if not LSP in masses:
            return None
        massG = 0.
        if pid in [ 1000023, 1000024, 1000006 ]:
            from base.constants import smMasses
            massG = smMasses[pid-1000000]
        else:
            return None
        dm = masses[pid] - masses[LSP]
        if dm > massG:
            return True
        return False

    def getAllowedParticles ( self, constraints ):
        allowed_particles = set()
        test_particles = { "W(": (24,), "Z(": (23,), "e": (11,),
            "mu": (13,), "tau": (15,), "l": (11,13), "L": (11,13,15),
            "W+(": ( 24,), "W-(": ( 24, ), "t(": (6,), "t+(": (6,),
            "t-(": (6,), "b": ( 5, ), "jet": ( 1,2,3,4 ), "q": ( 1,2,3,4 ),
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
        all_channels, counts = self.decaysAllowedBySLHA (  self.decaysForTxnames[txname] )
        allowed_particles = self.getAllowedParticles ( constraints )

        good_channels = {}
        self.debug ( f"for {txname} we filter:" ) 
        for mother,decays in all_channels.items():
            #if mother == 1000024:
            #    self.log ( f"debug allowed {allowed_particles} constraints {constraints}" )
            good_decays = decays
            if len(decays)>1: # no choice for == 1
                good_decays = set()
                for pids in decays:
                    if self.pidsAreGood ( pids, allowed_particles ):
                        good_decays.add ( pids )
            if len(good_decays)>0:
                good_channels[mother]=good_decays
            if good_decays == decays:
                self.debug ( f"  `- {mother}: {decays} (same)" )
            elif good_decays == set():
                self.debug ( f"  `- {mother}: {decays} -> empty" )
            else:
                self.debug ( f"  `- {mother}: {decays} -> {good_decays}" )
        return good_channels, counts

    def getDecaysForTxname ( self, result : Dict ) -> Dict:
        """ get some random decays starting points
        :param result: the whole result dictionary of the excess we are
        """
        txname = result["txn"]
        if not txname in self.decaysForTxnames:
            self.error ( f"we dont have any decays??" )
            sys.exit()

        if not "constraints" in result:
            print ( f"@@22 no constraint:  {result}" )
        constraints = result["constraints"]
        tmp, counts = self.getPossibleChannelsFromConstraints ( txname, constraints )

        decays = { ProtoModel.LSP: {} }
        for mother,daughters in tmp.items():
            nbr_tot = 0.
            ndaughters = len(daughters)
            if not mother in decays:
                decays[mother]={}
            for daughterpids in daughters:
                keys = tuple ( [ abs(dp) for dp in daughterpids ] )
                nbr = random.uniform(0.,1.)
                if not keys in decays[mother]:
                    decays[mother][keys]=nbr
                    nbr_tot += nbr
            for keys, nbr in decays[mother].items():
                decays[mother][keys]= nbr / nbr_tot
            # self.debug ( f"mother {mother} nbr_tot {nbr_tot} decays {decays}" )
        self.log ( f"we randomly choose the decays for {txname}:" )
        # self.log ( f"counts are {counts}" )
        for mother, daughters in decays.items():
            #if mother in counts:
            #    self.log ( f"for {mother} the counts are {counts[mother]}" )
            for dpid, nbr in daughters.items():
                # self.log ( f"dpid is {dpid}" )
                if mother in counts and dpid in counts[mother]:
                    ct = counts[mother][dpid]
                    decays[mother][dpid]/=ct
                    # daughters[dpid]/=ct
            self.log ( f"  `- {mother}: {daughters}" )
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

    def smearMasses ( self, masses : Dict )-> Dict:
        """ smear out the masses in dictionary """
        newmasses = {}
        if not LSP in masses: # in this case we just smear more simply
            for pid,mass in masses.items():
                nmass = mass * scipy.stats.norm.rvs ( loc = 1., scale =.2 )
                newmasses[pid]=nmass
            return newmasses
        oldlspmass = masses[LSP]
        lspmass = float ( masses[LSP]*scipy.stats.norm.rvs ( loc = 1., scale =.2 ) )
        self.log ( f"we randomly smear m(LSP/{LSP}): {masses[LSP]:.1f} -> {lspmass:.1f}" )
        newmasses[LSP]=lspmass
        for pid,mass in masses.items():
            if pid == LSP:
                continue
            tpid = self.tiePids ( pid, masses.keys() )
            if tpid is not None and tpid < pid:
                continue
            wasOnshell = self.isOnshell ( pid, masses )
            nmass = -1. 

            scale = .2
            while nmass < lspmass + 1: ## at least 1 gev distance to lsp
                deltam = max ( 1., mass - oldlspmass )
                nmass = lspmass + float ( deltam * scipy.stats.norm.rvs ( loc = 1., scale = scale) )
                newmasses[pid]=nmass
                isOnsh = self.isOnshell ( pid, newmasses )
                if wasOnshell == True and isOnsh == False:
                    nmass = -1. # continue
                if wasOnshell == False and isOnsh == True:
                    nmass = -1. # continue
                scale *= 1.3
            self.log ( f"we randomly smear m({pid}): {mass:.1f} -> {nmass:.1f}" )
            
            if tpid is not None:
                self.log ( f"we randomly smear m({tpid}): {mass:.1f} -> {nmass:.1f}" )
                newmasses[tpid]=nmass
                
        if 1000023 in newmasses and 1000024 in newmasses:
            p = np.random.uniform ( 0, 1 )
            if p < .2:
                newmasses[1000023] = newmasses[1000024]
            if p > .8:
                newmasses[1000024] = newmasses[1000023]
        return newmasses

    def getRandomSubmodelForTxname ( self, result : Dict ) -> Dict:
        """ given a txname, create a random submodel.
        :param txname: txname for which to create submodel
        :param result: the whole result dictionary of the excess we are
        exploiting
        """
        txname = result["txn"]
        if txname in self.mapTxnames:
            txname = self.mapTxnames[txname]
            self.log ( f"we switch txname from {result['txn']} to {txname}" )
        masses = self.smearMasses ( result["masses"] )
        decays = self.getDecaysForTxname ( result )
        ssms = self.getSSMsForTxname ( txname )
        model = { "masses": masses, "decays": decays, "ssmultipliers": ssms }
        return model

    def propose ( self ):
        """ propose a random initial model. """
        # choose a random txn
        # self.getRandomMassForLSP()
        self.dsIdsInProposal = set()
        submodels = []
        nmodels = np.random.choice ( [1,2,3] )
        self.log ( f"proposed model {self.dictfile} will consist of {RED}{nmodels} submodels{RESET}." )
        for i in range(nmodels):
            result = self.randomlyChooseOneResult()
            submodel = self.getRandomSubmodelForTxname ( result )
            if submodel == None:
                self.warn ( f"got none as random submodel" )
                continue
            submodels.append ( submodel )
        self.submodels = submodels
        model = self.mergeNModels ( submodels )
        
        if self.environ.allowN1N1Prod and model is not None:
            ## the factor of .2 below is because our reference xsecs for N1N1 are those of C1C1
            ## so too high
            ssm = float ( .2 * np.exp ( scipy.stats.norm.rvs() ) )
            model["ssmultipliers"][(1000022,1000022)]=ssm
        if True:
            from ptools.helpers import py_dumps
            ds = py_dumps ( model )
            self.log ( f"{RED}we propose the model:{RESET}" )
            self.log ( f"\n{ds}" )
        return model

    def predictForModel ( self, model : Union[dict,None] = None ) -> dict:
        """ call pr.predict for the given model """
        if model == None:
            model = self.propose()
            self.debug ( f"predictForModel {model}" )
        ma = Manipulator( model, walkerid = self.walkerid, 
                          environ = self.environ )

        self.pr = Predictor( self.walkerid, self.environ )
        self.log( f"predictForModel, run predict(1)" )
        self.pr.predict ( ma, keep_predictions = True, force_computation_K=True )
        # FIXME why do I have to call this twice?
        self.cr = Critic( self.walkerid, self.environ )
        crr = False, "critic not run"
        ctr = 0
        scale = 1.
        self.log( f"predictForModel, K={ma.M.K}, run critics" )
        while ctr < 3:
            crr = self.cr.predict_critic ( ma.M )
            self.log ( f"predictForModel asking the critic ({ctr}): {crr}" )
            verdict = crr[0]
            if verdict == True:
                break
            scale *= .6
            ma.M.rescaleXSecsBy ( scale ) ## FIXME do sth smarter here!
            ctr += 1
        self.log( f"predictForModel, run predict(2)" )
        self.pr.predict ( ma, keep_predictions = True, force_computation_K=True )

        # update the ssms, they might have changed
        model["ssmultipliers"]=ma.M.ssmultipliers
        ret = { "ma": ma, "pr": self.pr, "model": model, "cr": self.cr }
        ret["K"] = ma.M.K
        ret["crr"]=crr
        ret["TL"] = ma.M.TL
        return ret

    def bestOfN ( self, n : int = 5 ):
        """ get <n> initial contender models, return the best """
        models = {}
        for i in range(n):
            self.log ( f"{CYAN}best of N: {i+1}/{n}{RESET}: starting" )
            ret = self.predictForModel ( None )
            model = ret["model"]
            K = ret["K"]
            crr = ret["crr"]
            model["K"] = K
            model["TL"] = ret["TL"]
            model["crr"] = ret["crr"]
            models[ ret["K"] ] = model
            self.log ( f"{CYAN}best of N: {i+1}/{n}{RESET}: got K={formatObject(K,'.1f')} cr={crr}{RESET}" )
        # print ( "bestOfFive: {d}" )
        keys = [ k for k in models.keys() if k is not None ]
        skeys = [ f"{k:.1f}" for k in keys ]
        if len(keys) == 0:
            self.log ( f"best of {n} got us no good model" )
            return None
        maxK = max( keys )
        self.log ( f"best of {n} got: K=max({', '.join(skeys)})={maxK:.2f}" )
        return models[maxK]

    def interact ( self ):
        """ interactive shell, for debugging and development """
        print ( f"try e.g:" )
        print ( f'ret = self.predictForModel ( )' )
        print ( f'globals().update(ret)' )
        print ( f"print ( K, model )" )
        import IPython
        IPython.embed( colors = "neutral" )

if __name__ == "__main__":
    import argparse
    argparser = argparse.ArgumentParser(
            description='CLI of initialiser' )
    argparser.add_argument ( '-d', '--dictfile',
            help='input database dict file ["signal_database.dict"]',
            type=str, default="*_database.dict" )
    argparser.add_argument ( '-v', '--verbose',
            help='verbose', action="store_true" )
    argparser.add_argument ( '-r', '--recompute_cache',
            help='dont use pids.cache and xsecs.cache cache files', action="store_true" )
    ## 310.dict is also a good default, for the actual observations
    args = argparser.parse_args()
    if args.recompute_cache:
        if os.path.exists ( Initialiser.cachefile ):
            print ( f"[Initialiser] deleting old {Initialiser.cachefile}" )
            os.unlink ( Initialiser.cachefile )
        if os.path.exists ( Initialiser.effscachefile ):
            print ( f"[Initialiser] deleting old {Initialiser.effscachefile}" )
            os.unlink ( Initialiser.effscachefile )
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
    environ = RunEnviron( "run.dict" )
    ini = Initialiser( "ini", args.dictfile, environ, verbose = args.verbose )
    ini.interact()
