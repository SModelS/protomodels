#!/usr/bin/env python3

""" Class that encapsulates a BSM model. """

__all__ = [ "ProtoModel" ]

import random, tempfile, os, time, colorama, copy, sys, pickle, random
sys.path.insert(0,"../")
sys.path.insert(0,f"{os.environ['HOME']}/git/smodels/")
from smodels.tools.wrapperBase import WrapperBase
# the default tempdir of wrapper base is /tmp
# WrapperBase.defaulttempdir="./" ## keep the temps in our folder
# WrapperBase.defaulttempdir="/dev/shm" ## keep the temps in shared memory
from ptools.refxsecComputer import RefXSecComputer
from smodels.base.physicsUnits import TeV, fb
from ptools import helpers
from ptools.sparticleNames import SParticleNames
from typing import Union, List, Tuple

from smodels.base.smodelsLogging import setLogLevel
from builder.loggerbase import LoggerBase
setLogLevel ( "error" )


class ProtoModel ( LoggerBase ):
    """ encodes one theoretical model, i.e. the particles, their masses, their
        branchings, their signal strength modifiers.
    """
    # __slots__ = [ "walkerid", "keep_meta" ]

    LSP = 1000022 ## the LSP is hard coded
    SLHATEMPDIR = "/tmp/" # "./" where do i keep the temporary SLHA files?
    #SLHATEMPDIR = "/dev/shm/" # "./" where do i keep the temporary SLHA files?

    def __init__ ( self, walkerid : int = 0, keep_meta : bool = True, dbversion : str = "????" ):
        """
        :param keep_meta: If True, keep also all the data in best combo (makes
                          this a heavyweight object)
        :param walkerid: id of current walker
        :param dbversion: the version of the database, to track provenance
        """
        super(ProtoModel,self).__init__ ( walkerid )
        self.walkerid = walkerid
        self.keep_meta = keep_meta ## keep all meta info? big!
        self.version = 1 ## version of this class
        self.maxMass = 2400. ## maximum masses we consider
        self.step = 0 ## count the steps
        self.dbversion = dbversion ## keep track of the database version
        self.particles = [ 1000001, 2000001, 1000002, 2000002, 1000003, 2000003,
                  1000004, 2000004, 1000005, 2000005, 1000006, 2000006, 1000011,
                  2000011, 1000012, 1000013, 2000013, 1000014, 1000015, 2000015,
                  1000016, 1000021, 1000022, 1000023, 1000025, 1000035, 1000024,
                  1000037 ]
        self.onesquark = False ## only one light squark
        self.twosquark = False  ## a few squarks, but not all
        self.manysquark = True ## many squarks
        if self.onesquark:
            self.particles = [ 1000001, 1000005, 1000006, 1000011, 1000012,
                      1000013, 1000014, 1000015,  1000016, 1000021, 1000022,
                      1000023, 1000025, 1000024, 1000037 ]
            self.templateSLHA = "templates/template_1q.slha"
        if self.twosquark:
            self.particles = [ 1000001, 1000002, 1000004, 1000005, 1000006, 1000011,
                      1000012, 1000013, 1000014, 1000015, 1000016, 1000021, 1000022,
                      1000023, 1000025, 1000024, 1000037 ]
            self.templateSLHA = "templates/template_2q.slha"
        if self.manysquark:
            self.particles = [ 1000001, 1000002, 1000003, 1000004, 1000005, 1000006,
                      2000005, 2000006, 1000011, 1000012, 1000013, 1000014, 1000015,
                      1000016, 1000021, 1000022, 1000023, 1000025, 1000024, 1000037 ]
            self.templateSLHA = "templates/template1g.slha"
            if False:
                self.particles.append ( 2000021 )
                self.particles.append ( 3000006 )
                self.templateSLHA = "templates/template2g.slha"
            # self.templateSLHA = "templates/template_many.slha"
        self.templateSLHA = os.path.join ( os.path.dirname ( __file__ ), self.templateSLHA )
        self.computer = RefXSecComputer()
        self.codeversion = "2.0"
        self.initializeModel()

    def initializeModel(self):
        """Use the template SLHA file to store possible decays and initialize the LSP"""

        #Make sure the masses, decays and multipliers are empty
        self.tpList = [] ## store information about the theory predictions
        self.rvalues = [] ## store the r values of the exclusion attempt
        self.llhd=0.
        self.muhat = 1.
        self.mumax = None # the maximum mu allowed by the critic
        self.Z = 0.0
        self.K = None
        self.letters = ""
        self.description = ""
        self.bestCombo = None
        self.decays = {} ## the actual branchings
        self.masses = {}
        self.possibledecays = {} ## list all possible decay channels
        self.decay_keys = {} #list the key associated with each decay of a pid
        self._stored_xsecs = () #Store cross-sections. It should only be accesses through getXsecs()!
        self._xsecMasses = {} #Store the masses used for computing the cross-sections
        self._xsecSSMs = {} #Store the signal strenght multiplier used for computing the cross-sections
        self.ssmultipliers = {} ## signal strength multipliers
        ## Inititiaze LSP
        self.masses[ProtoModel.LSP]=random.uniform(200,500)
        self.decays[ProtoModel.LSP]= {}
        pids = [(self.LSP,self.LSP)]
        if self.hasAntiParticle(self.LSP):
            pids += [(self.LSP,-self.LSP),(-self.LSP,-self.LSP)]
        for pidpair in pids:
            self.ssmultipliers[tuple(sorted(pidpair))]= 1.0
        
        slha_decay_keys = []
        
        with open ( self.templateSLHA ) as slhaf:
            tmp = slhaf.readlines()
            for line in tmp:
                p = line.find("#" )
                if p > -1:
                    line = line[:p]
                if "D" in line and not "DECAY" in line:
                    slhaline = line.strip().split(" ")
                    slhaline = [l for l in slhaline if l!='']
                    #decay_key = [slhaline[0],slhaline[1:]]
                    slha_decay_keys.append(slhaline)
                            
        
        for p in self.particles:
            decays = []
            dkey = {}
            for key in slha_decay_keys:
                if f"D{p}" in key[0]:
                    dpid,dpid2,dpid3,dpd = None,None,None,None

                    if int(key[1]) == 2:                    #2body decay
                        dpid = abs(int(key[2]))
                        dpid2 = abs(int(key[-1]))
                        dpd = (dpid,dpid2)
                    elif int(key[1]) == 3:                  #3body decay
                        dpid = abs(int(key[2]))
                        dpid2 = abs(int(key[3]))
                        dpid3 = abs(int(key[4]))
                        dpd = (dpid,dpid2,dpid3)

                    decays.append ( dpd )
                    dkey.update({dpd: key[0]})
            
            self.possibledecays[p]=decays
            self.decay_keys[p] = dkey

    def __str__(self):
        """ return basic information on model
        """
        namer = SParticleNames ( susy=False )

        pNames = [namer.asciiName ( pid ) for pid in self.unFrozenParticles()]
        pNames = ','.join(pNames)
        pStr = f'ProtoModel ({pNames}):'
        if self.K:
            pStr += f' K = {self.K:1.2f}'
        else:
            pStr += f' K = {self.K}'
        if self.Z:
            pStr += f', Z = {self.Z:1.2f}'
        else:
            pStr += f', Z = {self.Z}'

        return pStr

    def __repr__(self):
        """ shortened version of __str__"""
        sK = str(self.K)
        import numpy as np
        if type(self.K) in [ float, np.float64 ]:
            sK=f"{self.K:1.2f}"
        sZ = str(self.Z)
        if type(self.Z) in [ float, np.float64 ]:
            sZ=f"{self.Z:1.2f}"
        pStr = f'ProtoModel ({sK}, {sZ})'
        return pStr

    def hasAntiParticle ( self, pid ):
        """ for a given pid, do i also have to consider its antiparticle
            -pid in the signal strength multipliers? """
        if pid in [ 1000021, 1000022, 1000023, 1000025, 1000035, 1000012,
                    1000014, 1000016, 2000012, 2000014, 2000016, 2000021 ]:
            return False
        return True

    def toTuple ( self, pid1 : int, pid2 : int ):
        """ turn pid1, pid2 into a sorted tuple """
        a=[pid1,pid2]
        a.sort()
        return tuple(a)

    def getXsecs(self) -> Tuple[List,str]:
        """
        Return the cross-sections.
        If they have already been computed (and stored in self._stored_xsecs)
        AND the masses and signal strength multipliers habe not been modified, return the stored value.
        Otherwise, re-compute the cross-sections.

        :return: list of cross-sections
        """

        #If xsecMasses has not been defined or differs from current masses,
        #recompute xsecs
        if self.masses == self._xsecMasses and self.ssmultipliers == self._xsecSSMs:
            if self._stored_xsecs and len(self._stored_xsecs)>0:
                return self._stored_xsecs

        self.delXSecs()
        #If something has changed, re-compute the cross-sections.
        #Xsecs are computed, self._xsecMasses and self._xsecSSM are updated.
        #The results are sored in the SLHA and self._stored_xsec.
        self.computeXSecs()

        return self._stored_xsecs

    def getOpenChannels(self,pid : int ):
        """get the list of open decay channels for particle pid. Open channels are
        the decays to unfrozen particles and to lighter particles.

        :param pid: PID for particle

        :return: List with the daughter pids for each decay channel
        """

        #Get list of possible decay channels:
        openChannels = set()
        unfrozen = self.unFrozenParticles()
        smMasses = {  6: 173., 24: 81., 23: 92., 25: 125.,
                         15: 1.77, 4: 1.2, 5: 5.0 }
        #Get all relevant masses
        allMasses = dict([[pid,mass] for pid,mass in self.masses.items()])
        allMasses.update(smMasses)

        for dpid in self.possibledecays[pid]:
            #Get the list of BSM particles in the decay:
            if isinstance(dpid,(list,tuple)):
                pidList = [abs(p) for p in dpid if abs(p) in self.particles]
            else:
                self.highlight ( "warn", f"a decay channel without the SM particle is specified in {pid}:{str(dpid)}" )
                pidList = [abs(dpid)]
            #Skip decays to unfrozen particles
            if not all([dp in unfrozen for dp in pidList]):
                continue
            #Get total daughter mass (it should only be a single mass)
            mdaughter = sum([allMasses[abs(p)] for p in dpid if abs(p) in allMasses])

            #Skip decays to heavier particles
            if mdaughter >= self.masses[pid]:
                continue

            openChannels.add ( dpid )

        openChannels = list(openChannels)
        
        return openChannels

    def frozenParticles ( self ):
        """ returns a list of all particles that can be regarded as frozen, i.e.
        are not in the unfrozen list."""

        unfrozen = self.unFrozenParticles()
        ret = [pid for pid in self.particles if not pid in unfrozen]
        return ret

    def cleanBestCombo ( self ):
        """ remove unneeded stuff before storing """
        if hasattr ( self, "keep_meta" ) and self.keep_meta:
            return ## dont remove best combo
        from tester.combiner import Combiner
        combiner = Combiner( self.walkerid )
        if hasattr ( self, "bestCombo" ) and self.bestCombo != None:
            self.bestCombo = combiner.removeDataFromBestCombo ( self.bestCombo )

    def almostSameAs ( self, other ):
        """ check if a model is essentially the same as <other> """

        if self.masses.keys() != other.masses.keys():
            return False

        massDiff = [abs(m-other.masses[pid])/m for pid,m in self.masses.items() if m]
        if max(massDiff) > 1e-5:
            return False

        ## now check ssmultipliers
        pidpairs = set ( self.ssmultipliers.keys() )
        pidpairs = pidpairs.union ( set ( other.ssmultipliers.keys() ) )
        for pidpair in pidpairs:
            ss = 1.
            if pidpair in self.ssmultipliers.keys():
                ss = self.ssmultipliers[pidpair]
            os = 1.
            if pidpair in other.ssmultipliers.keys():
                os = other.ssmultipliers[pidpair]
            if ss == 0.:
                if os == 0.:
                    continue
                else:
                    return False
            if abs ( ss - os ) / ss > 1e-6:
                return False
        ## now check decays
        pids = set ( self.decays.keys() )
        pids = pids.union ( set ( other.decays.keys() ) )
        for pid in pids:
            sdecays, odecays = {}, {}
            if pid in self.decays:
                sdecays = self.decays[pid]
            if pid in other.decays:
                odecays = other.decays[pid]
            dpids = set ( sdecays.keys() )
            dpid = dpids.union ( set ( odecays.keys() ) )
            for dpid in dpids:
                sbr, obr = 0., 0.
                if dpid in sdecays:
                    sbr = sdecays[dpid]
                if dpid in odecays:
                    obr = odecays[dpid]
                if sbr == 0.:
                    if obr < 1e-6:
                        continue
                    else:
                        return False
                if abs ( sbr - obr ) / sbr > 1e-6:
                    return False
        return True

    def unFrozenParticles ( self, withLSP : bool = True ):
        """ returns a list of all particles in self.masses with
            mass less than 100 TeV """

        ret = []
        for m,v in self.masses.items():
            if abs(v)<1e5:
                ret.append(m)
        if not withLSP and self.LSP in ret:
            ret.remove(self.LSP)
        return ret

    def printMasses( self ):
        """ convenience function to print masses with particle names """
        particles = []
        namer = SParticleNames ( susy=False )
        for pid,m in self.masses.items():
            if m > 99000:
                continue
            particles.append ( f"{namer.asciiName ( pid )}: {m}" )
        print ( ", ".join ( particles ) )

    def computeXSecs ( self, keep_slha : bool = False ):
        """ compute xsecs given the masses and signal strength multipliers of the model.
         The results are stored in self._stored_xsecs and should be accessed through getXsecs.
        :param keep_slha: if true, then keep slha file at the end

        """

        hasComputed = False
        countAttempts = 0
        while not hasComputed:
            tmpSLHA = ""
            try:
                xsecs = []
                #Create temporary file with the current model (without cross-sections)
                tmpSLHA = tempfile.mktemp( prefix=f".{self.walkerid}_xsecfile",
                                           suffix=".slha",dir=self.SLHATEMPDIR )
                tmpSLHA = self.createSLHAFile(tmpSLHA, addXsecs = False)
                for sqrts in [8, 13]:
                    self.computer.compute( sqrts, tmpSLHA, ssmultipliers = self.ssmultipliers )
                    # for x in self.computer.loXsecs:
                    #     xsecs.append ( x )
                    # self.computer.loXsecs = []
                    for x in self.computer.xsecs:
                        xsecs.append ( x )
                    self.computer.xsecs = []
                comment = f"produced at step {self.step}"
                pidsp = self.unFrozenParticles()
                pidsp.sort()
                namer = SParticleNames ( susy=False )
                prtcles = ", ".join ( map ( namer.asciiName, pidsp ) )
                self.log ( f"done computing {len(xsecs)} xsecs for pids {prtcles}" )
                self._stored_xsecs = ( xsecs, comment )
                self._xsecMasses = dict([[pid,m] for pid,m in self.masses.items()])
                self._xsecSSMs = dict([[pid,ssm] for pid,ssm in self.ssmultipliers.items()])
                hasComputed = True
                if not keep_slha and os.path.exists ( tmpSLHA ): ## remove
                    os.remove( tmpSLHA )
                break
                #Remove temp file
            except Exception as e:
                if not keep_slha and os.path.exists ( tmpSLHA ): ## remove
                    os.remove( tmpSLHA )
                countAttempts += 1
                if countAttempts > 1:
                    self.log( f"error computing cross-sections: {e}, attempt # {countAttempts}" )
                    self.pprint( f"error computing cross-sections: {e}, attempt # {countAttempts}" )
                    import traceback
                    trackback.print_stack()
                # helpers.cpPythia8()
                time.sleep ( random.uniform ( 5, 10 ) )
                if countAttempts > 5:
                    break

        if keep_slha:
            self.createSLHAFile( self.currentSLHA, addXsecs = True )

    def rescaleXSecsBy(self, s : float ):
        """rescale the stored cross-sections by a factor s"""

        #Before rescaling, make sure we get the latest cross-sections:
        x = self.getXsecs()
        xsecs = x[0]
        comment = x[1]
        for xsec in xsecs:
            xsec.value *= s
        for k,v in self.ssmultipliers.items():
            self.ssmultipliers[k] = v * s

        self._stored_xsecs = (xsecs,comment)
        self._xsecSSMs = dict([[pid,ssm] for pid,ssm in self.ssmultipliers.items()])

    def delCurrentSLHA ( self ):
        """ remove current slha file, if it exists """
        if hasattr ( self, "currentSLHA" ) and type(self.currentSLHA)==str and \
                os.path.exists ( self.currentSLHA ):
            # print ( "[protomodel] del", self.currentSLHA )
            os.unlink ( self.currentSLHA )

    def createNewSLHAFileName ( self, prefix : str = "cur" ) -> str:
        """ create a new SLHA file name. Needed when e.g. unpickling
        :returns: slha filename
        """
        self.delCurrentSLHA()
        self.currentSLHA = tempfile.mktemp( prefix=f".{prefix}{self.walkerid}_",
                    suffix=".slha",dir=self.SLHATEMPDIR)
        return self.currentSLHA

    def checkTemplateSLHA ( self ):
        if not os.path.exists ( self.templateSLHA ):
            if "/mnt/hephy/" in self.templateSLHA:
                trySLHA = self.templateSLHA.replace(f"{os.environ['CODEDIR']}/smodels-utils/protomodels/","./" )
                if os.path.exists ( trySLHA ):
                    self.templateSLHA = trySLHA
                    return

    def writeSLHAFile ( self, outputSLHA : str ):
        """ write the slha file, plug in protomodel params """
        #Get template data:
        with open( self.templateSLHA ) as f:
            lines=f.readlines()
        unfrozen = self.unFrozenParticles()

        with open(outputSLHA,'wt') as outF:
            for i,l in enumerate(lines):
                for pid in self.particles:
                    #Skip lines which have no mass or decay tags
                    if not f"M{pid}" in l and not f"D{pid}" in l:
                        continue

                    #Get information for particle
                    if pid in unfrozen:
                        mass = self.masses[pid]
                        decays = self.decays[pid]
                    else:
                        mass = 1e6 #decoupled mass
                        decays = {} #no decays for frozen particles

                    #Replace mass tag:
                    if f"M{pid}" in l:
                        l = l.replace( f"M{pid}", f"{mass:.1f}" )
                    else:
                        decayTag = l.strip().split()[0]
                        decayPids = decayTag.replace('D','').split('_')
                        dpids = tuple([int(p) for p in decayPids[1:]]) #daughter pids
                        if len(dpids) == 1:
                            dpids = dpids[0]
                        if dpids in decays:
                            br = decays[dpids]
                            l = l.replace(decayTag, f"{br:.5f}" )
                        else:
                            l = ""

                #Only write line if it is not empty
                if l:
                    outF.write(l)
            outF.close()

    def createSLHAFile ( self, outputSLHA : Union[str,None] = None,
                         addXsecs : bool = True ):
        """ Creates the SLHA file with the masses, decays and cross-sections stored in the model.

        :param outputSLHA: Name of the SLHA file to be created. If None a tempfile will be created and
                           its name will be stored in self.currentSLHA.
        :param addXsecs: If True, include cross-sections in the file, else only write spectrum and decays.

        :return: Name of the SLHA file created
        """
        self.delCurrentSLHA()

        #If output is not defined, create file and store in self.currentSLHA
        if outputSLHA is None:
            self.createNewSLHAFileName()
            outputSLHA = self.currentSLHA

        #Set template file (if not yet defined)
        self.checkTemplateSLHA()

        #Replace masses and decays with values for the unFrozenParticles:
        self.writeSLHAFile ( outputSLHA )

        ctAttempts = 0
        hasXSecs = False
        #Add cross-sections:
        if addXsecs:
            while not hasXSecs:
                # Cross-sections will be computed if something has changed
                xsecs = self.getXsecs()
                #print ( "[protomodels] adding xsecs to", outputSLHA )
                #for xsec in xsecs[0]:
                #    print ( "[protomodel] adding xsec", str(xsec) )
                if len(xsecs)>0:
                    if not os.path.exists ( outputSLHA ):
                        self.writeSLHAFile ( outputSLHA )
                    self.computer.addXSecToFile( xsecs[0], outputSLHA )
                    self.computer.addMultipliersToFile ( self.ssmultipliers, outputSLHA )
                    self.computer.addCommentToFile ( xsecs[1], outputSLHA )
                    hasXSecs = True
                else:
                    ctAttempts += 1
                    self.pprint ( f"empty cross section container at attempt {ctAttempts}? whats going on?"  )
                    self.delXSecs()
                    if ctAttempts > 5:
                        break
                    time.sleep ( random.uniform ( 0.5, 2.*ctAttempts ) )

        return outputSLHA

    def dict ( self ):
        """ return the dictionary that can be written out """
        tmp = self.getXsecs()[0]
        xsecs = {}
        for xsec in tmp:
            xsecs[(xsec.pid,xsec.info.sqrts.asNumber(TeV))]=xsec.value.asNumber(fb)
        return { "masses": self.masses, "ssmultipliers": self.ssmultipliers,
                 "decays": self.decays, "xsecs[fb]": xsecs }

    def relevantSSMultipliers ( self ):
        """ of all the ss mulipliers, return only the relevant ones,
            i.e. the ones for unfrozen particles and value != 1 """
        ret = {}
        frozen = self.frozenParticles()
        for pids,v in self.ssmultipliers.items():
            if abs ( v - 1. ) < 1e-5:
                continue
            isRelevant = True
            for pid in pids:
                if abs(pid) in frozen:
                    isRelevant = False
            if isRelevant:
                ret[pids]=v
        return ret

    def describe ( self ):
        """ describe a bit the protomodel """
        ndecays,nd = 0, 0
        for k,v in self.decays.items():
            if k == ProtoModel.LSP: ## dont count LSP
                continue
            ndecays += len(v)
            nd += 1
        nssms = len(self.ssmultipliers)
        print ( f"{len(self.masses)} masses, {ndecays}[{nd}] decays, {nssms} ss multipliers" )

    def delXSecs ( self ):
        """ delete stored cross section, if they exist """
        self._stored_xsecs = ()
        self._xsecMasses = {}
        self._xsecSSMs = {}

    def copy(self, cp_predictions : bool = False):
        """
        Create a copy of self. If cp_predictions the bestCombo and tpList attributes
        is copied using deepcopy.

        :returns: copy of protomodel
        """

        #Initialize empty model:
        newmodel = self.__class__( self.walkerid )

        #Copy information
        newmodel.keep_meta = self.keep_meta
        newmodel.maxMass = self.maxMass
        newmodel.step = self.step
        newmodel.dbversion = self.dbversion
        newmodel.codeversion = self.codeversion
        newmodel.particles = self.particles[:]
        newmodel.onesquark = self.onesquark ## only one light squark
        newmodel.twosquark = self.twosquark  ## a few squarks, but not all
        newmodel.manysquark = self.manysquark ## many squarks
        newmodel.templateSLHA = self.templateSLHA[:]
        newmodel.possibledecays = dict([[key,val] for key,val in self.possibledecays.items()])
        decayDict = {}
        for pid,dec in self.decays.items():
            decayDict[pid] = dict([[dpids,br] for dpids,br in dec.items()])
        newmodel.decays = decayDict
        newmodel.masses = dict([[pid,mass] for pid,mass in self.masses.items()])
        newmodel.ssmultipliers = dict([[pidPair,mass] for pidPair,mass in self.ssmultipliers.items()])
        newmodel.rvalues = self.rvalues[:]
        newmodel.llhd = self.llhd
        newmodel.muhat = self.muhat
        newmodel.mumax = self.mumax
        newmodel.Z = self.Z
        newmodel.K = self.K
        newmodel.letters = self.letters[:]
        newmodel.description = self.description[:]
        newmodel._stored_xsecs = copy.deepcopy(self._stored_xsecs)
        newmodel._xsecSSMs = dict([[pid,ssm] for pid,ssm in self._xsecSSMs.items()])
        newmodel._xsecMasses = dict([[pid,m] for pid,m in self._xsecMasses.items()])
        if cp_predictions:
            newmodel.tpList = copy.deepcopy(self.tpList)
            newmodel.bestCombo = copy.deepcopy(self.bestCombo)

        return newmodel

    def lightCopy(self,rmAttr=None):
        """Makes a light copy of the model using helpers.lightObjCopy.
        If rmAttr is None, it will remove the default attributes defined in
        helpers.lightObjCopy."""

        if rmAttr is not None:
            return helpers.lightObjCopy(self,rmAttr=rmAttr)
        else:
            return helpers.lightObjCopy(self)

if __name__ == "__main__":
    p = ProtoModel( 1 )
    p.createSLHAFile()
    p.computeXSecs()
