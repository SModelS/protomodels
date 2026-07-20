#!/usr/bin/env python3

"""Data class representing a BSM physics model.

Encodes one theoretical model: particles, their masses, branchings,
and signal strength modifiers.  Provides SLHA file creation and
cross-section computation.
"""

__all__ = ["ProtoModel"]

import copy
import os
import pickle
import sys
import tempfile
import time
from typing import Dict, List, Optional, Set, Tuple, Union

import numpy as np

sys.path.insert(0, "../")

from smodels.tools.wrapperBase import WrapperBase
from smodels.base.physicsUnits import TeV, fb
from smodels.base.smodelsLogging import setLogLevel

from base.loggerbase import LoggerBase
from base.runEnviron import RunEnviron
from ptools.refxsecComputer import RefXSecComputer
from ptools import helpers
from ptools.helpers import formatObject
from ptools.sparticleNames import SParticleNames

setLogLevel("error")


class ProtoModel(LoggerBase):
    """Encodes one BSM model: particle content, masses, branchings, and
    signal-strength modifiers.

    The class is deliberately kept data-centric; algorithms that modify
    models live in :class:`builder.manipulator.Manipulator`.
    """

    LSP: int = 1000022  # the LSP is hard coded
    SLHATEMPDIR: str = "/tmp/"

    def __init__(
        self,
        walkerid: Union[str, int] = 0,
        keep_meta: bool = True,
        environ: Optional[RunEnviron] = None,
    ) -> None:
        """
        :param keep_meta: If ``True``, keep all data in best combo (heavyweight).
        :param walkerid: ID of current walker.
        :param environ: A :class:`RunEnviron` instance (must not be ``None``).
        """
        assert environ is not None, "set environ!"
        super().__init__(walkerid)
        self.walkerid = walkerid
        self.keep_meta = keep_meta
        self.version: int = 1
        self.maxMass: float = 2400.0
        self.environ = environ
        self.step: int = 0
        # Cache invalidation markers for particle lists
        self._unfrozen_cache: Optional[List[int]] = None
        self._frozen_cache: Optional[List[int]] = None
        self._masses_snapshot: Optional[Dict] = None
        self.getParticleContent()
        self.computer = RefXSecComputer(
            verbose=False,
            allowN1N1Prod=environ.allowN1N1Prod,
            walkerid=walkerid,
        )
        self.protomodels_version = "2.0"
        self.initializeModel()

    @property ## convenience
    def allowN1N1Prod(self):
        return self.computer.allowN1N1Prod

    @allowN1N1Prod.setter ## convenience
    def allowN1N1Prod(self,flag : bool ):
        self.computer.allowN1N1Prod = flag

    def getParticleContent ( self ):
        """ for self.environ.templateSLHA, get its particle content as a list.
        save the content in self.particles.
        also, define potential forced_degeneracies
        """
        assert os.path.exists ( self.environ.templateSLHA ), \
                f"{self.environ.templateSLHA} does not exist"
        particles = set()
        mass_params = set()
        slha = ""
        ## this is a list of force degeneracies
        self.forced_degeneracies = []
        self.decaylessParticles = ( ProtoModel.LSP, )
        # decaylessParticles = [ ProtoModel.LSP, 1000023, 1000024 ]
        with open ( self.environ.templateSLHA, "rt" ) as f:
            lines = f.readlines()
            for line in lines:
                # stop at the decays
                if line.startswith ( "DECAY" ):
                    break
                slha += line
            for line in lines:
                if "DEGENERACY:" in line:
                    p1 = line.find("DEGENERACY:" )
                    token = line[p1+11:]
                    p2 = token.find("#")
                    if p2 > -1:
                        token = token[:p2]
                    token = token.strip()
                    try:
                        parsed = eval(token)
                        self.forced_degeneracies.append(parsed)
                    except (SyntaxError,ValueError) as e:
                        self.error ( f"cannot parse {token}: {e}" )
                if "DECAYLESSPARTICLES:" in line:
                    p1 = line.find("DECAYLESSPARTICLES:" )
                    token = line[p1+19:]
                    p2 = token.find("#")
                    if p2 > -1:
                        token = token[:p2]
                    token = token.strip()
                    try:
                        parsed = eval(token)
                        self.decaylessParticles=parsed
                    except (SyntaxError,ValueError) as e:
                        self.error ( f"cannot parse {token}: {e}" )

        # print ( "slha", slha )
        import pyslha
        f = pyslha.readSLHA ( slha )
        masses = f.blocks["MASS"]
        for pid,mass in masses.items():
            if type(mass) in [ str ] and mass.startswith("M"):
                mass_param = int ( mass.replace("M","") )
                particles.add ( pid )
                assert mass_param == pid, f"we assume that the mass parameter {mass} has the same number as the particle {pid}"
        self.particles = list ( particles ) # thats the particles

    def initializeModel(self) -> None:
        """Initialise model attributes from the template SLHA file.

        Sets up decays, mass parameters, and the LSP mass.
        """
        self._invalidateParticleCache()
        self.ul_critic_tpList = []
        self.rvalues = []
        self.llhd = 0.0
        self.muhat = 1.0
        self.mumax = None
        self.TL = 0.0
        self.K = None
        self.letters = ""
        self.description = ""
        self.bestCombo = None
        self.decays: Dict[int, Dict] = {}
        self.masses: Dict[int, float] = {}
        self.possibledecays: Dict[int, list] = {}
        self.decay_keys: Dict[int, Dict] = {}
        self.inv_decay_keys: Dict[int, Dict] = {}
        self.decay_tuples: Dict[int, Dict] = {}
        self._stored_xsecs: Tuple = ()
        self._xsecMasses: Dict = {}
        self._xsecSSMs: Dict = {}
        self.ssmultipliers: Dict = {}
        ## Inititiaze LSP
        if False and self.environ.allowN1N1Prod:
            ## if we allow this, we might also start with this
            self.ssmultipliers[(self.LSP,self.LSP)]=1
            self.log ( f"we start with LSP,LSP production" )
        self.masses[ProtoModel.LSP] = float(np.random.uniform(100,500))
        self.decays[ProtoModel.LSP]= {}
        #pids = [(self.LSP,self.LSP)]                   #No (LSP,LSP) pair production
        pids = []
        if self.hasAntiParticle(self.LSP):
            pids += [(self.LSP,-self.LSP),(-self.LSP,-self.LSP)]
        for pidpair in pids:
            from scipy.stats import lognorm
            self.ssmultipliers[tuple(sorted(pidpair))]= float(lognorm.rvs(s = 1.0, scale = 1.0))

        slha_decay_keys = []

        with open ( self.environ.templateSLHA ) as slhaf:
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
            dtuples = {}
            dinvkey = {}
            for key in slha_decay_keys:
                if f"D{p}" in key[0]:
                    dpid,dpid2,dpid3,dpd = None,None,None,None
                    dtuple=tuple(map(int,key[0].split("_")[1:]))

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
                    if not key[0] in dinvkey:
                        dinvkey[ key[0] ] = set()
                    dinvkey[ key[0] ].add ( dpd )
                    dtuples.update({dpd: dtuple})

            self.possibledecays[p]=decays
            self.decay_keys[p] = dkey
            self.inv_decay_keys[p] = dinvkey
            self.decay_tuples[p] = dtuples

    def __str__(self) -> str:
        """Return a short human-readable summary of the model.

        Includes unfrozen particle names, K value, and likelihood ratio.
        """
        namer = SParticleNames ( susy=False )

        pNames = [namer.asciiName ( pid ) for pid in self.unFrozenParticles()]
        pNames = ','.join(pNames)
        pStr = f'ProtoModel ({pNames}):'
        if self.K:
            pStr += f' K = {self.K:1.2f}'
        else:
            pStr += f' K = {self.K}'
        if self.TL:
            pStr += f', TL = {self.TL:1.2f}'
        else:
            pStr += f', TL = {self.TL}'

        return pStr

    def __repr__(self):
        """ shortened version of __str__"""
        sK, sTL = formatObject(self.K), formatObject(self.TL)
        pStr = f'ProtoModel ({sK}, {sTL})'
        return pStr

    def hasAntiParticle(self, pid: int) -> bool:
        """Determine if particle *pid* needs its antiparticle considered
        in signal-strength multipliers.

        Self-conjugate particles (gluino, neutralinos, sneutrinos, etc.)
        return ``False``.

        :param pid: PDG particle ID.
        :returns: ``True`` if the antiparticle must be tracked.
        """
        _SELF_CONJUGATE = frozenset({
            1000021, 1000022, 1000023, 1000025, 1000035,
            1000012, 1000014, 1000016,
            2000012, 2000014, 2000016, 2000021,
        })
        return abs(pid) not in _SELF_CONJUGATE

    def toTuple(self, pid1: int, pid2: int) -> Tuple[int, int]:
        """Return *pid1*, *pid2* as a canonically sorted tuple."""
        return tuple(sorted((pid1, pid2)))

    def getXsecs(self) -> list:
        """Return the cross-sections.

        If they have already been computed (and stored in
        ``self._stored_xsecs``) **and** the masses and signal-strength
        multipliers have not been modified, return the cached value.
        Otherwise re-compute the cross-sections.

        :return: list of cross-sections
        """
        if self.masses == self._xsecMasses and self.ssmultipliers == self._xsecSSMs:
            if self._stored_xsecs:
                return self._stored_xsecs

        self.delXSecs()
        # Recompute — updates _xsecMasses, _xsecSSMs, and _stored_xsecs.
        self.computeXSecs()

        return self._stored_xsecs

    def getAllowedProdModes(self, return_mass : bool = False ):
        """
        Get the list of allowed production modes for the protomodel
        :param return_mass: If True, return the mass of the particles in the
        prod modes along with the prod modes

        :returns: List of all allowed production modes for the protomodel
        """

        tmpSLHA = tempfile.mktemp( prefix=f".{self.walkerid}_xsecfile",
                                   suffix=".slha",dir=self.SLHATEMPDIR )
        slhafile = self.createSLHAFile(tmpSLHA, addXsecs=False)
        channels = self.computer.findOpenChannels(slhafile)

        if return_mass: return channels

        prodModes = []
        for modes in channels:
            prodModes.append(modes['pids'])
        if len(prodModes) == 0:
            print(f"huh? we have 0 prod modes? We have {len(channels)} channels.")
        return prodModes

    def getOpenChannels(self, pid : int ):
        """Get the list of open decay channels for particle *pid*.

        Open channels are decays to unfrozen particles and to lighter
        particles.

        :param pid: PDG particle ID
        :return: list of daughter-pid tuples for each open decay channel
        """
        from base.constants import smMasses, smWidths

        openChannels: list = []
        unfrozen = self.unFrozenParticles()
        allMasses = {**self.masses, **smMasses}

        offshell = False
        if pid == 1000023 and pid in self.masses and self.LSP in self.masses and \
                (self.masses[pid] - self.masses[self.LSP]) < (smMasses["Z"] + smWidths["Z"]):
            offshell = True
        elif pid == 1000024 and pid in self.masses and self.LSP in self.masses and \
                  (self.masses[pid] - self.masses[self.LSP]) < (smMasses["W"] + smWidths["W"]):
            offshell = True
        else: offshell = False

        for dpid in self.possibledecays[pid]:
            #Get the list of BSM particles in the decay:
            if isinstance(dpid,(list,tuple)):
                pidList = [abs(p) for p in dpid if abs(p) in self.particles]
            else:
                self.highlight ( "warn", f"a decay channel without the SM particle is specified in {pid}:{str(dpid)}" )
                pidList = [abs(dpid)]
            #Skip decays to unfrozen particles
            if not all(dp in unfrozen for dp in pidList):
                continue
            #Get total daughter mass (it should only be a single mass)
            mdaughter = sum(allMasses.get(abs(p), 0) for p in dpid)

            #Skip decays to heavier particles
            if pid not in self.masses or mdaughter >= self.masses[pid]:
                continue

            if not offshell and len(dpid) == 3 and pid in (1000023, 1000024):
                continue

            openChannels.append(dpid)

        #remove all decay channels associated with a dkey if one of them is
        # not present for offshell decays to ensure flavor democracy
        if offshell:
            dk_groups: dict = {}
            for dpid, dk in self.decay_keys[pid].items():
                dk_groups.setdefault(dk, []).append(dpid)
            for dpid, dk in self.decay_keys[pid].items():
                if dpid in openChannels:
                    dec_not_present = [dc for dc in dk_groups.get(dk, []) if dc not in openChannels]
                    if dec_not_present:
                        self.highlight("warn", f"{dec_not_present} not in the open channels {openChannels} -- it's probably not open. but {dpid} is open. For now we will remove {dpid} from the open channels, ok?")
                        openChannels.remove(dpid)
                        self.highlight("info", f"Open channels are now {openChannels}")

        return openChannels

    def _invalidateParticleCache(self) -> None:
        """Mark the cached particle lists as stale."""
        self._unfrozen_cache = None
        self._frozen_cache = None
        self._masses_snapshot = None

    def _ensure_cache_valid(self) -> None:
        """Rebuild cache if masses dict has changed."""
        if self._masses_snapshot is not self.masses:
            self._invalidateParticleCache()
            self._masses_snapshot = self.masses

    def frozenParticles(self) -> List[int]:
        """Return PIDs of all particles that are *not* in the unfrozen list.

        Cached per masses-dict identity to avoid repeated O(n) scans.

        :returns: List of frozen particle IDs.
        """
        self._ensure_cache_valid()
        if self._frozen_cache is not None:
            return self._frozen_cache
        unfrozen = set(self.unFrozenParticles())
        self._frozen_cache = [pid for pid in self.particles if pid not in unfrozen]
        return self._frozen_cache

    def cleanBestCombo(self) -> None:
        """Remove unneeded data from bestCombo before storing."""
        if hasattr(self, "keep_meta") and self.keep_meta:
            return
        from tester.combiner import Combiner
        combiner = Combiner(self.walkerid)
        if hasattr(self, "bestCombo") and self.bestCombo is not None:
            self.bestCombo = combiner.removeDataFromBestCombo(self.bestCombo)

    def almostSameAs(self, other: "ProtoModel") -> bool:
        """Check if a model is essentially the same as *other*.

        Compares masses (relative tolerance 1e-5), signal-strength
        multipliers, and branching ratios.

        :param other: Another ProtoModel to compare against.
        :returns: ``True`` if the two models are indistinguishable.
        """
        if self.masses.keys() != other.masses.keys():
            return False

        massDiff = [abs(m - other.masses[pid]) / m for pid, m in self.masses.items() if m]
        if massDiff and max(massDiff) > 1e-5:
            return False

        # Compare signal-strength multipliers
        pidpairs = set(self.ssmultipliers.keys()) | set(other.ssmultipliers.keys())
        for pidpair in pidpairs:
            ss = self.ssmultipliers.get(pidpair, 1.0)
            os_val = other.ssmultipliers.get(pidpair, 1.0)
            if ss == 0.0:
                if os_val == 0.0:
                    continue
                return False
            if abs(ss - os_val) / ss > 1e-6:
                return False

        # Compare decays
        pids = set(self.decays.keys()) | set(other.decays.keys())
        for pid in pids:
            sdecays = self.decays.get(pid, {})
            odecays = other.decays.get(pid, {})
            dpids = set(sdecays.keys()) | set(odecays.keys())
            for dpid in dpids:
                sbr = sdecays.get(dpid, 0.0)
                obr = odecays.get(dpid, 0.0)
                if sbr == 0.0:
                    if obr < 1e-6:
                        continue
                    return False
                if abs(sbr - obr) / sbr > 1e-6:
                    return False
        return True

    def unFrozenParticles(self, withLSP: bool = True) -> List[int]:
        """Return PIDs of all particles with mass < 100 TeV.

        Results are cached and invalidated when ``self.masses`` is mutated.

        :param withLSP: If ``False``, exclude the LSP from the result.
        :returns: List of unfrozen particle IDs.
        """
        self._ensure_cache_valid()
        if self._unfrozen_cache is not None:
            ret = self._unfrozen_cache
        else:
            ret = [pid for pid, v in self.masses.items() if abs(v) < 1e5]
            self._unfrozen_cache = ret
        if not withLSP and self.LSP in ret:
            return [pid for pid in ret if pid != self.LSP]
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
        """ compute xsecs given the masses and signal strength multipliers of the
        model. The results are stored in self._stored_xsecs and should be
        accessed through getXsecs.

        :param keep_slha: if true, then keep slha file at the end
        :returns: current slha file, if slha file is kept, else none
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
                    self.pprint( f"error computing cross-sections: {e}, attempt # {countAttempts}" )
                    import traceback
                    traceback.print_stack()
                # helpers.cpPythia8()
                time.sleep ( np.random.uniform ( 5, 10 ) )
                if countAttempts > 5:
                    break

        if keep_slha:
            self.createSLHAFile( self.currentSLHA, addXsecs = True )
            return self.currentSLHA

    def rescaleXSecsBy(self, s : float, excl : list = [], cap_ssm = 100):
        """
        Rescale the stored cross-sections by a factor s
        :param excl: do not rescale xsecs/ssms of the prod modes in the list
        """
        #
        #if s > 100: return

        #Before rescaling, make sure we get the latest cross-sections:
        x = self.getXsecs()
        xsecs = x[0]
        comment = x[1]
        force_ssm = {}
        for k,v in self.ssmultipliers.items():
            if k not in excl:
                if v*s > cap_ssm:
                    force_ssm.update({k:v})
                    self.ssmultipliers[k] = cap_ssm       #do not let ssm > 100(?)
                else: self.ssmultipliers[k] = v * s
        for xsec in xsecs:
            if xsec.pid not in excl:
                if xsec.pid in force_ssm.keys(): xsec.value *= cap_ssm/force_ssm[xsec.pid]
                else: xsec.value *= s
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
        if not os.path.exists ( self.environ.templateSLHA ):
            if "/mnt/hephy/" in self.environ.templateSLHA:
                trySLHA = self.environ.templateSLHA.replace(f"{os.environ['CODEDIR']}/smodels-utils/protomodels/","./" )
                if os.path.exists ( trySLHA ):
                    self.environ.templateName = trySLHA
                    return

    def _writeSLHAFile ( self, outputSLHA : os.PathLike ):
        """ write the slha file, plug in protomodel params. this method does
        however not add the xsecs, for this look at createSLHAFile.

        :param outputSLHA: name of slha file to write
        """
        #Get template data:
        with open( self.environ.templateSLHA ) as f:
            lines=f.readlines()
        unfrozen = self.unFrozenParticles()
        # in "covered" we log that everything in self.decays
        # is covered in the slha file. else we complain
        covered = copy.deepcopy ( self.decays )
        covered.pop ( self.LSP ) ## no need to check
        ## in "inSLHAFile" we take note of all decays that are mentioned
        ## in the template slha file
        inSLHAFile = {}
        totalBRs = {}

        with open(outputSLHA,'wt') as outF:
            for i,l in enumerate(lines):
                for pid in self.particles:
                    #Skip lines which have no mass or decay tags
                    if not f"M{pid}" in l and not f"D{pid}" in l:
                        continue
                    if not pid in inSLHAFile:
                        inSLHAFile[pid]=set()
                        totalBRs[pid]=[]

                    #Get information for particle
                    if pid in unfrozen:
                        mass = self.masses[pid]
                        decays = {}
                        if pid in self.decays:
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
                        inSLHAFile[pid].add ( dpids )
                        if len(dpids) == 1:
                            dpids = dpids[0]
                        if dpids in decays:
                            if decayTag in self.inv_decay_keys[pid]:
                                for m_dpids in self.inv_decay_keys[pid][decayTag]:
                                    if m_dpids in covered[pid]:
                                        covered[pid].pop ( m_dpids )
                            #if dpids in covered[pid]:
                            #    print ( f"popping {dpids} from covered[{pid}]" )
                            #    covered[pid].pop ( dpids )
                            br = decays[dpids]
                            totalBRs[pid].append ( (dpids, br ) )
                            l = l.replace(decayTag, f"{br:.5f}" )
                        else:
                            l = ""

                #Only write line if it is not empty
                if l:
                    outF.write(l)
            remains = {}
            for pid,decays in covered.items():
                if len(decays)>0:
                    if not pid in remains:
                        remains[pid]={}
                    for dkey, dvalue in decays.items():
                        remains[pid][dkey]=dvalue
            if len(remains)>0:
                for pid,decays in remains.items():
                    for decay in decays:
                        self.error ( f"Protomodel lists a decay {pid} -> {decay}, but no equivalent found in template slha file!" )
                        import itertools
                        hasProposal = False
                        for d in itertools.permutations ( decay ):
                            if d in inSLHAFile[pid]:
                                self.error ( f"did you mean {pid} -> {d}?" )
                                hasProposal = True
                        if not hasProposal:
                            dpds = ", ".join ( map ( str, inSLHAFile[pid] ) )
                            self.error ( f"channels mentioned in template file: {dpds}" )
                sys.exit(-1)
            for pid, allbrs in totalBRs.items():
                totalbr = sum( [ x[1] for x in allbrs ])
                if abs(totalbr-0.) > 1e-3 and abs(totalbr-1.) > 1e-3:
                    self.error ( f"total brs for {pid} add up to {totalbr:.3f} != 1." )
                    self.error ( f"contributions are:" )
                    for dpids_br in totalBRs[pid]:
                        self.error ( f"{dpids_br[0]}: {dpids_br[1]}" )
                    sys.exit(-1)
            frozen = self.frozenParticles()
            # now make the frozen particles stable (to quench smodels warnings,
            # nothing else)
            outF.write ( "\n" )
            for fpid in frozen:
                line = f"DECAY   {fpid}    {0E+00}\n\n"
                outF.write ( line )
            outF.close()

    def createSLHAFile ( self, outputSLHA : Union[str,None] = None,
                         addXsecs : bool = True ) -> str:
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
        self._writeSLHAFile ( outputSLHA )

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
                        self._writeSLHAFile ( outputSLHA )
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
                    time.sleep ( np.random.uniform ( 0.5, 2.*ctAttempts ) )

        return outputSLHA

    def dict(self, sort_dict: bool = False) -> Dict:
        """Return a JSON-serialisable dictionary of the model state.

        :param sort_dict: If ``True``, sort keys for deterministic output.
        :returns: Dictionary with ``masses``, ``ssmultipliers``, ``decays``,
            and ``xsecs[fb]``.
        """
        xsecs = {}
        tmp = self.getXsecs()
        if len(tmp) > 0:
            for xsec in tmp[0]:
                xsecs[(xsec.pid, xsec.info.sqrts.asNumber(TeV))] = xsec.value.asNumber(fb)
        if sort_dict:
            return {
                "masses": {pid: self.masses[pid] for pid in sorted(self.masses)},
                "ssmultipliers": {pp: self.ssmultipliers[pp] for pp in sorted(self.ssmultipliers)},
                "decays": {
                    pid: {dpid: self.decays[pid][dpid] for dpid in sorted(self.decays[pid])}
                    for pid in sorted(self.decays)
                },
                "xsecs[fb]": xsecs,
            }
        return {
            "masses": self.masses,
            "ssmultipliers": self.ssmultipliers,
            "decays": self.decays,
            "xsecs[fb]": xsecs,
        }

    def relevantSSMultipliers(self) -> Dict:
        """Return only signal-strength multipliers for unfrozen particles
        with values that deviate from unity.

        :returns: Filtered SSM dictionary.
        """
        frozen = set(self.frozenParticles())
        return {
            pids: v
            for pids, v in self.ssmultipliers.items()
            if abs(v - 1.0) >= 1e-5 and not any(abs(pid) in frozen for pid in pids)
        }

    def describe(self) -> None:
        """Print a brief summary of the model contents."""
        ndecays, nd = 0, 0
        for k, v in self.decays.items():
            if k == ProtoModel.LSP:
                continue
            ndecays += len(v)
            nd += 1
        nssms = len(self.ssmultipliers)
        print(f"{len(self.masses)} masses, {ndecays}[{nd}] decays, {nssms} ss multipliers")

    def delXSecs ( self ):
        """ delete stored cross section, if they exist """
        self._stored_xsecs = ()
        self._xsecMasses = {}
        self._xsecSSMs = {}

    def copy(self, cp_predictions: bool = False) -> "ProtoModel":
        """Create a deep copy of this model.

        :param cp_predictions: If ``True``, also copy ``bestCombo`` and
            ``ul_critic_tpList`` via ``deepcopy``.
        :returns: A new ProtoModel instance with identical state.
        """
        newmodel = self.__class__(self.walkerid, self.keep_meta, self.environ)
        newmodel.keep_meta = self.keep_meta
        newmodel.maxMass = self.maxMass
        newmodel.step = self.step
        newmodel.protomodels_version = self.protomodels_version
        newmodel.forced_degeneracies = self.forced_degeneracies
        newmodel.decaylessParticles = self.decaylessParticles
        newmodel.particles = self.particles[:]
        newmodel.environ = self.environ
        newmodel.possibledecays = dict(self.possibledecays)
        newmodel.decays = {pid: dict(dec) for pid, dec in self.decays.items()}
        newmodel.masses = dict(self.masses)
        newmodel.ssmultipliers = dict(self.ssmultipliers)
        newmodel.rvalues = self.rvalues[:]
        newmodel.llhd = self.llhd
        newmodel.muhat = self.muhat
        newmodel.mumax = self.mumax
        newmodel.TL = self.TL
        newmodel.K = self.K
        newmodel.letters = self.letters[:]
        newmodel.description = self.description[:]
        newmodel._stored_xsecs = copy.deepcopy(self._stored_xsecs)
        newmodel._xsecSSMs = dict(self._xsecSSMs)
        newmodel._xsecMasses = dict(self._xsecMasses)
        if cp_predictions:
            newmodel.bestCombo = copy.deepcopy(self.bestCombo)
            newmodel.ul_critic_tpList = copy.deepcopy(self.ul_critic_tpList)
            newmodel.llhd_critic_preds = copy.deepcopy(self.llhd_critic_preds)
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
