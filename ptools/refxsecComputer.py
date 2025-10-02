#!/usr/bin/env python3

"""
.. module:: refxsecComputer
   :synopsis: code to retrieve cross sections from LHC Xsec group's
              references, published at
              https://twiki.cern.ch/twiki/bin/view/LHCPhysics/SUSYCrossSections

.. moduleauthor:: Wolfgang Waltenberger <wolfgang.waltenberger@gmail.com>

"""

from __future__ import print_function
from smodels.tools.wrapperBase import WrapperBase
from smodels.tools import wrapperBase
# smodels v2
# from smodels.tools.smodelsLogging import logger, setLogLevel
# from smodels.tools.physicsUnits import fb, pb, GeV, TeV, mb, unum
# from smodels.theory.crossSection import LO, NLO, NLL, NNLL
# from smodels.theory import crossSection
# from smodels.theory.exceptions import SModelSTheoryError as SModelSError
from smodels.base.physicsUnits import fb, pb, GeV, TeV, mb, unum
from smodels.base.crossSection import LO, NLO, NLL, NNLL
from smodels.base import crossSection
from smodels.base.exceptions import SModelSBaseError as SModelSError
from smodels import installation as smodelsinstallation
import os, sys, io, shutil, pyslha
from base.loggerbase import LoggerBase
from typing import Union


class RefXSecComputer ( LoggerBase ):
    """
    The xsec computer that simply looks up reference cross sections,
    and interpolates them.
    """
    version = "1.0" ## make sure we can trace changes in the tables
    hasWarned = { "omitted": 0 }

    def __init__( self, verbose : int = 1, allowN1N1Prod : bool = False, walkerid : Union[None,int] = 0 ):
        """
        :param verbose: turn on verbose mode, for debugging
        :param allowN1N1Prod: if true, then allow also N1 N1 production
        """
        super ( RefXSecComputer, self ).__init__ ( walkerid )
        self.verbose = verbose
        self._allowN1N1Prod = allowN1N1Prod
        codedir = "../"
        try:
            import csetup
            codedir = csetup.setup()
        except ImportError as e:
            pass
        self.shareDir = f"{codedir}/ptools/xsecTables/"
        # productions of same-sign-pid pairs when the particle is within reach
        self.samesignmodes = [ 1000012, 1000014, 1000016, 1000021, 1000023, 1000025 ]
        self.addN1N1ProdMode()
        # production of opposite-sign-pid pairs when the particle is within reach
        self.oppositesignmodes = ( 1000001, 1000002, 1000003, 1000004, 1000005, 
                2000005, 1000006, 2000006, 1000011, 1000013, 1000015, 1000024, 
                1000037 )

        # associate production
        #IMP! When adding prod modes, make sure pid[0] < pid[1]!!
        self.associateproduction = ( ( 1000001, 1000021 ), ( 1000002, 1000021),
                ( 1000003, 1000021 ), ( 1000004, 1000021 ), ( 1000005, 1000021),
                ( 1000021, 2000005 ), ( 1000006, 1000021 ), ( 1000021, 2000006),
                ( -1000012, 1000011), ( -1000014, 1000013), (-1000016, 1000015),
                ( -1000011, 1000012 ), ( -1000013, 1000014 ), ( -1000015, 1000016),             
                #removing charge for charginos in chargino-neutralino prod as doesnt matter
                ( 1000022, 1000023 ), ( 1000023, 1000025 ),
                # ( 1000022, 1000024 ), #( -1000024, 1000022),
                # ( 1000023, 1000024 ),  #( -1000024, 1000023 ),
                #( 1000023, 1000037 ), ( 1000024, 1000025 ), ( -1000024, 1000037 ),
                #( -1000037, 1000024 ), ( 1000025, 1000037 ) ) #, ( -1000037, 1000023 ),
                )
                #( -1000024, 1000025 ), ( -1000037, 1000025 )

        ## these are productions where we explicitly disallow the negative sign versions,
        ## e.g. C1+,N2 ## this is based on the assumption that our database of chargino results
        ## is entirely charge symmetrical!! this might change with e.g. protomodels v3
        self.associatesamesignproduction = ( ( 1000022, 1000024 ), 
                ( 1000023, 1000024 ), ( 1000024, 1000025 ), 
                ( 1000022, 1000037 ), ( 1000023, 1000037 ), 
                ( 1000025, 1000037 ) )
                 
        # self.schannel = ( 35, 55, )
        self.schannel = tuple()

    @property
    def allowN1N1Prod(self):
        return self._allowN1N1Prod

    @allowN1N1Prod.setter ## convenience
    def allowN1N1Prod(self,flag : bool ):
        self._allowN1N1Prod = flag
        self.addN1N1ProdMode()

    def addN1N1ProdMode ( self ):
        """ if allowN1N1Prod, add the mode """
        if not self.allowN1N1Prod:
            return
        self.samesignmodes = tuple ( list ( self.samesignmodes ) + [ 1000022 ] )

    def warn ( self, *txt ):
        stxt=str(*txt)
        if not stxt in self.hasWarned:
            self.hasWarned[stxt]=0
        self.hasWarned[stxt]+=1
        if self.hasWarned[stxt]<3:
            self.warning ( *txt )
        if self.hasWarned[stxt]==3:
            self.hasWarned["omitted"]+=1
            if self.hasWarned["omitted"]<2:
                self.warning ( "(omitted similar msgs)" )

    def checkFileExists(self, inputFile):
        """
        Check if file exists, raise an IOError if it does not.
        :returns: absolute file name if file exists.

        """
        nFile = self.absPath(inputFile)
        if not os.path.exists(nFile):
            raise IOError(f"file {nFile} does not exist")
        return nFile

    def dictToXSection ( self, D ):
        """ create an XSection object from dictionary D """
        xsec = crossSection.XSection()
        xsec.value = D["xsec"]*pb
        xsec._pid = D["pids"]
        xsec.info = crossSection.XSectionInfo ( D["sqrts"]*TeV, D["order"], D["label"] )
        xsec.comment = D["comment"]
        return xsec

    def lastLineShouldBeEmpty ( self, slhafile ):
        """ the last line in the slhafile should be an empty line.
            append one, if need be. """
        with open ( slhafile, "rt" ) as f:
            lines=f.readlines()
            f.close()
        lastline = lines[-1].strip()
        if len(lastline)>1:
            with open ( slhafile, "at" ) as f:
                f.write ( "\n" )
                f.close()

    def addXSecToFile( self, xsecs, slhafile, complain=True):
        """
        Write cross sections to an SLHA file, analogous to xsecComputer

        :param xsecs: a XSectionList object containing the cross sections
        :param slhafile: target file for writing the cross sections in SLHA format
        :param complain: complain if there are already cross sections in file

        """

        if not os.path.isfile(slhafile):
            line = f"SLHA file {slhafile} not found."
            self.error( line )
            raise SModelSError( line )
        if len(xsecs) == 0:
            self.debug(f"No cross sections available for {slhafile}.")
            return False
        self.debug ( f"I have xsecs: {', '.join(map(str,xsecs))}" )
        # Check if file already contain cross section blocks
        xSectionList = crossSection.getXsecFromSLHAFile(slhafile)
        if xSectionList and complain:
            self.info("SLHA file already contains XSECTION blocks. Adding "
                           "only missing cross sections.")

        # Write cross sections to file, if they do not overlap any cross section in
        # the file
        self.lastLineShouldBeEmpty ( slhafile )

        outfile = open(slhafile, 'a')
        nxsecs = 0
        for xsec in xsecs:
            xseccomment = f"reference xsecs v{self.version} [pb]"
            writeXsec = True
            # print ( "in addXSecToFile comment", xsec, hasattr ( xsec, "comment" ) )
            if hasattr ( xsec, "comment" ) and xsec.comment not in [ None, "", "None", " (None)" ]:
                xseccomment += f" {xsec.comment}"
            for oldxsec in xSectionList:
                if oldxsec.info == xsec.info and set(oldxsec.pid) == set(xsec.pid):
                    writeXsec = False
                    break
            if writeXsec:
                nxsecs += 1
                outfile.write( f"{self.xsecToBlock(xsec, (2212, 2212), xseccomment)}\n")
        outfile.close()

        return nxsecs

    def xsecToBlock( self, xsec, inPDGs=(2212, 2212), comment=None, xsecUnit = pb):
        """
        Generate a string for a XSECTION block in the SLHA format from a XSection
        object.

        :param inPDGs: defines the PDGs of the incoming states
                       (default = 2212,2212)

        :param comment: is added at the end of the header as a comment
        :param xsecUnit: unit of cross sections to be written (default is pb).
                         Must be a Unum unit.

        """
        if type(xsec) != type(crossSection.XSection()):
            self.error("Wrong input")
            raise SModelSError()
        # Sqrt(s) in GeV
        header = f"XSECTION  {xsec.info.sqrts / GeV!s}"
        for pdg in inPDGs:
            # PDGs of incoming states
            header += f" {pdg!s}"
        # Number of outgoing states
        xsec.pid = tuple ( [ x for x in xsec.pid if x != None ] )
        header += f" {len(xsec.pid)!s}"
        for pid in xsec.pid:
            # PDGs of outgoing states
            header += f" {pid!s}"
        if comment:
            header += f" # {comment!s}"  # Comment
        entry = f"  0  {xsec.info.order!s}  0  0  0  0  {float(xsec.value / xsecUnit):16.8E} SModelSv{smodelsinstallation.version()}"

        return f"\n{header}\n{entry}"

    def computeForOneFile ( self, sqrtses, inputFile,
                 tofile, ssmultipliers = None, comment = None,
                 ignore_pids = None, ewk = None ):
        """
        Compute the cross sections for one file.

        :param sqrtses: list of sqrt{s} tu run pythia, as a unum (e.g. [7*TeV])
        :param inputFile: input SLHA file to compute xsecs for
        :param tofile: False, True, "all": write results to file,
                       if "all" also write lower xsecs to file.
        :param ssmultipliers: optionally supply signal strengh multipliers,
                          given as dictionary of the tuple of the mothers' pids as keys and
                          multipliers as values, e.g { (1000001,1000021):1.1 }.
        :param comment: an optional comment that gets added to the slha file.
        :param ignore_pids: if not None, pids to ignore when computing xsecs
        :param ewk: if hino and C1C1 is asked, use hino xsecs

        :returns: number of xsections that have been computed
        """
        self.xsecs = crossSection.XSectionList()
        nXSecs = 0
        if tofile:
            self.info(f"Computing SLHA cross section from {inputFile}, adding to SLHA file." )
            complain = True ## dont complain about already existing xsecs,
            # if we were the ones writing them
            for s in sqrtses:
                # ss = s*TeV
                self.compute( s, inputFile, ssmultipliers = ssmultipliers,
                              ignore_pids = ignore_pids, ewk = ewk )
                if tofile == "all" and hasattr ( self, "loXSecs" ):
                    nXSecs += self.addXSecToFile(self.loXsecs, inputFile, complain )
                    complain = False
                if tofile != False:
                    ## FIXME check if higher orders are already in from
                    ## ref xsecs
                    # print ( "xsecs", self.xsecs, hasattr ( self.xsecs[0], "comment" ) )
                    nXSecs += self.addXSecToFile( self.xsecs, inputFile, complain)
                    complain = False
            if nXSecs > 0: ## only add if we actually added xsecs
                self.addMultipliersToFile ( ssmultipliers, inputFile )
            self.addCommentToFile ( comment, inputFile )
            self.cleanSLHAFile ( inputFile )
        else:
            self.info(f"Computing SLHA cross section from {inputFile}." )
            print()
            print( "     Cross sections:" )
            print( "=======================" )
            for s in sqrtses:
                ss = s*TeV
                self.compute( ss, inputFile, unlink=unlink, loFromSlha=lOfromSLHA,
                              ssmultipliers = ssmultipliers )
                for xsec in self.xsecs:
                    nXSecs += 1
                    print( f"{xsec.info.label} {xsec.pid:>20}:  {xsec.value / pb:.3e} pb" )
            print()
        return nXSecs

    def cleanSLHAFile ( self, slhafile ):
        """ clean up F, if needed. remove double newlines """
        f=open(slhafile,"rt")
        lines=f.readlines()
        f.close()
        f=open(slhafile,"wt")
        newline = False
        for line in lines:
            #if "Signal strength multipliers" in line:
            #    continue
            if newline and line == "\n":
                continue
            if line == "\n":
                newline = True
            else:
                newline = False
            f.write ( line )
        f.close()

    def computeForBunch ( self, sqrtses, inputFiles, tofile, ssmultipliers=None ):
        """ compute xsecs for a bunch of slha files """
        for inputFile in inputFiles:
            self.debug ( f"computing xsec for {inputFile}" )
            self.computeForOneFile ( sqrtses, inputFile, tofile,
                                     ssmultipliers = ssmultipliers )

    def addCommentToFile ( self, comment, slhaFile ):
        """ add the optional comment to file """
        if comment in [ None, "" ]:
            return
        if not os.path.isfile(slhaFile ):
            self.error(f"SLHA file {slhaFile} not found." )
            raise SModelSError()
        outfile = open(slhaFile, 'a')
        outfile.write ( f"\n# {comment}\n" )
        outfile.close()

    def addMultipliersToFile ( self, ssmultipliers, slhaFile ):
        """ add the signal strength multipliers to the SLHA file """
        if ssmultipliers in [ None, {} ]:
            return
        if not os.path.isfile(slhaFile ):
            self.error(f"SLHA file {slhaFile} not found." )
            raise SModelSError()
        tokens = []
        for k,v in ssmultipliers.items():
            tokens.append ( f"{k}:{v:.4g}" )
        newline = f"# Signal strength multipliers: {', '.join(tokens)}"
        with open(slhaFile, 'r' ) as r:
            lines = r.readlines()
            r.close()

        rewrite = []
        for line in lines:
            if "Signal strength multipliers" in line:
                if ( line.strip() == newline ):
                    self.debug ( "Signal strength multipliers have alread been applied." )
                else:
                    self.error ( "Different signal strength multipliers have alread been applied!!!" )
                    rewrite.append ( f"{line} ERROR inconsistent!" )
            else:
                if not "produced at step" in line:
                    rewrite.append ( line )
        outfile = open(slhaFile, 'w')
        for line in rewrite:
            outfile.write ( line )
        if line != "\n": ## last line not an empty newline?
            outfile.write ( "\n" )
        outfile.write ( newline )
        outfile.write ( "\n" )
        outfile.close()

    def isIn ( self, pids, ignores ):
        """ are pids in the ignores? """
        anyTuples = False
        if type(ignores) in [ tuple, list ]:
            for i in ignores:
                if type(i)==tuple:
                    anyTuples = True
                    break
        if anyTuples:
            # we have a list of ignores, needs to be in any one
            for i in ignores:
                if self.isIn ( pids, i ):
                    return True
            return False
        if type(pids) == list:
            pids = tuple(pids)
        if type(ignores ) == list:
            ignores = tuple(ignores)
        if type(ignores) == int:
            ignores = ( ignores, )
        if len(ignores)==1:
            if ignores[0] in pids or -ignores[0] in pids:
                return True
            return False
        if pids == ignores:
            return True
        ## reverse?
        if ( pids[1], pids[0]) == ignores:
            return True
        ## sign flip
        if ( -pids[0], pids[1]) == ignores:
            return True
        if ( pids[0], -pids[1]) == ignores:
            return True
        return False

    def selectChannels ( self, channels, ignore ):
        """ from channels, filter out all with pids that are in ignore """
        if ignore in [ None, [], "" ]:
            return channels
        ret = []
        try:
            myignores = eval(ignore)
        except Exception as e:
            self.error ( f"dont understand the ignore arg {ignore}: {e}" )
            sys.exit()
        for c in channels:
            tobeignored = False
            isin = self.isIn ( c["pids"], myignores )
            if isin:
                tobeignored = True
            if not tobeignored:
                self.info ( f"selecting {c}" )
                ret.append ( c )
        return ret

    def getPythiaCard( self, pids ):
        """
        Gives path to the pythia config file, modified to compute only the relevant channels
        :param pids: IDs of particles to produce
        :returns: path to the pythia config file
        """
        import tempfile
        codedir = "../"
        try:
            import csetup
            codedir = csetup.setup()
        except ImportError as e:
            pass
        originalConfig = f"{codedir}/ptools/pythia8.cfg"
        temporaryConfig = tempfile.mktemp( prefix=f"pythia8_", suffix=".cfg", dir="/tmp/" )
        with open(originalConfig,'r') as file:
            lines = file.readlines()
        for i,line in enumerate(lines):
            if line.split() and line.split()[0] == "SUSY:idA":
                lines[i] = f'SUSY:idA = {int(abs(pids[0]))}           ! 0: all\n'
            if line.split() and line.split()[0] == "SUSY:idB":
                lines[i] = f'SUSY:idB = {int(abs(pids[1]))}           ! 0: all\n'
        with open(temporaryConfig,'w') as file:
            for line in lines:
                file.write(line)
        return temporaryConfig

    def compute( self, sqrts, slhafile, ssmultipliers = None,
                 ignore_pids = None, ewk = "wino" ):
        """
        Retrieve cross sections
        :param sqrts: center of mass energies
        :param slhafile: SLHA file
        :param ssmultipliers: optionally supply signal strengh multipliers,
                given as dictionary of the tuple of the mothers' pids as keys and
                multipliers as values, e.g { (1000001,1000021):1.1 }.
        :param ignore_pids: ignore pids for production
        :param ewk: 
        :returns: List of cross sections to be added
        """
        channels = self.findOpenChannels ( slhafile )
        #print ( "open", [ c["pids"] for c in channels ] )
        channels = self.selectChannels ( channels, ignore_pids )
        #print ( "selected ", channels )
        xsecs = crossSection.XSectionList()
        for channel in channels:
            # obtain xsecs for all masses, but for the given channel
            pids = channel["pids"]
            if pids[1]!=None and pids[1] < pids[0]:
                pids = ( pids[1], pids[0] )
            # obtain xsecs for all masses, but for the given channel
            xsecall,order,comment = self.getXSecsFor ( pids[0], pids[1],
                    sqrts, ewk, channel["masses"] )
            # print ( f"for channel {pids}: {str(xsecall)[:10]}" )
            ## interpolate for the mass that we are looking for
            if xsecall == None:
                continue
            xsec = self.interpolate ( channel["masses"], xsecall )
            if xsecall is None:
                self.debug (f"No cross section table found for {channel}")
            if xsec == None:
                try:
                    self.warning (f'No cross section obtained for channel {channel} at {sqrts} TeV. Will try with PYTHIA.')
                    from smodels.tools.xsecComputer import XSecComputer, NLL
                    xsecComputer = XSecComputer ( NLL, 5000, pythiaVersion=8, maycompile=True )
                    pythia = xsecComputer.getPythia()
                    pythia.tempdir = None # reset, to make sure it works in parallel mode
                    pythia.nevents = xsecComputer.nevents
                    pythia.sqrts = sqrts
                    pythia.pythiacard = self.getPythiaCard ( pids )
                    loXsecs = pythia.run( slhafile, unlink=True )
                    xsecComputer.loXsecs = loXsecs
                    xsecComputer.loXsecs.sort()
                    xsecComputer.xsecs = xsecComputer.addHigherOrders ( sqrts * TeV, slhafile )
                    xsecComputer.xsecs.sort()
                    for x in xsecComputer.loXsecs: # If no NLL xsec, use LO xsec instead
                        if set(pids) == set(x.pid):
                            xsec = x.value.asNumber(pb)
                    for x in xsecComputer.xsecs:
                        if set(pids) == set(x.pid):
                            xsec = x.value.asNumber(pb)
                    comment += " Computed with PYTHIA"
                    os.unlik(pythia.pythiacard)
                    pythia.pythiacard = None
                    self.warning (f'PYTHIA cross section: {xsec}')
                except Exception as e:
                    if pythia.pythiacard != None:
                        os.unlink(pythia.pythiacard)
                    pythia.pythiacard = None
                    self.error (f'PYTHIA computation failed for channel {channel} at {sqrts} TeV, the reason is the following: {e}')
            if xsec == None:
                self.debug ('Failed also with PYTHIA. No cross section returned')
                continue
            if ssmultipliers != None and ( pids[1], pids[0] ) in ssmultipliers:
                pids = ( pids[1], pids[0] )
            if ssmultipliers != None and pids in ssmultipliers:
                ssm = ssmultipliers[pids]
                channel["ssm"] = ssm
                xsec = xsec * ssm
                channel["xsec"] = xsec
                channel["sqrts"] = sqrts
                channel["order"] = order
                channel["comment"] = comment
                orderStr = crossSection.orderToString(order,False,False)
                channel["label"] = f"{int(sqrts)} TeV ({orderStr})"
                a = self.dictToXSection ( channel )
                a.comment = comment
                # print ( "adding", a, hasattr ( a, "comment" ) )
                xsecs.add ( a )
            else:
                self.debug (f"No signal strength multiplier for {pids}")
                self.debug (pids,pids in ssmultipliers)
                self.debug (ssmultipliers)
            
        # print ( "xdding", xsecs, hasattr ( xsecs[0], "comment" ) )
        self.xsecs = xsecs
        #if len(self.xsecs)>0:
        #    self.xsecs[0].comment = comment
        # print ( "xdding", self.xsecs, hasattr ( self.xsecs[0], "comment" ) )

    def findOpenChannels ( self, slhafile ):
        slhadata = pyslha.readSLHAFile ( slhafile )
        masses = slhadata.blocks["MASS"]
        # print ( "findOpenChannels" )
        channels = []

        for pid,mass in masses.items():
            if pid < 999999 and pid not in self.schannel:
                continue
            if type(mass) not in [ float, int ]:
                self.error ( f"I found a mass of {mass} in {slhafile}, do not know what to do with it." )
                sys.exit(-1)
            if mass > 5000:
                continue
            if pid in self.schannel:
                channels.append ( { "pids": (pid,None), "masses": ( mass,None ) } )

            if pid in self.samesignmodes:
                channels.append ( { "pids": (pid,pid), "masses": ( mass, mass ) } )
            if pid in self.oppositesignmodes:
                channels.append ( { "pids": (-pid,pid), "masses": ( mass, mass ) } )
            for jpid, jmass in masses.items():
                if pid >= jpid or jpid < 999999 or jmass > 5000:
                    continue
                if (pid,jpid) in self.associateproduction:
                    channels.append ( { "pids": (pid,jpid), "masses": (mass, jmass ) } )
                if -pid < jpid and (-pid,jpid) in self.associateproduction:
                    channels.append ( { "pids": (-pid,jpid), "masses": (mass, jmass ) } )
                if jpid < -pid and (jpid,-pid) in self.associateproduction:
                    channels.append ( { "pids": (jpid,-pid), "masses": (jmass, mass ) } )
                    ppair = (-pid,jpid); masspair = (mass,jmass)
                if pid < -jpid and (pid,-jpid) in self.associateproduction:
                    channels.append ( { "pids": (pid,-jpid), "masses": (mass, jmass ) } )
                if -jpid < pid and (-jpid,pid) in self.associateproduction:
                    channels.append ( { "pids": (-jpid,pid), "masses": (jmass, mass ) } )
                if (abs(pid),abs(jpid)) in self.associatesamesignproduction:
                    channels.append ( { "pids": (abs(pid),abs(jpid)), "masses": (mass, jmass ) } )

        if len(channels)==0:
            print ( f"[refxsecComputer] found no open channels for {slhafile}" )
        return channels

    def cast ( self, xseckeys ):
        xsecs = []
        for k in xseckeys:
            if type(k) == float:
                xsecs.append ( k )
            if type(k) == str:
                xsecs.append ( float(k) )
            if type(k) == tuple:
                t = []
                for v in k:
                    if type(v)==float:
                        t.append ( v )
                    if type(v)==str:
                        t.append ( float(v) )
                xsecs.append ( tuple(t) )
        return xsecs

    def outOfBounds ( self, mass, xsecs ):
        """ check if masses are out of bounds """
        if type(mass) in [ int, float ]:
            if mass > max(xsecs):
                self.info ( f"mass {int(mass)}>{int(max(xsecs))} too high to interpolate, leave it as is." )
                return True
            if mass < min(xsecs):
                self.info ( f"mass {int(mass)}<{int(min(xsecs))} too low to interpolate, leave it as is." )
                return True
            return False
        ## masses are tuple
        for i,mi in enumerate(mass):
            if type(xsecs[0]) in [ float, int ]:
                xi = xsecs
            else:
                xi = [ x[i] for x in xsecs ]
            if mi < min(xi):
                self.info ( f"{i}st mass {mi}<{min(xi)} too low to interpolate, leave it as is." )
            if mi > max(xi):
                self.info ( f"{i}st mass {mi}>{max(xi)} too high to interpolate, leave it as is." )
                return True
        return False

    def collapse ( self, mass ):
        """ ??? """
        if type(mass) in [ int, float ]:
            return mass
        if type(mass) in [ str ]:
            return float(mass)
        if type(mass) in [ tuple, list ]:
            mass = [ x for x in mass if x != None ]
        # for i in range(len(mass)-1):
        #     smass = mass[i]+mass[i+1]
        #     if smass > 1e-6 and abs (mass[i]-mass[i+1]) / smass > 1e-3:
        #         return mass
        # return mass[0]
        return sum(mass)/len(mass) # Use mean of masses (if we don't use non-degenerate cross sections)

    def interpolate ( self, mass, xsecs ):
        """ interpolate between masses """
        if len ( xsecs ) == 0:
            self.error ( f"cannot interpolate empty set" )
            return None
        mass = self.collapse ( mass )
        if mass in xsecs:
            return xsecs[mass]
        xsecks = self.cast ( xsecs.keys() )
        if self.outOfBounds ( mass, xsecks ):
            return None
        if type(mass)==float:
            from scipy.interpolate import interp1d
            return interp1d ( xsecks, list(xsecs.values()) )( mass )
        from scipy.interpolate import griddata
        ret = griddata ( xsecks, list(xsecs.values()), mass )
        return ret

    def getXSecsFrom ( self, path, pb = True, columns={"mass":0,"xsec":1 } ):
        """ retrieve xsecs from filename
        :param pb: xsecs given in pb
        :param indices: the indices of the columns in the table, for mass and xsec
        """
        ret = {}
        if not os.path.exists ( path ):
            self.info ( f"could not find {path}" )
            return ret
        self.info ( f"getting xsecs from {path}" )
        f = open ( path, "rt" )
        lines=f.readlines()
        f.close()
        # print ( "get", columns )
        for line in lines:
            if line.startswith("#"):
                continue
            tokens = line.split ()
            if len(tokens)<2:
                continue
            m = columns["mass"]
            if type(m) in [ int ]:
                mass = float ( tokens[ m ] )
            if type(m) in [ list, tuple ]:
                mass = tuple( tokens[x] for x in m )
                if len(mass) == 1:
                    mass = float(mass)
            xsec = float(tokens[ columns["xsec"] ].replace("GeV","") )
            if not pb:
                xsec = xsec / 1000.
            ret[ mass ] = xsec
        return ret

    def getXSecsFor ( self, pid1, pid2, sqrts, ewk, masses ):
        """ get the xsec dictionary for pid1/pid2, sqrts
        :param ewk: specify the ewkino process (hino, or wino, or None)
        """
        self.debug ( f"asking for cross sections for pids={pid1,pid2}, {sqrts} TeV" )
        filename = None
        order = 0
        pb = True
        columns = { "mass": 0, "xsec": 1 }
        isEWK = False
        comment = ""
        # comment="refxsec [pb]"
        if pid1 in [ 35 ] and pid2 == None:
            #if sqrts == 8: # we dont have xsecs for that
            #    return None, None, None
            filename = f"xsecScalar{sqrts}.txt"
            columns["xsec"]=1
            isEWK=False
            order = LO
        if pid1 in [ 55 ] and pid2 == None:
            #if sqrts == 8: # we dont have xsecs for that
            #    return None, None, None
            filename = f"xsecVector{sqrts}.txt"
            columns["xsec"]=1
            isEWK=False
            order = LO
        if pid1 in [ 1000021 ] and pid2 == pid1: # Gluino pair production
            filename = f"xsecgluino{int(sqrts)}.txt"
            columns["xsec"] = 2
            order = NNLL # 4
            if masses[0] < 500:
                order = NLL
        if pid1 in [ 1000001, 1000002, 1000003, 1000004, 1000005, 2000005, 1000006, 2000006 ] and pid2 in [ 1000021 ]: # Gluino-squark productions
            filename = f"xsecsquark{int(sqrts)}.txt"
            columns["xsec"] = 2
            order = NNLL
        if pid1 in [ -1000001, -1000002, -1000003, -1000004 ] and pid2 == -pid1: # Squark-antisquark productions (light quarks)
            filename = f"xsecsquark{int(sqrts)}.txt"
            columns["xsec"] = 2
            order = NNLL
        if pid1 in [ -1000005, -2000005, -1000006, -2000006 ] and pid2 == -pid1: # Squark-antisquark productions (heavy quarks)
            filename = f"xsecstop{int(sqrts)}.txt"
            columns["xsec"] = 2
            order = NNLL
        if pid1 in [ -1000011, -1000012, -1000013, -1000014, -1000015, -1000016 ] and pid2 in [ 1000011, 1000012, 1000013, 1000014, 1000015, 1000016 ]: # Slepton pair productions
            filename = f"xsecslepslep{int(sqrts)}.txt"
            if sqrts == 8:
                pb == False
        if pid1 in [ -1000024, -1000037 ] and pid2 in [ 1000022, 1000023, 1000025 ]: # Charginos(-) neutralinos production
            filename = f"xsecN2C1m{int(sqrts)}.txt"
            order = NLL
            isEWK = True
            if sqrts == 13:
                pb = False
        if (pid1 in [ 1000022, 1000023, 1000025 ] and pid2 in [ 1000024, 1000037 ]) or (pid1 == 1000024 and pid2 == 1000025): # Charginos(+) neutralinos productions -- 'or' condition because pid1 < pid2
            filename = f"xsecN2C1p{int(sqrts)}.txt"
            order = NLL
            if sqrts == 13:
                pb = False
            isEWK = True
        N1N2N3 = ( 1000022, 1000023, 1000025 )
        N2N3 = ( 1000023, 1000025 )
        neutralinos = N2N3
        if self.allowN1N1Prod:
            neutralinos = N1N2N3
        if pid1 in N1N2N3 and pid2 in neutralinos: # Neutralinos pair production
            if masses[1]+masses[0] == 0.:
                self.info ( f"Asking for massless {(pid1,pid2)} production. Will return None." )
                return None, None, None
            self.info ( f"Asking for {(pid1,pid2)} production. Will use C1 C1 xsecs (No N1 N2 cross sections at 8 TeV)." )
            filename = f"xsecC1C1{int(sqrts)}.txt"
            #filename = "xsecN2N1p%d.txt" % sqrts
            order = NLL
            pb = False
            isEWK = True
        if pid1 in [ -1000024, -1000037 ] and pid2 in [ 1000024, 1000037 ]: # Charginos pair productions
            filename = f"xsecC1C1{int(sqrts)}.txt"
            order = NLL #3
            pb = False
            isEWK = True
        if filename == None:
            self.info ( f"could not identify filename for xsecs for {pid1,pid2}" )
            # self.info ( "seems like we dont have ref xsecs for the pids %d/%d?" % ( pid1, pid2 ) )
            return None, None, None
            # sys.exit()
        if ewk == "hino":
            filename = filename.replace(".txt","hino.txt" )
        if isEWK:
            if comment == "":
                comment = f"({ewk})"
        path = os.path.join ( self.shareDir, filename )
        self.info ( f"[refxsecComputer] will query {filename}" )
        if not os.path.exists ( path ):
            self.error ( f"{path} missing" )
            sys.exit(-1)
        xsecs = self.getXSecsFrom ( path, pb, columns )
        return xsecs,order,comment

if __name__ == "__main__":
    import argparse
    argparser = argparse.ArgumentParser( description = "compute xsecs by looking up reference cross sections" )
    argparser.add_argument ( "-f", "--inputfile",
            help="slha file(s) [./simplyGluino.slha]", nargs="*",
            type=str, default="./simplyGluino.slha" )
    argparser.add_argument ( '-s', '--sqrts',
            help='center-of-mass energies [8 13]',
            type=int, nargs="*", default=None )
    argparser.add_argument ( "-i", "--ignore_pids",
            help="ignore pids", type=str, default=None )
    args = argparser.parse_args()
    sqrts = args.sqrts
    if sqrts == None:
        sqrts = [ 8, 13 ]
    tool = RefXSecComputer( verbose = 2, allowN1N1Prod = False, walkerid = 0 )
    slhapaths = args.inputfile
    ssmultipliers = { (1000021,1000021):2. }
    ssmultipliers = None
    ignores= args.ignore_pids
    ignores= eval ( args.ignore_pids )
    for slhapath in slhapaths:
        tool.computeForOneFile ( sqrtses=sqrts, inputFile = slhapath, tofile=True,
                      ssmultipliers = ssmultipliers, ignore_pids = ignores )
