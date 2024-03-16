#!/usr/bin/env python3

import pickle, os, sys, subprocess, time, glob, colorama, math, scipy
import numpy as np
sys.path.insert(0,"../../")
from protomodels.csetup import setup
setup()
from ptools import hiscoreTools
from builder.manipulator import Manipulator
from tester.predictor import Predictor
from tester.combiner import Combiner
from builder import protomodel
from builder.protomodel import ProtoModel
from smodels.base.physicsUnits import fb, TeV
from smodels.matching.theoryPrediction import TheoryPrediction
from smodels.base import runtime
from smodels_utils.plotting import rulerPlotter, decayPlotter
from smodels_utils.helper.bibtexTools import BibtexWriter
from ptools.sparticleNames import SParticleNames
from smodels.base.smodelsLogging import logger
logger.setLevel("ERROR")
from os import PathLike
from colorama import Fore
from typing import Union, Dict, TextIO, List
from ptools.helpers import computeP, computeZFromP
from smodels_utils.helper.prettyDescriptions import prettyTexAnalysisName
from builder.loggerbase import LoggerBase
        
namer = SParticleNames ( susy = False )
bibtex = BibtexWriter()

# runtime._experimental = True

class HiscorePlotter ( LoggerBase ):
    def __init__ ( self ):
        super ( HiscorePlotter, self ).__init__ ( 0 )
        self.url = "https://smodels.github.io/"

    def gitCommit ( self, dest, upload, wanted : bool ):
        """ if wanted, then git commit and git push to smodels.githuib.io
        :param dest: e.g. "/scratch-cbe/users/w.w/git/smodels.github.io/protomodels/latest/"
        :param upload: e.g. "latest"
        :param wanted: were we even asked to git-commit?
        """
        if not wanted:
            return False
        comment = f"automated update by plotHiscore.py to {upload}:\n    K={self.protomodel.K:.3f} Z={self.protomodel.Z:.2f} walkerid={self.protomodel.walkerid}"
        destdir = dest
        destdir = destdir.replace(upload,"")
        destdir = destdir.replace("//","")
        cmd = f"cd {destdir}; git pull; git add '{upload}'; git commit -m '{comment}'; git push"
        self.pprint ( f"exec: {cmd}" )
        out = subprocess.getoutput ( cmd )
        if out != "":
            self.pprint ( f"{out}" )
        return True

    def discussPredictions ( self ):
        print ( "How the K comes about. Best combo:" )
        combo = self.protomodel.bestCombo
        self.pprint ( "best combination:" )
        self.pprint ( "=================" )
        for pred in combo:
            self.pprint ( f" - {pred.expResult.globalInfo.id}:{','.join ( set ( map ( str, pred.txnames ) ) )}" )

    def getExtremeSSMs ( self, ssm : Dict, largest : bool , nm : int = 7 ):
        """ get latex code describing the most extreme signal strength multipliers.
        :param largest: if true, describe the largest, else the smallest.
        :param nm: number of ssms to describe. -1 means "all"
        :returns: latex code
        """
        if len(ssm) == 0:
            ssm = { 0: "\\mathrm{none}" }
        keys = list ( ssm.keys() )
        keys.sort( reverse=largest )
        extreme = "smallest s"
        if largest:
            extreme = "largest s"
        if len(ssm)<8:
            extreme = "S"
        if nm == -1:
            extreme = "S"
            nm = len(keys)
        s = ""
        for ctr,k in enumerate(keys[:nm]):
            if k == 0 and ssm[k] == "\\mathrm{none}":
                s += "\\mathrm{none}  "
            else:
                s += f"{ssm[k]}={k}; "
        if len(s)>2:
            s = s[:-2]
        ret = f"{extreme}ignal strength multipliers: ${s}$"
        return ret

    def hasSignals ( self ):
        """ are there signals stored in the theory predictions? """
        for tp in self.protomodel.bestCombo:
            if hasattr ( tp.dataset, "dataInfo" ) and hasattr(tp.dataset.dataInfo,"sigN" ):
                return "did", True
            if hasattr ( tp.dataset, "txnameList" ):
                for txn in tp.dataset.txnameList:
                    if hasattr ( txn, "sigmaN" ):
                        return "did", True
        return "did not", False

    def significanceOfTP ( self, tp : TheoryPrediction ) -> Union[None,float]:
        """ compute the significance of one theory prediction 

        :param tp: the theory prediction to compute the significance for
        :returns: significance of theory prediction, None if unsuccessful
        """
        dtype = tp.dataType()
        if dtype == "efficiencyMap":
            dI = tp.dataset.dataInfo
            eBG = dI.expectedBG
            bgErr = dI.bgError
            p = computeP ( dI.observedN, eBG, bgErr )
            Z = computeZFromP ( p )
            return Z
        if dtype == "upperLimit":
            oUL = tp.getUpperLimit ( expected = False ).asNumber(fb)
            eUL = tp.getUpperLimit ( expected = True )
            if type(eUL)==type(None):
                return None ## cannot compute a significance
            eUL = eUL.asNumber(fb)
            sigma_exp = eUL / 1.96 # the expected scale, sigma
            Z = ( oUL - eUL ) / sigma_exp
            return Z
        if dtype == "combined":
            llhd = tp.likelihood( expected=False, return_nll=True )
            l0 = tp.lsm ( return_nll = True )
            chi2 = 2 * ( l0 - llhd )
            p = 1. - scipy.stats.chi2.cdf ( chi2, df=1 )
            Z = computeZFromP ( p )
            return Z
        return None

    def oneEntryHtmlRawNumbers ( self, tp : TheoryPrediction, f : TextIO,
           isInBestCombo : bool, Zvalue : Union[None,float] ) -> str:
        """ write one entry in rawnumbers.html into <f>
        :param tp: the theory prediction to be written the entry for.
        :param f: the file handle for rawnumbers.html
        :param isInBestCombo: indicate whether the theory prediction
        is part of the best combination. we will hilight these.

        :returns: 'analysisid:datatype' as a string
        """
        self.debug ( f"adding {tp.dataset.globalInfo.id} Z={Zvalue:.2f} to rawnumbers.html" )
        didordidnot,hassigs = self.hasSignals ( )
        anaId = tp.analysisId()
        idAndUrl = self.anaNameAndUrl ( tp )
        dtype = tp.dataType()
        ltopos = list(set(map(str,tp.txnames ) ))
        ltopos.sort()
        topos = ", ".join ( ltopos )
        S = "?"
        dt = { "upperLimit": "ul", "efficiencyMap": "em", "combined": "comb" }
        f.write ( f"<tr><td>{idAndUrl}</td><td>{dt[dtype]}</td> " )
        Z = self.significanceOfTP ( tp )
        if dtype == "efficiencyMap":
            dI = tp.dataset.dataInfo
            did = dI.dataId # .replace("_","\_")
            maxLen=9
            maxLen=18
            if len(did)>maxLen:
                did=did[:maxLen-3]+" ..."
            eBG = dI.expectedBG
            if eBG == int(eBG):
                eBG=int(eBG)
            bgErr = dI.bgError
            if bgErr == int(bgErr):
                bgErr=int(bgErr)
            S = "N/A"
            if Z != None:
                S = f"{Z:.1f} &sigma;" 
            pids = self.combiner.getAllPidsOfTheoryPred ( tp )
            obsN = dI.observedN
            if ( obsN - int(obsN) ) < 1e-6:
                obsN=int(obsN)
            particles = namer.htmlName ( pids, addSign = False,
                                          addBrackets = False )
            sobsN = str(obsN)
            if type(obsN) == float:
                sobsN = f"{obsN:.2f}" 
            f.write ( f'<td>{did}</td><td>{topos}</td><td>{sobsN}</td><td>{eBG:.2f} +/- {bgErr:.2f}</td><td style="text-align:right">{S}</td><td style="text-align:right">{particles}</td>' )
            if hassigs:
                sig = "-"
                if hasattr ( dI, "sigN" ):
                    sig = f"{dI.sigN}" 
                    if type(dI.sigN) in [ float, np.float64 ]:
                        sig = f"{dI.sigN:.2f}"
                f.write ( f'<td style="text-align:right">{sig}</td>' )
        if dtype == "upperLimit":
            S = "?"
            llhd = tp.likelihood( expected=False )
            oUL = tp.getUpperLimit ( expected = False ).asNumber(fb)
            eUL = tp.getUpperLimit ( expected = True )
            seUL = "N/A"
            S = "N/A"
            if Z != None:
                S = f"{Z:.1g} &sigma;"
            if type(eUL)!=type(None):
                eUL = eUL.asNumber(fb)
                seUL = f"{eUL:.1g} fb"
            pids = self.combiner.getAllPidsOfTheoryPred ( tp )
            particles = namer.htmlName ( pids, addSign = False, addBrackets = False )
            f.write ( f'<td>-</td><td>{topos}</td>' ) 
            f.write ( f'<td> {oUL:.1g} </td>' )
            f.write ( f'<td> {seUL} fb</td>' ) 
            f.write ( f'<td style="text-align:right">{S}</td>' )
            f.write ( f'<td style="text-align:right">{particles}</td>'  )
            if hassigs:
                sig = "-"
                for txn in tp.txnames:
                # for txn in tp.dataset.txnameList:
                    if hasattr ( txn, "sigmaN" ):
                        sig = f"{txn.sigmaN:.2f} fb"
                f.write ( f'<td style="text-align:right">{sig}</td>' )
        if dtype == "combined":
            S = "?"
            llhd = tp.likelihood( expected=False, return_nll=True )
            eUL = tp.getUpperLimit ( expected = True ).asNumber(fb)
            oUL = tp.getUpperLimit ( expected = False ).asNumber(fb)
            if Z != None:
                S = f"{Z:.1f} &sigma;"
            pids = self.combiner.getAllPidsOfTheoryPred ( tp )
            particles = namer.htmlName ( pids, addSign = False, addBrackets = False )
            f.write ( f'<td>-</td><td>{topos}</td>' ) 
            f.write ( f'<td> {oUL:.1g} fb </td><td> {eUL:.1g} fb</td>' )
            f.write ( f'<td style="text-align:right">{S}</td>' ) 
            f.write ( f'<td style="text-align:right">{particles}</td>')
            if hassigs:
                sig = "-"
                for txn in tp.txnames:
                # for txn in tp.dataset.txnameList:
                    if hasattr ( txn, "sigmaN" ):
                        sig = f"{txn.sigmaN:.2f} fb"
                f.write ( f'<td style="text-align:right">{sig}</td>' )
        f.write ( '</tr>\n' )
        return f"{anaId}:{dtype}"

    def createTpDictionaryWithZs ( self, 
            predictions : List[TheoryPrediction] ) -> Dict:
        """ create a dictionary with significances as keys,
        and predictions as values """
        tpAndZ = {} # dict of theory predictions with significances as keys
        for tp in predictions:
            Z = self.significanceOfTP ( tp )
            if Z in [ None ] or np.isnan ( Z ) or np.isinf(Z):
                Z = -100
            if np.isneginf(Z):
                Z = -100
            while Z in tpAndZ:
                Z+=1e-10
            tpAndZ[Z]=tp
        return tpAndZ

    def writeRawNumbersHtml ( self ):
        """ write out the raw numbers of the excess, as html """
        f=open("rawnumbers.html","wt")
        f.write("<table>\n" )
        f.write("<tr><th>Analysis Name</th><th>Type</th><th>Dataset</th><th>Topos</th><th>Observed</th><th>Expected</th><th>Approx &sigma;</th><th>Particles</th>" )
        didordidnot,hassigs = self.hasSignals ( )
        self.pprint ( f"protomodel's database {didordidnot} have fake signals." )
        if hassigs:
            f.write("<th>Signal</th>" )
        f.write("\n</tr>\n" )
        hasListed = []
        tpAndZ = self.createTpDictionaryWithZs ( self.protomodel.bestCombo )
        Zvalues = list ( tpAndZ.keys() )
        # Zvalues.sort( key = lambda x: -x if type(x) in [ float, int ] else 100 )
        Zvalues.sort( reverse = True )
        for Zvalue in Zvalues:
            tp = tpAndZ[Zvalue]
            idfer = self.oneEntryHtmlRawNumbers ( tp, f, True, Zvalue )
            hasListed.append ( idfer )
        f.write("</table>\n" )
        f.write("<table style='color:grey'>\n" )
        f.write("<tr><th>Analysis Name</th><th>Type</th><th>Dataset</th><th>Topos</th><th>Observed</th><th>Expected</th><th>Approx &sigma;</th><th>Particles</th>" )
        tpAndZ = self.createTpDictionaryWithZs ( self.predictor.predictions )
        Zvalues = list ( tpAndZ.keys() )
        # Zvalues.sort( key = lambda x: -x if type(x) in [ float, int ] else 100 )
        Zvalues.sort( reverse = True )
        for Zvalue in Zvalues:
            #if Zvalue < -10:
            #    continue
            tp = tpAndZ[Zvalue]
            anaId = tp.analysisId()
            dtype = tp.dataType()
            idfer = f"{anaId}:{dtype}"
            if not idfer in hasListed:
                anaType = self.oneEntryHtmlRawNumbers ( tp, f, False, Zvalue )
        f.write("</table>\n" )
        f.close()

    def writeRawNumbersLatex ( self, usePrettyNames : bool = True ):
        """ write out the raw numbers of the excess, in latex,
        to rawnumbers.tex

        :param usePrettyNames: use the pretty names, not analysis ids
        """
        print ( "raw numbers of excess" )
        print ( "=====================" )
        f=open("rawnumbers.tex","wt")
        f.write("\\begin{tabular}{l|c|r|r|c|r|r}\n" )
        f.write("\\bf{Analysis Name} & \\bf{Dataset} & \\bf{Obs} & \\bf{Expected} & \\bf{Z} & \\bf{Particle} & \\bf{Signal} \\\\\n" )
        f.write("\\hline\n" )
        for tp in self.protomodel.bestCombo:
            self.oneEntryTexRawNumbers ( tp, f, usePrettyNames )
        f.write("\end{tabular}\n" )
        f.close()

    def oneEntryTexRawNumbers ( self, tp : TheoryPrediction, f : TextIO,
           usePrettyNames : bool ) -> str:
        """ write one tex entry for the given theory prediction
        :param tp: write the entry for that theory prediction
        :param usePrettyNames: use the pretty names, not analysis ids

        :returns: anaid:datatype
        """
        anaId = tp.analysisId()
        ananame = anaId
        if usePrettyNames:
            ananame = prettyTexAnalysisName ( None, anaid = anaId )
        dtype = tp.dataType()
        print ( f"[plotHiscore] item {anaId} ({dtype})" )
        dt = { "upperLimit": "ul", "efficiencyMap": "em" }
        ref = bibtex.query ( anaId )
        f.write ( f"{ananame}~\\cite{{ref}} & " )
        Z = self.significanceOfTP ( tp )
        if dtype == "efficiencyMap":
            dI = tp.dataset.dataInfo
            obsN = dI.observedN
            if ( obsN - int(obsN) ) < 1e-6:
                obsN=int(obsN)
            print ( f"  `- {dI.dataId}: observedN {obsN}, bg {dI.expectedBG} +/- {dI.bgError}" )
            did = dI.dataId.replace("_","\_")
            if len(did)>9:
                did=did[:6]+" ..."
            eBG = dI.expectedBG
            if eBG == int(eBG):
                eBG=int(eBG)
            bgErr = dI.bgError
            if bgErr == int(bgErr):
                bgErr=int(bgErr)
            S = "N/A"
            if Z != None:
                S = f"{Z:.1f} $\sigma$"
            # pids = tp.PIDs
            pids = self.combiner.getAllPidsOfTheoryPred ( tp )
            particles = namer.texName ( pids, addDollars=True, addSign = False,
                                          addBrackets = False )
            obs = dI.observedN
            if obs == 0.:
                obs = 0
            else:
                if abs ( obs - int(obs) ) / obs < 1e-6:
                    obs = int ( obs )
            sigN = tp.xsection.asNumber(fb) * tp.dataset.globalInfo.lumi.asNumber(1/fb)
            sigmapred=f"{sigN:.2f}"
            f.write ( f"{did} & {obs} & {eBG} $\\pm$ {bgErr} & {S} & {particles} & {sigmapred} \\\\ \n" )
        if dtype in [ "upperLimit", "combined" ]:
            S = "?"
            llhd = tp.likelihood ( expected=False, return_nll = True )
            eUL = tp.getUpperLimit ( expected = True ).asNumber(fb)
            oUL = tp.getUpperLimit ( expected = False ).asNumber(fb)
            sigma_exp = eUL / 1.96 # the expected scale, sigma
            if Z is not None:
                S = f"{Z:.1f} $\sigma$"
            pids = self.combiner.getAllPidsOfTheoryPred ( tp )
            particles = namer.texName ( pids, addDollars=True, addSign = False,
                                        addBrackets = False )
            sigmapred=f"{tp.xsection.asNumber(fb)} fb"
            print ( f"  `- observed {oUL:.2f}*fb, expected {eUL:.2f}*fb {Z:.1f} sigma" )
            f.write ( f" & {oUL:.1f} fb & {eUL:.1f} fb & {S} & {particles} & {sigmapred} \\\\ \n" )
        return f"{anaId}:{dtype}"

    def findXSecOfPids ( self, xsecs, pids ):
        """ find the cross sections for pids
        :returns: xsec, as unum object
        """
        for xsec in xsecs:
            sqrts = xsec.info.sqrts.asNumber(TeV)
            if sqrts < 10:
                continue
            order = xsec.info.order
            if order > 0:
                continue
            if pids == xsec.pid:
                return xsec
        return 0.*fb

    def writeTex ( self, keep_tex : bool ):
        """ write the comment about ss multipliers and particle contributions, in tex.
        Creates texdoc.png.
        :param keep_tex: keep tex source of texdoc.png
        """
        cpids = {}
        frozen = self.protomodel.frozenParticles()
        xsecs = self.protomodel.getXsecs()[0]
        ssms = self.getUnfrozenSSMs ( frozen, includeOnes=True )
        for pids,v in ssms.items():
            xsec = self.findXSecOfPids ( xsecs, pids )
            if xsec < 0.001 * fb: ## only for xsecs we care about
                continue
            sv = f"{v:.2g}"
            if v < .1:
                sv = f"{v:.1g}"
            if not sv in cpids:
                cpids[sv]=[]
            cpids[sv].append ( pids )

        ssm = {}
        for v,pids in cpids.items():
            sp = []
            for pairs in pids:
                pname = namer.texName ( pairs, addSign = True, addBrackets=True )
                sp.append ( pname )
            ssm[v] = ", ".join ( sp )

        particleContributionList = ""
        if hasattr ( self.protomodel, "particleContributions" ):
            print ( "[plotHiscore] contributions-by-particle are defined" )
            #particleContributionList+="\\\\Contributions by particles: $"
            particleContributionList+="\\\\"
            particleContributionList+="Contributions by particles: $"
            totalcont = 0. ## to normalize contributions
            for k,v in self.protomodel.particleContributions.items():
                totalcont += (self.protomodel.K - v)
            tok = {}
            for k,v in self.protomodel.particleContributions.items():
                if v in tok.keys():
                    v+=1e-6
                perc = 100.
                if totalcont != 0.:
                    perc = round(100.*(self.protomodel.K - v)/totalcont )
                tok[v] = f"{namer.texName(k)}: K_\mathrm{{without}}={v:.2f} ({perc}%)"
            keys = list ( tok.keys() )
            keys.sort()
            for v in keys:
                particleContributionList+= tok[v] + ", "
            if len(keys)>0:
                particleContributionList = particleContributionList[:-2]
            #particleContributionList+= ", ".join ( tok )
            particleContributionList+="$"
        else:
            print ( "[plotHiscore] protomodel has no ``particleContributions'' defined." )

        from plotting import tex2png
        src = self.getExtremeSSMs ( ssm, largest=True, nm = 7 )
        src += "\\\\"
        nsmallest = min( 7, len(ssm)-7 )
        if nsmallest>0:
            src += self.getExtremeSSMs ( ssm, largest=False, nm = nsmallest )
        src += particleContributionList
        if keep_tex:
            with open("texdoc.tex","wt") as f:
                f.write ( src+"\n" )
                f.close()
            print ( f"[plotHiscore] wrote {os.getcwd()}/texdoc.tex" )
        try:
            p = tex2png.Latex ( src, 600 ).write()
            f = open ( "texdoc.png", "wb" )
            f.write ( p[0] )
            f.close()
        except Exception as e:
            print ( f"[plotHiscore] Exception when latexing: {e}" )
        return src

    def anaNameAndUrl ( self, ana : Union[str,TheoryPrediction], 
            forPdf : bool = False ) -> str:
        """ given analysis, return analysis name and URL,
        as html code or pdf hyperref
        :param ana: ExpRes or TheoryPred or str object
        :param forPdf: if True, create for Pdf hyperref, else for html

        :returns: html or hyperref string
        """
        if forPdf:
            ## FIXME implement!
            return ana.analysisId()
        if type(ana)==str:
            for i in self.protomodel.bestCombo:
                url = ""
                if i.analysisId() == ana and hasattr ( i.dataset.globalInfo, "url" ):
                    url = i.dataset.globalInfo.url
                    break
            return f"<a href={url}>{ana}</a>"
        if type(ana)==TheoryPrediction:
            if not hasattr ( ana.dataset.globalInfo, "url" ):
                return ( ana.analysisId() )
            return f"<a href={ana.dataset.globalInfo.url}>{ana.analysisId()}</a>"

    def getPrettyName ( self, rv ):
        if hasattr ( rv.dataset.globalInfo, "prettyName" ):
            prettyName = rv.dataset.globalInfo.prettyName
            return prettyName
        anaId = rv.analysisId()
        pnames = { "CMS-SUS-19-006": "0L + jets, MHT",
                   "CMS-SUS-16-033": "0L + jets + Etmiss (using MHT)",
                   "CMS-SUS-16-036": "0L + jets + Etmiss (using MT2)",
                   "CMS-SUS-16-049": "All hadronic stop"
        }
        if anaId in pnames:
            return pnames[anaId]
        return anaId

    def writeRValuesTex ( self, rvalues, usePrettyNames = True ):
        """ write out the leading rvalues of the critic, in latex
        :param usePrettyNames: use the pretty names, not analysis ids
        """
        g=open("rvalues.tex","wt")
        g.write ( "\\begin{tabular}{l|c|c|c|c|c}\n" )
        g.write ( "\\bf{Analysis Name} & \\bf{Production} & $\sigma_{XX}$ (fb) & $\sigma^\mathrm{UL}_\mathrm{obs}$ (fb) & $\sigma^\mathrm{UL}_\mathrm{exp}$ (fb) & $r$ \\\\\n" )
        #g.write ( "\\begin{tabular}{l|c|r|r}\n" )
        #g.write ( "\\bf{Analysis Name} & \\bf{Topo} & $r_{\mathrm{obs}}$ & $r_{\mathrm{exp}}$ \\\\\n" )
        g.write ( "\\hline\n" )
        for rv in rvalues[:5]:
            srv="N/A"
            if type(rv['rexp']) in [ float, np.float64, np.float32 ]:
                srv="%.2f" % rv['rexp']
            else:
                srv=str(rv['rexp'])
            anaId = rv['tp'].analysisId()
            prettyName = getPrettyName( rv['tp'] )
            prettyName = prettyTexAnalysisName ( prettyName, anaid = anaId )
            ref = bibtex.query ( anaId )
            txnames = ",".join ( map(str,rv['tp'].txnames) )
            allpids = rv['tp'].PIDs
            pids=[]
            for p in allpids:
                tmp=[]
                for b in p:
                    if type(b[0])==int:
                        tmp.append(b[0])
                    else:
                        tmp.append(b[0][0])
                pids.append(tuple(tmp))
            prod = []
            for p in pids:
                tmp = namer.texName ( p, addDollars=True, addSign = False,
                                       addBrackets = True )
                prod.append (tmp)
            prod = "; ".join(prod)
            sigmapred = "20.13"
            sigmapred = rv['tp'].xsection.asNumber(fb)
            sigmaexp = "--"
            if type(rv['tp'].getUpperLimit ( expected = True )) != type(None):
                sigmaexp = "%.2f" % rv['tp'].getUpperLimit ( expected=True ).asNumber(fb)
            sigmaobs = rv['tp'].getUpperLimit().asNumber(fb)
            g.write ( f"{prettyName}~\\cite{{{ref}}} & {prod} & {sigmapred:.2f} & {sigmaobs:.2f} & {sigmaexp} & {rv['obs']:.2f}\\\\\n" )
        g.write ( "\\end{tabular}\n" )
        g.close()

    def getDatabaseVersion ( self, dbpath = "default.pcl" ):
        dbver = "???"
        if hasattr ( self.protomodel, "dbversion" ):
            dbver = self.protomodel.dbversion
            if not "???" in dbver:
                return dbver
        if os.path.exists ( dbpath ):
            ## try to get db version from db file
            from smodels.experiment.databaseObj import Database
            db = Database ( dbpath )
            dbver = db.databaseVersion
        return dbver

    def writeIndexTex ( self, texdoc ):
        """ write the index.tex file
        :param texdoc: the source that goes into texdoc.png
        """
        ssm = []
        for k,v in self.protomodel.ssmultipliers.items():
            if abs(v-1.)<1e-3:
                continue
            ssm.append ( f"{namer.texName(k,addSign=True)}: {v:.2f}" )
        f=open("index.tex","w")
        f.write ( f"Our current winner has a score of \\K={self.protomodel.K:.2f}, ")
        strategy = "aggressive"
        dbver = self.getDatabaseVersion ( )
        dotlessv = dbver.replace(".","")
        f.write ( f" it was produced with database {{\\tt v{dotlessv}}}, combination strategy {{\\tt {strategy}}} walker {self.protomodel.walkerid} in step {self.protomodel.step}." )
        f.write ( "\n" )
        if hasattr ( protomodel, "tpList" ):
            rvalues=protomodel.tpList
            rvalues.sort(key=lambda x: x['robs'],reverse=True )
            writeRValuesTex ( rvalues )
            #writeRValuesTexOld ( rvalues )
        else:
            print ( "[plotHiscore] protomodel has no r values!" )

        if hasattr ( protomodel, "analysisContributions" ):
            print ( "[plotHiscore] contributions-per-analysis are defined" )
            # f.write ( "Contributions per analysis:\n\\begin{itemize}[noitemsep,nolistsep]\n" )
            f.write ( "The contributions per analysis are given in Tab.~\\ref{tab:analysiscontributions}.\n" )
            f.write ( "\\begin{table}\n" )
            f.write ( "\\begin{center}\n" )
            f.write ( "\\begin{tabular}{l|c|c}\n" )
            f.write ( "\\bf{Analysis Name} & \\bf{\\K(without)} & \\bf{Contribution} \\\\\n" )
            f.write ( "\\hline\n" )
            conts = []
            dKtot = 0.
            Ktot = protomodel.K
            contributions = protomodel.analysisContributions
            for k,v in contributions.items():
                conts.append ( ( v, k ) )
                dKtot += Ktot - v
            # dKtot = dKtot # * protomodel.K
            Ks, Kwo = {}, {}
            """
            for k,v in contributions:
                Ks[k] = ( protomodel.K - dKtot * v)
                Kwo[k] = protomodel.K
            """
            conts.sort( reverse=True )
            for v,k in conts:
                Kwithout= contributions[k]
                cont = ( Ktot - Kwithout ) / dKtot
                f.write (f"{k} & {Kwithout:.2f} & {int(round(100.*cont))}% \\\\ \n")
            f.write ( "\\end{tabular}\n" )
            f.write ( "\\end{center}\n" )
            f.write ( "\\caption{Contributions to the test statistic \\K. \\K(without) denotes the \\K value obtained in absence of the particular analysis.}\n" )
            f.write ( "\\label{tab:analysiscontributions}\n" )
            f.write ( "\\end{table}\n" )
        else:
            print ( "[plotHiscore] contributions are not defined" )

        height = 32*int((len(ssm)+3)/4)
        if ssm == []:
            height = 32
        if hasattr ( protomodel, "particleContributions" ):
            height += 32
        contrs = texdoc.replace(":"," are " ).replace("S","The s").replace(";",", " )
        contrs = contrs.replace( "\\\\\\\\Contr", "; the contr" )
        f.write ( contrs + "\n" )
        f.close()
        print ( "[plotHiscore] Wrote index.tex" )

    def getUnfrozenSSMs ( self, frozen, includeOnes=False, dropLSPLSP=True ):
        """ of all SSMs, leave out the ones with frozen particles
        :param frozen: list of pids of frozen particles
        :param includeOnes: if False, then also filter out values close to unity
        :param dropLSPLSP: if True, dont include LSP,LSP production
        :returns: dictionary of SSMs without frozen particles
        """
        # ssms = protomodel.ssmultipliers
        ma = Manipulator ( self.protomodel )
        ssms = ma.simplifySSMs ( threshold = .01 * fb )
        D={}
        for pids,v in ssms.items():
            hasFrozenParticle = False
            if dropLSPLSP:
                if pids == ( ProtoModel.LSP, ProtoModel.LSP ):
                    continue

            for pid in pids:
                if pid in frozen or -pid in frozen:
                    hasFrozenParticle = True
            if hasFrozenParticle:
                continue
            if not includeOnes and abs(v-1.)<1e-2:
                continue
            D[pids]=v
        return D

    def addLikelihoodPlots ( self, handle : TextIO ):
        """ add links to the likelihood plots """
        files = glob.glob ( "llhd*.png" )
        if len(files)==0:
            return
        handle.write ( "llhd plots:")
        for aFile in files:
            dt = int ( time.time() - 1709120000 )
            fname = aFile.replace(".png","").replace("llhd","")
            line = f" <a href={aFile}?{dt}>{fname}</a>"
            handle.write ( line )
        handle.write ( "\n" )
            

    def addSPlots ( self, handle : TextIO ):
        """ add links to the teststatistic plots """
        files = glob.glob ( "M*.png" )
        if len(files)==0:
            return
        handle.write ( "S plots:")
        for aFile in files:
            dt = int ( time.time() - 1709120000 )
            fname = aFile.replace(".png","").replace("M","")
            line = f" <a href={aFile}?{dt}>{fname}</a>"
            handle.write ( line )
        handle.write ( "\n" )

    def writeIndexHtml ( self ):
        """ write the index.html file, see e.g.
            https://smodels.github.io/protomodels/
        """
        ssm = []
        frozen = self.protomodel.frozenParticles()
        ssms = self.getUnfrozenSSMs ( frozen, includeOnes=True )
        for k,v in ssms.items():
            ssm.append ( f"{namer.htmlName(k,addSign=True) }: {v:.2g}" )
        f=open("index.html","w")
        f.write ( "<html>\n" )
        f.write ( "<body>\n" )
        f.write ( "<center>\n" )
        f.write ( f"<table><td><h1>" )
        f.write ( f"Current best protomodel: <i>K</i>={self.protomodel.K:.2f}" )
        f.write ( f", <i>Z</i>={self.protomodel.Z:.2f}" )
        f.write ( f"</h1><td>" )
        f.write ( f"<img height=60px src={self.url}/protomodels/logos/protomodel_lego.png>" )
        f.write ( "</table>\n" )
        f.write ( "</center>\n" )
        dbver = self.getDatabaseVersion (  )
        strategy = "aggressive"
        dotlessv = dbver.replace(".","")
        dt = int ( time.time() - 1593000000 )
        f.write ( "<b><a href=./hiscore.slha>ProtoModel</a> <a href=./pmodel.dict>(dict)</a> " )
        f.write ( f"produced with <a href={self.url}/docs/Validation{dotlessv}>database v{dbver}</a>" )
        f.write ( f", combination strategy <a href=./matrix.png>{strategy}</a> in walker {self.protomodel.walkerid} step {self.protomodel.step}.</b> " )
        if hasattr ( self.protomodel, "particleContributions" ):
            f.write ( f"<i>K</i> plots for: <a href=./M1000022.png?{dt}>{namer.htmlName(1000022)}</a>" )
            for k,v in self.protomodel.particleContributions.items():
                f.write ( ", " )
                f.write ( f"<a href=./M{k}.png?{dt}>{namer.htmlName(k)}</a>" )
            f.write ( ". HPD plots for: " )
            first = True
            for k,v in self.protomodel.particleContributions.items():
                if not first:
                    f.write ( ", " )
                f.write ( f"<a href=./llhd{k}.png?{dt}>{namer.htmlName(k)}</a>" )
                first = False
        # fixme replace with some autodetection mechanism
        # take out all frozen ssm plots
        self.addLikelihoodPlots ( f )
        self.addSPlots ( f )
        """
        ossms = { (-1000006,1000006), (1000021,1000021), (-2000006,2000006) }
        for fname in glob.glob("ssm_*_*.png" ):
            pids = fname.replace("ssm_","").replace(".png","")
            pids = tuple ( map ( str, pids.split("_") ) )
            ossms.add ( pids )
        frozen = self.protomodel.frozenParticles()
        ssms = set()
        for pids in ossms:
            hasFrozenPid=False
            for pid in pids:
                if pid in frozen: #  or -pid in frozen:
                    hasFrozenPid = True
                    break
            if not hasFrozenPid:
                ssms.add ( pids )


        if len(ssms)>0:
        f.write ( " SSM plots for: " )
        first = True
        for pids in ssms:
            if not first:
                f.write ( ", " )
            f.write ( f"<a href=./ssm_{pids[0]}_{pids[1]}.png?{dt}>({pids[0]},{pids[1]})</a>"  )
            first = False
        """
        f.write ( "<br>\n" )
        f.write ( "<table width=80%>\n<tr><td>\n" )
        if hasattr ( self.protomodel, "tpList" ):
            rvalues=self.protomodel.tpList
            rvalues.sort(key=lambda x: x['robs'],reverse=True )
            f.write ( f"<br><b>{len(rvalues)} predictions available. Highest r values are:</b><br><ul>\n" )
            for rv in rvalues[:4]:
                srv="N/A"
                if type(rv['rexp']) in [ float, np.float64, np.float32 ]:
                    srv= f"{rv['rexp']:.2f}"
                elif type(rv['rexp']) != type(None):
                    srv=str(rv['rexp'])
                f.write ( f"<li>{self.anaNameAndUrl ( rv['tp'] )}:{rv['tp'].dataType(short=True)}:{','.join ( set (map(str,rv['tp'].txnames) ) )} r={rv['robs']:.2f}, r<sub>exp</sub>={srv}<br>\n" )
            f.write("</ul>\n")
        else:
            print ( "[plotHiscore] protomodel has no r values!" )

        if hasattr ( self.protomodel, "analysisContributions" ):
            print ( "[plotHiscore] contributions-per-analysis are defined" )
            f.write ( "<td><br><b>contributions to K:</b><br>\n<ul>\n" )
            conts = []
            Ktot = self.protomodel.K
            dKtot = 0.
            contributions = self.protomodel.analysisContributions
            for k,v in contributions.items():
                conts.append ( ( v, k ) )
                dKtot += Ktot - v
            conts.sort( reverse=True )
            for v,k in conts:
                Kwithout= contributions[k]
                cont = ( Ktot - Kwithout ) / dKtot
                nameAndUrl = self.anaNameAndUrl ( k )
                kv = str(v)
                if type(v) in [ float, np.float64 ]:
                    kv = f"{v:.2f} ({int(round(100.*cont))}%)"
                f.write ( f"<li> {nameAndUrl}: {kv}\n" )
            # f.write ( "</table>\n" )
        else:
            print ( "[plotHiscore] analysis-contributions are not defined" )

        height = 32*int((len(ssm)+3)/4)
        if ssm == []:
            height = 32
        if hasattr ( protomodel, "particleContributions" ):
            height += 32
        t0 = int(time.time())
        f.write ( f"<td><img width=600px src=./texdoc.png?{t0}>\n" )
        f.write ( f"<br><font size=-1>Last updated: {time.asctime()}</font>\n" )
        f.write ( "</table>" )
        f.write ( '<table style="width:80%">\n' )
        f.write ( "<td width=45%>" )
        f.write ( f"<img height=580px src=./ruler.png?{t0}>" )
        f.write ( "<td width=55%>" )
        f.write ( f"<img height=220px src=./decays.png?{t0}>\n" )
        f.write ( f'<font size=-3><iframe type="text/html" height="300px" width="100%" frameborder="0" src="./rawnumbers.html?{t0}"></iframe></font>\n' )
        f.write ( "</table>\n" )
        f.write ( "</body>\n" )
        f.write ( "</html>\n" )
        f.close()
        print ( "[plotHiscore] Wrote index.html" )

    def copyFilesToGithub( self ):
        files = [ "hiscore.slha", "index.html", "decays.png",
                  "ruler.png", "texdoc.png", "pmodel.dict", "rawnumbers.html" ]
        # files += [ "matrix.png" ]
        for f in files:
            if not os.path.exists ( f ):
                continue
            O = subprocess.getoutput ( f"cp {f} ../../smodels.github.io/protomodels/" )
            if len(O)>0:
                print ( "[plotHiscore.py] when copying files: %s" % O )

    def getPIDsOfTPred ( self, tpred, ret, integrateDataType=True, integrateSRs=True ):
        """ get the list of PIDs that the theory prediction should be assigned to
        :param tpred: theory prediction
        :param ret: results of a previous run of this function, so we can add iteratively
        :param integrateDataType: if False, then use anaid:dtype (eg CMS-SUS-19-006:ul) as values
        :param integrateSRs: if False, then use anaid:SR (eg CMS-SUS-19-006:SRC) as values.
                             takes precedence over integrateDataType
        :returns: dictionary with pid as key and sets of ana ids as value
        """
        LSP = 1000022
        name = tpred.analysisId()
        if not integrateSRs:
            SR = tpred.dataId()
            name = name + ":" + str(SR)
        elif not integrateDataType:
            dType = "em"
            if tpred.dataId() in [ "None", None ]:
                dType = "ul"
            name = name + ":" + dType
        pids = self.combiner.getAllPidsOfTheoryPred ( tpred )
        for pid in pids:
            if pid == LSP:
                continue
            if not pid in ret:
                ret[pid]=set()
            ret[pid].add ( name )
        return ret

    def plotRuler( self, verbosity : str, horizontal : bool ):
        """ plot the ruler plot, given protomodel.
        :param verbosity: one of: error, warning, info, debug
        :param horizontal: plot horizontal ruler plot, if true. Else, vertical.
        """
        fname = "ruler.png"
        if horizontal:
            fname = "horizontal.png"
        print ( f"[plotHiscore] now draw {fname}" )
        resultsForPIDs = {}
        for tpred in self.protomodel.bestCombo:
            resultsForPIDs = self.getPIDsOfTPred ( tpred, resultsForPIDs )
        resultsFor = {}
        for pid,values in resultsForPIDs.items():
            if pid in self.protomodel.masses:
                resultsFor[ self.protomodel.masses[pid] ] = values
            else:
                print ( "[plotHiscore] why is pid %s not in mass dict %s?" % ( pid, str(protomodel.masses) ) )

        if verbosity == "debug":
            print ( f'[plotHiscore] ../smodels-utils/smodels_utils/plotting/rulerPlotter.py -o ruler.png --hasResultsFor "{str(resultsFor)}" {self.protomodel.currentSLHA}'  )

        plotter = rulerPlotter.RulerPlot ( self.protomodel.currentSLHA, fname,
                                           Range=(None, None), mergesquark = False,
                                           drawdecays = False,
                                           hasResultsFor = resultsFor,
                                           trim = True )
        if horizontal:

            plotter.drawHorizontal()
        else:
            plotter.drawVertical()

    def plotDecays ( self, verbosity : str, outfile : str = "decays.png" ):
        try:
            import pygraphviz
        except ModuleNotFoundError as e:
            print ( f"[plotHiscore] skipping decays, no pygraphviz found!" )
            return
        verbosity = verbosity.lower().strip()
        print ( f"[plotHiscore] now draw {outfile}" )
        options = { "tex": True, "color": True, "dot": True, "squarks": True,
                    "weakinos": True, "sleptons": True, "neato": True,
                    "separatecharm": True,
                    "integratesquarks": False, "leptons": True }
        options["rmin"] = 0.005
        ## FIXME add cross sections.
        if verbosity in [ "debug", "info" ]:
            soptions = ""
            for k,v in options.items():
                if v==True:
                    soptions += f"--{k} "
            ma = Manipulator ( self.protomodel )
            ssms = ma.simplifySSMs()
            # soptions+=' --ssmultipliers "%s"' % ssms
            print ( f"{Fore.GREEN}../smodels-utils/smodels_utils/plotting/decayPlotter.py -f {self.protomodel.currentSLHA} -o {outfile} {soptions}{Fore.RESET}" )
        decayPlotter.draw ( self.protomodel.currentSLHA, outfile, options,
                            verbosity = verbosity,
                            ssmultipliers = self.protomodel.ssmultipliers )

    def plot ( self, number, verbosity, hiscorefile, options, dbpath ):
        ## plot hiscore number "number"
        pm = hiscoreTools.obtainHiscore ( number, hiscorefile )
        self.protomodel = pm
        self.combiner = Combiner ( self.protomodel.walkerid )
        self.predictor = Predictor ( 0, dbpath, do_srcombine = True )

        protoslha = self.protomodel.createSLHAFile ()
        subprocess.getoutput ( f"cp {protoslha} hiscore.slha" )
        m = Manipulator ( self.protomodel )
        print ( "[plotHiscore] now write pmodel.dict" )
        m.writeDictFile()
        opts = [ "ruler", "decays", "predictions", "copy", "html" ]
        for i in opts:
            if not i in options:
                options[i]=True

        plotruler = options["ruler"]
        horizontal = False
        if "horizontal" in options and options["horizontal"]:
            horizontal = True
        if plotruler:
            self.plotRuler ( verbosity, horizontal )
        plotdecays = options["decays"]
        if plotdecays:
            self.plotDecays ( verbosity )

        if options["predictions"]:
            self.discussPredictions ( )
        if options["html"] or options["tex"]:
            texdoc = self.writeTex ( options["keep"] )
            if options["html"]:
                self.writeIndexHtml ( )
            if options["tex"]:
                self.writeIndexTex( texdoc )
        self.predictor.predict ( self.protomodel, keep_predictions = True )
        self.writeRawNumbersLatex ( )
        self.writeRawNumbersHtml ( )
        if options["keep"]:
            print ( f"[plotHiscore] keeping {self.protomodel.currentSLHA}" )
        else:
            self.protomodel.delCurrentSLHA()

    def compileTestText( self ):
        subprocess.getoutput ( "pdflatex test.tex" )

def runPlotting ( args ):
    if args.destinations:
        print ( "Upload destinations: " )
        print ( "      none: no upload" )
        print ( "       gpu: upload to GPU server, afs space." )
        print ( "            Result can be seen at http://smodels.github.io/protomodels/" )
        print ( "    github: upload to github git directory." )
        print ( f"            Result can be seen at {self.url}/protomodels" )
        print ( "interesting: upload to github git directory, 'interesting' folder." )
        print ( f"             Result can be seen at {self.url}/protomodels/interesting" )
        print ( "latest: upload to github git directory, 'latest' folder." )
        print ( f"             Result can be seen at {self.url}/protomodels/latest" )
        print ( "anomaly: upload to github git directory, 'anomaly' folder." )
        print ( f"             Result can be seen at {self.url}/protomodels/anomaly" )
        return
    upload = args.upload# .lower()
    if upload in [ "none", "None", "" ]:
        upload = None

    options = { "ruler": args.ruler, "decays": args.decays,
                "predictions": args.predictions, "html": args.html,
                "keep": args.keep, "tex": args.tex,
                "horizontal": args.horizontal }

    hiplt = HiscorePlotter()
    hiplt.plot ( args.number, args.verbosity, args.hiscorefile, options,
                   args.dbpath )
    if upload is None:
        return
    F = "decays.png ruler.png texdoc.png pmodel.dict hiscore.slha index.html rawnumbers.html"
    dest = ""
    destdir = f'{os.environ["HOME"]}/git'
    dest = f"{destdir}/smodels.github.io/protomodels/{upload}/"
    if upload == "github":
        dest = f"{destdir}/smodels.github.io/protomodels/"
    if upload in [ "interesting", "anomaly", "latest" ]:
        dest = f"{destdir}/smodels.github.io/protomodels/{upload}/"
    if "paper" in upload:
        dest = f"{destdir}/smodels.github.io/protomodels/{upload}/"

    if dest != "":
        if not os.path.exists ( dest ):
            print ( f"[plotHiscore] directory {dest} does not exist. Trying to mkdir it." )
            cmd = f"mkdir {dest}" ## -p ?
            a = subprocess.getoutput ( cmd )
            if a != "":
                print ( f"[plotHiscore] error when mkdir: {a}" )
        print ( f"[plotHiscore] copying to {dest}:", end="" )
        first=True
        for f in F.split(" "):
            if os.path.exists ( f ):
                if not first:
                    print ( ",", end="" )
                print ( " "+f, end="" )
                cmd = "cp"
                cmd = f"{cmd} {f} {dest}"
                a = subprocess.getoutput ( cmd )
                if a != "":
                    print ( "error: %s" % a )
                    sys.exit()
                first = False
        print ( "." )
        cleanUp = True ## clean up after moving
        if cleanUp:
            toClean = [ "decays.tex", "decays.log", "decays.aux", "decays.pdf",
                        "rawnumbers.tex", "decays.dot", "decays.svg" ]
            for f in toClean:
                if os.path.exists ( f ):
                    os.unlink ( f )
        r = hiplt.gitCommit( dest, upload, args.commit )
        if not r:
            destdir = dest
            destdir = destdir.replace(upload,"")
            destdir = destdir.replace("//","")
            print ( "[plotHiscore] done. now please do yourself: " )
            print ( f"cd {destdir}" )
            print ( f"git add {upload}" )
            print ( "git commit -m 'manual update'" )
            print ( "git push" )
        return

    print ( f"error, dont know what to do with upload sink '{upload}'" )


def main ():
    import argparse
    argparser = argparse.ArgumentParser(
            description='hiscore proto-model plotter')
    argparser.add_argument ( '-n', '--number',
            help='which hiscore to plot [0]',
            type=int, default=0 )
    argparser.add_argument ( '-f', '--hiscorefile',
            help='pickle file to draw from [<rundir>/hiscores.cache]',
            type=str, default="default"  )
    argparser.add_argument ( '-v', '--verbosity',
            help='verbosity -- debug, info, warn, err [info]',
            type=str, default="warn" )
    argparser.add_argument ( '-H', '--html',
            help='produce index.html',
            action="store_true" )
    argparser.add_argument ( '-R', '--ruler',
            help='produce ruler plot',
            action="store_true" )
    argparser.add_argument ( '-D', '--decays',
            help='produce decays plot',
            action="store_true" )
    argparser.add_argument ( '-P', '--predictions',
            help='list all predictions',
            action="store_true" )
    argparser.add_argument ( '-T', '--tex',
            help='produce the latex version',
            action="store_true" )
    argparser.add_argument ( '-A', '--all',
            help='produce everything (equivalent to -H -R -D -P -T)',
            action="store_true" )
    argparser.add_argument ( '-t', '--test',
            help='produce test.pdf file',
            action="store_true" )
    argparser.add_argument ( '-k', '--keep',
            help='keep cruft files',
            action="store_true" )
    argparser.add_argument ( '--horizontal',
            help='horizontal, not vertical ruler plot?',
            action="store_true" )
    argparser.add_argument ( '-u', '--upload',
            help='upload to one of the following destinations: none, gpu, github, anomaly, latest, interesting [none]. run --destinations to learn more',
            type=str, default="" )
    argparser.add_argument ( '-c', '--commit',
            help='also commit and push to smodels.github.io (works only with -u github, anomaly, latest, or interesting)',
            action="store_true" )
    argparser.add_argument ( '--rundir',
            help='override the default rundir [None]',
            type=str, default=None )
    argparser.add_argument ( '--dbpath',
            help='path to database [<rundir>/default.pcl]',
            type=str, default="<rundir>/default.pcl" )
    argparser.add_argument ( "--destinations",
            help="learn more about the upload destinations", action="store_true" )
    args = argparser.parse_args()
    rundir = setup( args.rundir )
    if args.all:
        args.html = True
        args.ruler = True
        args.decays = True
        args.predictions = True
        args.tex = True
    if args.hiscorefile == "default":
        args.hiscorefile = f"{rundir}/hiscores.cache"
    runPlotting ( args )
    if args.test:
        compileTestText()

if __name__ == "__main__":
    main()
