#!/usr/bin/env python3

import pickle, os, sys, subprocess, time, glob, colorama, math, numpy
sys.path.insert(0,"../")
sys.path.insert(0,"/scratch-cbe/users/wolfgan.waltenberger/git/protomodels/")
sys.path.insert(0,"/scratch-cbe/users/wolfgan.waltenberger/git/smodels-utils/prototools/")
from ptools import hiscoreTools
from builder.manipulator import Manipulator
from builder import protomodel
from builder.protomodel import ProtoModel
from csetup import setup
from smodels.tools.physicsUnits import fb, TeV
from smodels.theory.theoryPrediction import TheoryPrediction
from smodels.tools import runtime
from smodels_utils.plotting import rulerPlotter, decayPlotter
from smodels_utils.helper.bibtexTools import BibtexWriter
from ptools.sparticleNames import SParticleNames
from smodels.tools.smodelsLogging import logger
logger.setLevel("ERROR")

runtime._experimental = True

def obtain ( number, picklefile ):
    """ obtain hiscore number <number>
    :returns: model
    """
    if not os.path.exists ( picklefile ):
        print ( "[plotHiscore] %s does not exist. Trying to produce now with ./hiscore.hi" % picklefile )
        from argparse import Namespace
        args = Namespace()
        args.detailed = False
        args.print = False
        args.outfile = "hiscore.hi"
        args.infile = picklefile
        args.fetch = False
        args.maxloss = 0.005
        hiscoreTools.main ( args )

    with open( picklefile,"rb" ) as f:
        #fcntl.flock( f, fcntl.LOCK_EX )
        hiscores = pickle.load ( f )
        #fcntl.flock( f, fcntl.LOCK_UN )
    if len(hiscores)<number-1:
        print ( "[plotHiscore] asked for hiscore %d, but I only have %d" % \
                ( number, len(hiscores) ) )
    Z = hiscores[number].Z
    K = hiscores[number].K
    print ( "[plotHiscore] obtaining #%d: K=%.3f" % (number, K ) )
    return hiscores[ number ]

def gitCommit ( dest, upload, wanted ):
    """ if wanted, then git commit and git push to smodels.githuib.io
    :param dest: e.g. "/scratch-cbe/users/w.w/git/smodels.github.io/protomodels/latest/"
    :param upload: e.g. "latest"
    :param wanted: were we even asked to git-commit?
    """
    if not wanted:
        return False
    destdir = dest
    destdir = destdir.replace(upload,"")
    destdir = destdir.replace("//","")

    comment = "automated update by plotHiscore.py to %s" % destdir
    print ( "destdir", destdir, "upload", upload, "wanted", wanted )
    cmd = "cd %s; git pull; git add '%s'; git commit -m '%s'; git push" % \
           ( destdir, upload, comment )
    print ( "[plotHiscore] now git-commit: %s" % cmd )
    out = subprocess.getoutput ( cmd )
    if out != "":
        print ( "[plotHiscore] %s" % out )
    return True

def discussPredictions ( protomodel ):
    print ( "How the K comes about. Best combo:" )
    combo = protomodel.bestCombo
    for pred in combo:
        print ( "theory pred: %s:%s" % ( pred.expResult.globalInfo.id, ",".join ( map ( str, pred.txnames ) ) ) )
        # print ( "     `- ", pred.expResult.globalInfo.id, "ana", pred.analysis, "masses", pred.mass, "txnames", pred.txnames, "type", pred.dataType() )

def getExtremeSSMs ( ssm, largest, nm = 7 ):
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
            s += "%s=%s; " % ( ssm[k], k )
    if len(s)>2:
        s = s[:-2]
    ret = "%signal strength multipliers: $%s$" % ( extreme, s )
    return ret

def hasSignals ( protomodel ):
    """ are there signals stored in the theory predictions? """
    for tp in protomodel.bestCombo:
        if hasattr ( tp.dataset, "dataInfo" ) and hasattr(tp.dataset.dataInfo,"sigN" ):
            return "did", True
        if hasattr ( tp.dataset, "txnameList" ):
            for txn in tp.dataset.txnameList:
                if hasattr ( txn, "sigmaN" ):
                    return "did", True
    return "did not", False

def writeRawNumbersHtml ( protomodel ):
    """ write out the raw numbers of the excess, as html """
    f=open("rawnumbers.html","wt")
    f.write("<table>\n" )
    f.write("<tr><th>Analysis Name</th><th>Type</th><th>Dataset</th><th>Topos</th><th>Observed</th><th>Expected</th><th>Approx &sigma;</th><th>Particles</th>" )
    didordidnot,hassigs = hasSignals ( protomodel )
    print ( f"[plotHiscore] protomodel's database {didordidnot} have fake signals." )
    if hassigs:
        f.write("<th>Signal</th>" )
    f.write("\n</tr>\n" )
    namer = SParticleNames ( susy = False )
    for tp in protomodel.bestCombo:
        anaId = tp.analysisId()
        idAndUrl = anaNameAndUrl ( tp )
        dtype = tp.dataType()
        ltopos = list(map(str,tp.txnames ) )
        ltopos.sort()
        topos = ", ".join ( ltopos )
        S = "?"
        dt = { "upperLimit": "ul", "efficiencyMap": "em", "combined": "comb" }
        f.write ( "<tr><td>%s</td><td>%s</td> " % ( idAndUrl, dt[dtype] ) )
        #f.write ( "<tr><td>%s</td><td>%s</td> " % ( anaId, dt[dtype] ) )
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
            toterr = math.sqrt ( bgErr**2 + eBG )
            if toterr > 0.:
                S = "%.1f &sigma;" % ( (dI.observedN - eBG ) / toterr )
            pids = set()
            for prod in tp.PIDs:
                for branch in prod:
                    for pid in branch:
                        if type(pid) == int and abs(pid)!=1000022:
                            pids.add ( abs(pid) )
                        if type(pid) in [ list, tuple ] and abs(pid[0])!=1000022:
                            pids.add ( abs(pid[0]) )
            obsN = dI.observedN
            if ( obsN - int(obsN) ) < 1e-6:
                obsN=int(obsN)
            particles = namer.htmlName ( pids, addSign = False,
                                          addBrackets = False )
            sobsN = str(obsN)
            if type(obsN) == float:
                sobsN = "%.2f" % obsN
            f.write ( '<td>%s</td><td>%s</td><td>%s</td><td>%s +/- %s</td><td style="text-align:right">%s</td><td style="text-align:right">%s</td>' % \
                      ( did, topos, sobsN, eBG, bgErr, S, particles ) )
            if hassigs:
                sig = "-"
                if hasattr ( dI, "sigN" ):
                    sig = "%s" % dI.sigN
                    if type(dI.sigN) in [ float ]:
                        sig = "%.2f" % dI.sigN
                f.write ( '<td style="text-align:right">%s</td>' % sig )
        if dtype == "upperLimit":
            S = "?"
            llhd = tp.likelihood( expected=False )
            eUL = tp.getUpperLimit ( expected = True ).asNumber(fb)
            oUL = tp.getUpperLimit ( expected = False ).asNumber(fb)
            sigma_exp = eUL / 1.96 # the expected scale, sigma
            Z = ( oUL - eUL ) / sigma_exp
            # Z = math.sqrt ( chi2 )
            S = "%.1f &sigma;" % Z
            # S = "%.2g l" % llhd
            # print ( "llhd,chi2,Z", llhd,chi2,Z )
            # p = 1. - scipy.stats.chi2.cdf ( chi2, df=1 )
            pids = set()
            for prod in tp.PIDs:
                for branch in prod:
                    for pid in branch:
                        if type(pid) == int and abs(pid)!=1000022:
                            pids.add ( abs(pid) )
                        if type(pid) in [ list, tuple ] and abs(pid[0])!=1000022:
                            pids.add ( abs(pid[0]) )
            #particles = helpers.toHtml ( pids, addSign = False,
            #                              addBrackets = False )
            particles = namer.htmlName ( pids, addSign = False, addBrackets = False )
            f.write ( '<td>-</td><td>%s</td><td> %.1f fb </td><td> %.1f fb</td><td style="text-align:right">%s</td><td style="text-align:right">%s</td>' % \
                    ( topos, tp.getUpperLimit().asNumber(fb), tp.getUpperLimit ( expected = True ).asNumber(fb),
                      S, particles ) )
            if hassigs:
                sig = "-"
                for txn in tp.txnames:
                # for txn in tp.dataset.txnameList:
                    if hasattr ( txn, "sigmaN" ):
                        sig = "%.2f fb" % txn.sigmaN
                f.write ( '<td style="text-align:right">%s</td>' % sig )
        if dtype == "combined":
            S = "?"
            llhd = tp.likelihood( expected=False )
            eUL = tp.getUpperLimit ( expected = True ).asNumber(fb)
            oUL = tp.getUpperLimit ( expected = False ).asNumber(fb)
            sigma_exp = eUL / 1.96 # the expected scale, sigma
            Z = ( oUL - eUL ) / sigma_exp
            # Z = math.sqrt ( chi2 )
            S = "%.1f &sigma;" % Z
            # S = "%.2g l" % llhd
            # print ( "llhd,chi2,Z", llhd,chi2,Z )
            # p = 1. - scipy.stats.chi2.cdf ( chi2, df=1 )
            pids = set()
            for prod in tp.PIDs:
                for branch in prod:
                    for pid in branch:
                        if type(pid) == int and abs(pid)!=1000022:
                            pids.add ( abs(pid) )
                        if type(pid) in [ list, tuple ] and abs(pid[0])!=1000022:
                            pids.add ( abs(pid[0]) )
            particles = namer.htmlName ( pids, addSign = False, addBrackets = False )
            f.write ( '<td>-</td><td>%s</td><td> %.1f fb </td><td> %.1f fb</td><td style="text-align:right">%s</td><td style="text-align:right">%s</td>' % \
                    ( topos, tp.getUpperLimit().asNumber(fb), tp.getUpperLimit ( expected = True ).asNumber(fb),
                      S, particles ) )
            if hassigs:
                sig = "-"
                for txn in tp.txnames:
                # for txn in tp.dataset.txnameList:
                    if hasattr ( txn, "sigmaN" ):
                        sig = "%.2f fb" % txn.sigmaN
                f.write ( '<td style="text-align:right">%s</td>' % sig )
        f.write ( '</tr>\n' )
    f.write("</table>\n" )
    f.close()

def writeRawNumbersLatex ( protomodel, usePrettyNames = True ):
    """ write out the raw numbers of the excess, in latex
    :param usePrettyNames: use the pretty names, not analysis ids
    """
    print ( "raw numbers of excess" )
    print ( "=====================" )
    f=open("rawnumbers.tex","wt")
    f.write("\\begin{tabular}{l|c|r|r|c|r|r}\n" )
    f.write("\\bf{Analysis Name} & \\bf{Dataset} & \\bf{Obs} & \\bf{Expected} & \\bf{Z} & \\bf{Particle} & \\bf{Signal} \\\\\n" )
    f.write("\\hline\n" )
    namer = SParticleNames ( susy = False )
    bibtex = BibtexWriter()
    from smodels_utils.helper.prettyDescriptions import prettyTexAnalysisName
    for tp in protomodel.bestCombo:
        anaId = tp.analysisId()
        ananame = anaId
        if usePrettyNames:
            ananame = prettyTexAnalysisName ( None, anaid = anaId )
        dtype = tp.dataType()
        print ( "[plotHiscore] item %s (%s)" % ( anaId, dtype ) )
        dt = { "upperLimit": "ul", "efficiencyMap": "em" }
        # f.write ( "%s & %s & " % ( anaId, dt[dtype] ) )
        ref = bibtex.query ( anaId )
        f.write ( "%s~\\cite{%s} & " % ( ananame, ref ) )
        if dtype == "efficiencyMap":
            dI = tp.dataset.dataInfo
            obsN = dI.observedN
            if ( obsN - int(obsN) ) < 1e-6:
                obsN=int(obsN)
            print ( "  `- %s: observedN %s, bg %s +/- %s" % \
                    ( dI.dataId, obsN, dI.expectedBG, dI.bgError ) )
            did = dI.dataId.replace("_","\_")
            if len(did)>9:
                did=did[:6]+" ..."
            eBG = dI.expectedBG
            if eBG == int(eBG):
                eBG=int(eBG)
            bgErr = dI.bgError
            if bgErr == int(bgErr):
                bgErr=int(bgErr)
            toterr = math.sqrt ( bgErr**2 + eBG )
            if toterr > 0.:
                S = "%.1f $\sigma$" % ( (dI.observedN - eBG ) / toterr )
            # pids = tp.PIDs
            pids = set()
            for prod in tp.PIDs:
                for branch in prod:
                    for pid in branch:
                        if type(pid) == int and abs(pid)!=1000022:
                            pids.add ( abs(pid) )
                        if type(pid) in [ list, tuple ]:
                            p = abs(pid[0])
                            if p!=1000022:
                                pids.add ( p )
            particles = namer.texName ( pids, addDollars=True, addSign = False,
                                          addBrackets = False )
            obs = dI.observedN
            if obs == 0.:
                obs = 0
            else:
                if abs ( obs - int(obs) ) / obs < 1e-6:
                    obs = int ( obs )
            sigN = tp.xsection.value.asNumber(fb) * tp.dataset.globalInfo.lumi.asNumber(1/fb)
            #sigmapred="%.2f fb" % ( tp.xsection.value.asNumber(fb) )
            sigmapred="%.2f" % sigN
            f.write ( "%s & %s & %s $\\pm$ %s & %s & %s & %s \\\\ \n" % \
                      ( did, obs, eBG, bgErr, S, particles, sigmapred ) )
        if dtype == "upperLimit":
            S = "?"
            llhd = tp.likelihood ( expected=False )
            eUL = tp.getUpperLimit ( expected = True ).asNumber(fb)
            oUL = tp.getUpperLimit ( expected = False ).asNumber(fb)
            sigma_exp = eUL / 1.96 # the expected scale, sigma
            Z = ( oUL - eUL ) / sigma_exp
            # Z = math.sqrt ( chi2 )
            S = "%.1f $\sigma$" % Z
            pids = set()
            for prod in tp.PIDs:
                for branch in prod:
                    for pid in branch:
                        if type(pid)==int and abs(pid)!=1000022:
                            pids.add ( abs(pid) )
                        if type(pid) in [ tuple, list ]:
                            for p in pid:
                                if type(p)==int and abs(p)!=1000022:
                                    pids.add ( abs(p) )
            particles = namer.texName ( pids, addDollars=True, addSign = False,
                                        addBrackets = False )
            sigmapred="%.2f fb" % ( tp.xsection.value.asNumber(fb) )
            print ( f"  `- observed {oUL}, expected {eUL}" )
            f.write ( f" & {oUL:.1f} fb & {eUL:.1f} fb & {S} & {particles} & {sigmapred} \\\\ \n" )
    f.write("\end{tabular}\n" )
    f.close()

def findXSecOfPids ( xsecs, pids ):
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
            return xsec.value
    return 0.*fb

def writeTex ( protomodel, keep_tex ):
    """ write the comment about ss multipliers and particle contributions, in tex.
    Creates texdoc.png.
    :param keep_tex: keep tex source of texdoc.png
    """
    cpids = {}
    frozen = protomodel.frozenParticles()
    xsecs = protomodel.getXsecs()[0]
    ssms = getUnfrozenSSMs ( protomodel, frozen, includeOnes=True )
    namer = SParticleNames ( susy = False )
    for pids,v in ssms.items():
        xsec = findXSecOfPids ( xsecs, pids )
        if xsec < 0.001 * fb: ## only for xsecs we care about
            continue
        sv = "%.2g" % v
        if v < .1:
            sv = "%.1g" % v
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
    if hasattr ( protomodel, "particleContributions" ):
        print ( "[plotHiscore] contributions-by-particle are defined" )
        #particleContributionList+="\\\\Contributions by particles: $"
        particleContributionList+="\\\\"
        particleContributionList+="Contributions by particles: $"
        totalcont = 0. ## to normalize contributions
        for k,v in protomodel.particleContributions.items():
            totalcont += (protomodel.K - v)
        tok = {}
        for k,v in protomodel.particleContributions.items():
            if v in tok.keys():
                v+=1e-6
            perc = 100.
            if totalcont != 0.:
                perc = round(100.*(protomodel.K - v)/totalcont )
            tok[v] = "%s: K_\mathrm{without}=%.2f (%d%s)" % ( namer.texName(k), v, perc, "\%" )
            # tok[v] = "%s = (%.2f) %d%s" % ( namer.texName(k), v, perc, "\%" )
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
    src = getExtremeSSMs ( ssm, largest=True, nm = 7 )
    src += "\\\\"
    nsmallest = min( 7, len(ssm)-7 )
    if nsmallest>0:
        src += getExtremeSSMs ( ssm, largest=False, nm = nsmallest )
    src += particleContributionList
    if keep_tex:
        with open("texdoc.tex","wt") as f:
            f.write ( src+"\n" )
            f.close()
        print ( "[plotHiscore] wrote %s/texdoc.tex" % os.getcwd() )
    try:
        p = tex2png.Latex ( src, 600 ).write()
        f = open ( "texdoc.png", "wb" )
        f.write ( p[0] )
        f.close()
    except Exception as e:
        print ( "[plotHiscore] Exception when latexing: %s" % e )
    return src

def anaNameAndUrl ( ana, forPdf=False, protomodel=None ):
    """ given analysis, return analysis name and URL,
    as html code or pdf hyperref
    :param ana: ExpRes or TheoryPred or str object
    :param forPdf: if True, create for Pdf hyperref, else for html
    :param protomodel: needed only if ana is str
    """
    if forPdf:
        ## FIXME implement!
        return ana.analysisId()
    if type(ana)==str:
        for i in protomodel.bestCombo:
            url = ""
            if i.analysisId() == ana and hasattr ( i.dataset.globalInfo, "url" ):
                url = i.dataset.globalInfo.url
                break
        return "<a href=%s>%s</a>" % \
               ( url, ana )
    if type(ana)==TheoryPrediction:
        if not hasattr ( ana.dataset.globalInfo, "url" ):
            return ( ana.analysisId() )
        return "<a href=%s>%s</a>" % \
               ( ana.dataset.globalInfo.url, ana.analysisId() )

def getPrettyName ( rv ):
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

def writeRValuesTex ( rvalues, usePrettyNames = True ):
    """ write out the leading rvalues of the critic, in latex
    :param usePrettyNames: use the pretty names, not analysis ids
    """
    namer = SParticleNames ( False )
    g=open("rvalues.tex","wt")
    g.write ( "\\begin{tabular}{l|c|c|c|c|c}\n" )
    g.write ( "\\bf{Analysis Name} & \\bf{Production} & $\sigma_{XX}$ (fb) & $\sigma^\mathrm{UL}_\mathrm{obs}$ (fb) & $\sigma^\mathrm{UL}_\mathrm{exp}$ (fb) & $r$ \\\\\n" )
    #g.write ( "\\begin{tabular}{l|c|r|r}\n" )
    #g.write ( "\\bf{Analysis Name} & \\bf{Topo} & $r_{\mathrm{obs}}$ & $r_{\mathrm{exp}}$ \\\\\n" )
    g.write ( "\\hline\n" )
    bibtex = BibtexWriter()
    from smodels_utils.helper.prettyDescriptions import prettyTexAnalysisName
    for rv in rvalues[:5]:
        srv="N/A"
        if type(rv[1]) in [ float, numpy.float64, numpy.float32 ]:
            srv="%.2f" % rv[1]
        else:
            srv=str(rv[1])
        anaId = rv[2].analysisId()
        prettyName = getPrettyName( rv[2] )
        prettyName = prettyTexAnalysisName ( prettyName, anaid = anaId )
        ref = bibtex.query ( anaId )
        txnames = ",".join ( map(str,rv[2].txnames) )
        allpids = rv[2].PIDs
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
        sigmapred = rv[2].xsection.value.asNumber(fb)
        sigmaexp = "--"
        if type(rv[2].getUpperLimit ( expected = True )) != type(None):
            sigmaexp = "%.2f" % rv[2].getUpperLimit ( expected=True ).asNumber(fb)
        sigmaobs = rv[2].getUpperLimit().asNumber(fb)
        g.write ( "%s~\\cite{%s} & %s & %.2f & %.2f & %s & %.2f\\\\\n" % \
                ( prettyName, ref, prod, sigmapred, sigmaobs, sigmaexp, rv[0] ) )
    g.write ( "\\end{tabular}\n" )
    g.close()

def getDatabaseVersion ( protomodel, dbpath = "default.pcl" ):
    dbver = "???"
    if hasattr ( protomodel, "dbversion" ):
        dbver = protomodel.dbversion
        if not "???" in dbver:
            return dbver
    if os.path.exists ( dbpath ):
        ## try to get db version from db file
        from smodels.experiment.databaseObj import Database
        db = Database ( dbpath )
        dbver = db.databaseVersion
    return dbver

def writeIndexTex ( protomodel, texdoc ):
    """ write the index.tex file
    :param texdoc: the source that goes into texdoc.png
    """
    ssm = []
    namer = SParticleNames ( False )
    for k,v in protomodel.ssmultipliers.items():
        if abs(v-1.)<1e-3:
            continue
        ssm.append ( "%s: %.2f" % ( namer.texName(k,addSign=True),v) )
    f=open("index.tex","w")
    f.write ( "Our current winner has a score of \\K=%.2f, " % \
              ( protomodel.K ) )
    strategy = "aggressive"
    dbver = getDatabaseVersion ( protomodel )
    dotlessv = dbver.replace(".","")
    f.write ( " it was produced with database {\\tt v%s}, combination strategy {\\tt %s} walker %d in step %d." % \
            ( dotlessv, strategy, protomodel.walkerid, protomodel.step ) )
    f.write ( "\n" )
    if hasattr ( protomodel, "tpList" ):
        rvalues=protomodel.tpList
        rvalues.sort(key=lambda x: x[0],reverse=True )
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
            f.write ( "%s & %.2f & %s%s \\\\ \n" % ( k, Kwithout, int(round(100.*cont)), "\\%" ) )
            # f.write ( "\item %s: %s%s\n" % ( k, int(round(100.*v)), "\\%" ) )
        # f.write ( "\end{itemize}\n" )
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
    # f.write ( "<td><img width=600px src=./texdoc.png>\n" ) #  % height )
    # f.write ( "\small{Last updated: %s}\n" % time.asctime() )
    # f.write ( "% include decays.png\n" )
    contrs = texdoc.replace(":"," are " ).replace("S","The s").replace(";",", " )
    contrs = contrs.replace( "\\\\\\\\Contr", "; the contr" )
    f.write ( contrs + "\n" )
    f.close()
    print ( "[plotHiscore] Wrote index.tex" )

def getUnfrozenSSMs ( protomodel, frozen, includeOnes=False, dropLSPLSP=True ):
    """ of all SSMs, leave out the ones with frozen particles
    :param protomodel: the, well, protomodel
    :param frozen: list of pids of frozen particles
    :param includeOnes: if False, then also filter out values close to unity
    :param dropLSPLSP: if True, dont include LSP,LSP production
    :returns: dictionary of SSMs without frozen particles
    """
    # ssms = protomodel.ssmultipliers
    ma = Manipulator ( protomodel )
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


def writeIndexHtml ( protomodel ):
    """ write the index.html file, see e.g.
        https://smodels.github.io/protomodels/
    """
    ssm = []
    namer = SParticleNames ( susy = False )
    frozen = protomodel.frozenParticles()
    ssms = getUnfrozenSSMs ( protomodel, frozen, includeOnes=True )
    for k,v in ssms.items():
        ssm.append ( "%s: %.2g" % ( namer.htmlName(k,addSign=True),v) )
    f=open("index.html","w")
    f.write ( "<html>\n" )
    f.write ( "<body>\n" )
    f.write ( "<center>\n" )
    f.write ( "<table><td><h1>Current best protomodel: <i>K</i>=%.2f</h1><td><img height=60px src=https://smodels.github.io/pics/banner.png></table>\n" % ( protomodel.K ) )
    f.write ( "</center>\n" )
    dbver = getDatabaseVersion ( protomodel )
    strategy = "aggressive"
    dotlessv = dbver.replace(".","")
    dt = int ( time.time() - 1593000000 )
    f.write ( "<b><a href=./hiscore.slha>ProtoModel</a> <a href=./pmodel.py>(dict)</a> produced with <a href=https://smodels.github.io/docs/Validation%s>database v%s</a>, combination strategy <a href=./matrix_%s.png>%s</a> in walker %d step %d.</b> " % \
            ( dotlessv, dbver, strategy, strategy, protomodel.walkerid, protomodel.step ) )
    if hasattr ( protomodel, "particleContributions" ):
        f.write ( "<i>K</i> plots for: <a href=./M1000022.png?%d>%s</a>" % \
                  ( dt, namer.htmlName(1000022) ) )
        for k,v in protomodel.particleContributions.items():
            f.write ( ", " )
            f.write ( "<a href=./M%d.png?%d>%s</a>" % ( k, dt, namer.htmlName(k) ) )
        f.write ( ". HPD plots for: " )
        first = True
        for k,v in protomodel.particleContributions.items():
            if not first:
                f.write ( ", " )
            f.write ( "<a href=./llhd%d.png?%d>%s</a>" % ( k, dt, namer.htmlName(k) ) )
            first = False
    # fixme replace with some autodetection mechanism
    ossms = { (-1000006,1000006), (1000021,1000021), (-2000006,2000006) }
    for fname in glob.glob("ssm_*_*.png" ):
        pids = fname.replace("ssm_","").replace(".png","")
        pids = tuple ( map ( int, pids.split("_") ) )
        ossms.add ( pids )
    frozen = protomodel.frozenParticles()
    # take out all frozen ssm plots
    ssms = set()
    for pids in ossms:
        hasFrozenPid=False
        for pid in pids:
            if pid in frozen or -pid in frozen:
                hasFrozenPid = True
                break
        if not hasFrozenPid:
            ssms.add ( pids )


    f.write ( ". SSM plots for: " )
    first = True
    for pids in ssms:
        if not first:
            f.write ( ", " )
        f.write ( "<a href=./ssm_%d_%d.png?%d>(%s,%s)</a>" % \
                  ( pids[0],pids[1], dt, namer.htmlName(pids[0],addSign=True),
                    namer.htmlName(pids[1],addSign=True) ) )
        first = False
    f.write ( "<br>\n" )
    f.write ( "<table width=80%>\n<tr><td>\n" )
    if hasattr ( protomodel, "tpList" ):
        rvalues=protomodel.tpList
        rvalues.sort(key=lambda x: x[0],reverse=True )
        f.write ( f"<br><b>{len(rvalues)} predictions available. Highest r values are:</b><br><ul>\n" )
        for rv in rvalues[:5]:
            srv="N/A"
            if type(rv[1]) in [ float, numpy.float64, numpy.float32 ]:
                srv="%.2f" % rv[1]
            elif type(rv[1]) != type(None):
                srv=str(rv[1])
            f.write ( f"<li>{anaNameAndUrl ( rv[2] )}:{','.join ( set (map(str,rv[2].txnames) ) )} r={rv[0]:.2f}, r<sub>exp</sub>={srv}<br>\n" )
        f.write("</ul>\n")
    else:
        print ( "[plotHiscore] protomodel has no r values!" )

    if hasattr ( protomodel, "analysisContributions" ):
        print ( "[plotHiscore] contributions-per-analysis are defined" )
        f.write ( "<td><br><b>K values without analyses: (contribution)</b><br>\n<ul>\n" )
        conts = []
        Ktot = protomodel.K
        dKtot = 0.
        contributions = protomodel.analysisContributions
        for k,v in contributions.items():
            conts.append ( ( v, k ) )
            dKtot += Ktot - v
        conts.sort( reverse=True )
        for v,k in conts:
            Kwithout= contributions[k]
            cont = ( Ktot - Kwithout ) / dKtot
            nameAndUrl = anaNameAndUrl ( k, protomodel=protomodel )
            kv = str(v)
            if type(v) in [ float, numpy.float64 ]:
                kv = "%.2f (%d%s)" % ( v,int(round(100.*cont)), "%" )
            f.write ( "<li> %s: %s\n" % ( nameAndUrl, kv ) )
        # f.write ( "</table>\n" )
    else:
        print ( "[plotHiscore] analysis-contributions are not defined" )

    height = 32*int((len(ssm)+3)/4)
    if ssm == []:
        height = 32
    if hasattr ( protomodel, "particleContributions" ):
        height += 32
    t0 = int(time.time())
    f.write ( "<td><img width=600px src=./texdoc.png?%d>\n" % ( t0 ) )
    f.write ( "<br><font size=-1>Last updated: %s</font>\n" % time.asctime() )
    f.write ( "</table>" )
    f.write ( '<table style="width:80%">\n' )
    f.write ( "<td width=45%>" )
    f.write ( "<img height=580px src=./ruler.png?%d>" % ( t0 ) )
    f.write ( "<td width=55%>" )
    f.write ( "<img height=380px src=./decays.png?%d>\n" % ( t0 ) )
    f.write ( '<font size=-3><iframe type="text/html" height="220px" width="100%s" frameborder="0" src="./rawnumbers.html?%d"></iframe></font>\n' % ( "%s", t0 ) )
    f.write ( "</table>\n" )
    # f.write ( "<br><font size=-1>Last updated: %s</font>\n" % time.asctime() )
    f.write ( "</body>\n" )
    f.write ( "</html>\n" )
    f.close()
    print ( "[plotHiscore] Wrote index.html" )

def copyFilesToGithub():
    files = [ "hiscore.slha", "index.html", "matrix_aggressive.png", "decays.png",
              "ruler.png", "texdoc.png", "pmodel.py", "rawnumbers.html" ]
    for f in files:
        if not os.path.exists ( f ):
            continue
        O = subprocess.getoutput ( "cp %s ../../smodels.github.io/protomodels/" % f )
        if len(O)>0:
            print ( "[plotHiscore.py] when copying files: %s" % O )

def getPIDsOfTPred ( tpred, ret, integrateDataType=True, integrateSRs=True ):
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
    for pids in tpred.PIDs:
        for br in pids:
            for pid in br:
                if type(pid) in [ list ]:
                    for pp in pid:
                        apid = abs(pp)
                        if not apid in ret and not apid == LSP:
                            ret[apid]=set()
                        if not apid == LSP:
                            ret[apid].add ( name )
                else:
                    apid = abs(pid)
                    if not apid in ret and not apid == LSP:
                        ret[apid]=set()
                    if not apid == LSP:
                        ret[apid].add ( name )
    return ret

def plotRuler( protomodel, verbosity, horizontal ):
    """ plot the ruler plot, given protomodel.
    :param horizontal: plot horizontal ruler plot, if true. Else, vertical.
    """
    fname = "ruler.png"
    if horizontal:
        fname = "horizontal.png"
    print ( "[plotHiscore] now draw %s" % fname )
    resultsForPIDs = {}
    for tpred in protomodel.bestCombo:
        resultsForPIDs =  getPIDsOfTPred ( tpred, resultsForPIDs )
        # print ( "p", pidsofpred )
        # resultsForPIDs.union ( getPIDsOfTPred ( tpred ) )
    resultsFor = {}
    for pid,values in resultsForPIDs.items():
        if pid in protomodel.masses:
            resultsFor[ protomodel.masses[pid] ] = values
        else:
            print ( "[plotHiscore] why is pid %s not in mass dict %s?" % ( pid, str(protomodel.masses) ) )

    if verbosity == "debug":
        print ( '[plotHiscore] ../smodels-utils/smodels_utils/plotting/rulerPlotter.py -o ruler.png --hasResultsFor "%s" %s' % \
                ( str(resultsFor), protomodel.currentSLHA ) )

    plotter = rulerPlotter.RulerPlot ( protomodel.currentSLHA, fname,
                                       Range=(None, None), mergesquark = False,
                                       drawdecays = False,
                                       hasResultsFor = resultsFor,
                                       trim = True )
    if horizontal:

        plotter.drawHorizontal()
    else:
        plotter.drawVertical()

def plotDecays ( protomodel, verbosity, outfile="decays.png" ):
    print ( "[plotHiscore] now draw %s" % outfile )
    options = { "tex": True, "color": True, "dot": True, "squarks": True,
                "weakinos": True, "sleptons": True, "neato": True,
                "separatecharm": True,
                "integratesquarks": False, "leptons": True }
    options["rmin"] = 0.005
    ## FIXME add cross sections.
    if verbosity == "debug":
        soptions = ""
        for k,v in options.items():
            if v==True:
                soptions += "--%s " % k
        ma = Manipulator ( protomodel )
        ssms = ma.simplifySSMs()
        # soptions+=' --ssmultipliers "%s"' % ssms
        print ( "%s../smodels_utils/plotting/decayPlotter.py -f %s -o %s %s%s" % \
                ( colorama.Fore.GREEN, protomodel.currentSLHA, outfile, soptions,
                  colorama.Fore.RESET ) )
    decayPlotter.draw ( protomodel.currentSLHA, outfile, options,
                        verbosity = verbosity,
                        ssmultipliers = protomodel.ssmultipliers )

def plot ( number, verbosity, picklefile, options, dbpath ):
    ## plot hiscore number "number"
    protomodel = obtain ( number, picklefile )

    protoslha = protomodel.createSLHAFile ()
    subprocess.getoutput ( "cp %s hiscore.slha" % protoslha )
    m = Manipulator ( protomodel )
    print ( "[plotHiscore] now write pmodel.py" )
    m.writeDictFile()
    opts = [ "ruler", "decays", "predictions", "copy", "html" ]
    for i in opts:
        if not i in options:
            options[i]=True

    plotruler = options["ruler"]
    horizontal = False
    if "horizontal" in options and options["horizontal"]:
        horizontal = True
    # print ( "[plotHiscore] slha file exists?", os.path.exists ( protomodel.currentSLHA ), protomodel.currentSLHA )
    if plotruler:
        plotRuler ( protomodel, verbosity, horizontal )
    plotdecays = options["decays"]
    if plotdecays:
        plotDecays ( protomodel, verbosity )

    if options["predictions"]:
        discussPredictions ( protomodel )
    if options["html"] or options["tex"]:
        texdoc = writeTex ( protomodel, options["keep_tex"] )
        if options["html"]:
            writeIndexHtml ( protomodel )
        if options["tex"]:
            writeIndexTex( protomodel, texdoc )
    writeRawNumbersLatex ( protomodel )
    writeRawNumbersHtml ( protomodel )
    protomodel.delCurrentSLHA()

def runPlotting ( args ):
    if args.destinations:
        print ( "Upload destinations: " )
        print ( "      none: no upload" )
        print ( "       gpu: upload to GPU server, afs space." )
        print ( "            Result can be seen at http://www.hephy.at/user/wwaltenberger/protomodels/" )
        print ( "    github: upload to github git directory." )
        print ( "            Result can be seen at https://smodels.github.io/protomodels" )
        print ( "interesting: upload to github git directory, 'interesting' folder." )
        print ( "             Result can be seen at https://smodels.github.io/protomodels/interesting" )
        print ( "latest: upload to github git directory, 'latest' folder." )
        print ( "             Result can be seen at https://smodels.github.io/protomodels/latest" )
        print ( "anomaly: upload to github git directory, 'anomaly' folder." )
        print ( "             Result can be seen at https://smodels.github.io/protomodels/anomaly" )
        return
    upload = args.upload# .lower()
    if upload in [ "none", "None", "" ]:
        upload = None

    options = { "ruler": args.ruler, "decays": args.decays,
                "predictions": args.predictions, "html": args.html,
                "keep_tex": args.keep, "tex": args.tex,
                "horizontal": args.horizontal }

    plot ( args.number, args.verbosity, args.picklefile, options, args.dbpath )
    if upload is None:
        return
    F = "*.png pmodel.py hiscore.slha index.html rawnumbers.html"
    dest = ""
    destdir = "%s/git" % os.environ["HOME"]
    dest = "%s/smodels.github.io/protomodels/%s/" % ( destdir, upload )
    if upload == "github":
        dest = "%s/smodels.github.io/protomodels/" % destdir
    if upload in [ "interesting", "anomaly", "latest" ]:
        dest = "%s/smodels.github.io/protomodels/%s/" % ( destdir, upload )
    if "paper" in upload:
        dest = "%s/smodels.github.io/protomodels/%s/" % ( destdir, upload )

    if dest != "":
        if not os.path.exists ( dest ):
            print ( f"[plotHiscore] directory {dest} does not exist. Trying to mkdir it." )
            cmd = f"mkdir {dest}" ## -p ?
            a = subprocess.getoutput ( cmd )
            if a != "":
                print ( "[plotHiscore] error when mkdir: %s" % a )

        print ( "[plotHiscore] copying to %s" % dest )
        cmd = "cp %s %s" % ( F, dest )
        a = subprocess.getoutput ( cmd )
        if a != "":
            print ( "error: %s" % a )
            sys.exit()
        r = gitCommit( dest, upload, args.commit )
        if not r:
            destdir = dest
            destdir = destdir.replace(upload,"")
            destdir = destdir.replace("//","")
            print ( "[plotHiscore] done. now please do yourself: " )
            print ( "cd %s" % destdir )
            print ( "git add %s" % upload )
            print ( "git commit -m 'manual update'" )
            print ( "git push" )
        return

    if upload == "gpu":
        import socket
        hostname = socket.gethostname()
        D = "/afs/hephy.at/user/w/wwaltenberger/www/protomodels"
        ## first the backup
        if "gpu" in hostname:
            ## make backup
            cmd = "cp %s/* %s/backup/" % ( D, D )
        else:
            cmd = "ssh hepgpu01.hephy.oeaw.ac.at cp %s/* %s/backup/" % ( D, D )
        print ( cmd )
        # now the new stuff
        O = subprocess.getoutput ( cmd )
        if len(O)>0:
            print ( "[plotHiscore.py] when uploading files: %s" % O )

        if "gpu" in hostname:
            ## make backup
            cmd = "cp %s/* %s/backup/" % ( D, D )
            subprocess.getoutput ( cmd )
            cmd = "cp %s %s" % (F, D )
        else:
            cmd = "scp %s hepgpu01.hephy.oeaw.ac.at:%s" % ( F, D )
        print ( cmd )
        O = subprocess.getoutput ( cmd )
        if len(O)>0:
            print ( "[plotHiscore.py] when uploading files: %s" % O )
        return
    print ( "error, dont know what to do with upload sink '%s'" % upload )

def compileTestText():
    subprocess.getoutput ( "pdflatex test.tex" )

def main ():
    import argparse
    argparser = argparse.ArgumentParser(
            description='hiscore proto-model plotter')
    argparser.add_argument ( '-n', '--number',
            help='which hiscore to plot [0]',
            type=int, default=0 )
    argparser.add_argument ( '-f', '--picklefile',
            help='pickle file to draw from [<rundir>/hiscore.hi]',
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
            help='keep latex files',
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
    if args.picklefile == "default":
        args.picklefile = "%s/hiscore.hi" % rundir
    runPlotting ( args )
    if args.test:
        compileTestText()

if __name__ == "__main__":
    main()
