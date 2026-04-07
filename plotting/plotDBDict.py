#!/usr/bin/env python3

""" plot the meta statistics of database.dict, see e.g.
    https://smodels.github.io/validation/300/significances.png """

import sys,os
dirname = os.path.dirname ( os.path.abspath ( __file__ ) )
dirname = dirname.replace("/plotting","").replace("/protomodels","")
sys.path.insert(0,dirname)

from smodels_utils.plotting import mpkitty as plt
from smodels_utils.helper import prettyDescriptions
import numpy as np
import os, glob, sys, math, fnmatch
from copy import deepcopy as cp
import scipy.stats
import matplotlib.mlab as mlab
from typing import Union, List, Dict

from ptools.helpers import computeP
from ptools.moreHelpers import namesForSetsOfTopologies
from protomodels.base.loggerbase import LoggerBase
from smodels_utils.helper.various import hasLLHD, removeAnaIdSuffices
from smodels_utils.helper.terminalcolors import GREEN, YELLOW, RESET, RED

def isSelectedCustomFunction ( anaid, txns : List ) -> bool:
    """ a generic hook to make it easy to select based on any requirement.
    just write the python code here. """
    mets = {'CMS-SUS-16-034','ATLAS-SUSY-2018-22','CMS-SUS-19-009','CMS-SUS-19-011','ATLAS-SUSY-2018-10','CMS-SUS-16-048','CMS-SUS-13-019','ATLAS-SUSY-2013-20','CMS-SUS-19-010','CMS-SUS-16-035','ATLAS-SUSY-2013-08','ATLAS-SUSY-2013-15','ATLAS-SUSY-2016-28','ATLAS-SUSY-2018-05','ATLAS-SUSY-2013-23','CMS-EXO-19-010','CMS-SUS-13-013','CMS-SUS-16-033','CMS-SUS-19-013','CMS-SUS-12-028','CMS-SUS-16-041','CMS-SUS-16-047','ATLAS-SUSY-2017-02','ATLAS-SUSY-2016-15','CMS-SUS-17-005','CMS-SUS-13-006','CMS-SUS-16-009','ATLAS-SUSY-2016-14','CMS-SUS-21-002','CMS-EXO-19-001','ATLAS-SUSY-2014-03','ATLAS-SUSY-2015-09','ATLAS-SUSY-2018-04','CMS-SUS-16-046','ATLAS-SUSY-2019-02','ATLAS-SUSY-2013-18','ATLAS-SUSY-2015-01','CMS-PAS-SUS-16-052','CMS-SUS-20-002','CMS-SUS-14-010','ATLAS-SUSY-2013-12','CMS-SUS-16-037','CMS-SUS-18-007','CMS-SUS-16-051','ATLAS-SUSY-2013-05','ATLAS-SUSY-2016-27','ATLAS-SUSY-2013-09','ATLAS-SUSY-2016-16','ATLAS-SUSY-2018-23','CMS-SUS-13-011','ATLAS-SUSY-2016-24','CMS-PAS-SUS-13-023','CMS-SUS-17-004','ATLAS-SUSY-2013-02','CMS-SUS-20-001','ATLAS-SUSY-2018-42','CMS-SUS-20-004','ATLAS-SUSY-2018-41','ATLAS-SUSY-2018-14','ATLAS-EXOT-2018-06','ATLAS-SUSY-2013-21','ATLAS-SUSY-2013-04','CMS-SUS-16-043','CMS-PAS-SUS-13-018','ATLAS-SUSY-2016-06','ATLAS-SUSY-2018-31','ATLAS-SUSY-2013-19','CMS-EXO-20-004','ATLAS-SUSY-2013-16','CMS-SUS-16-032','CMS-SUS-16-050','ATLAS-SUSY-2017-03','CMS-SUS-19-006','CMS-SUS-14-021','CMS-SUS-18-004','ATLAS-SUSY-2015-02','ATLAS-SUSY-2019-08','ATLAS-SUSY-2016-08','CMS-SUS-16-042','CMS-SUS-18-002','CMS-SUS-16-045','CMS-SUS-17-009','ATLAS-SUSY-2016-33','ATLAS-SUSY-2018-09','CMS-SUS-13-002','ATLAS-SUSY-2018-40','CMS-SUS-19-008','CMS-SUS-17-010','CMS-SUS-13-007','ATLAS-SUSY-2013-11','ATLAS-SUSY-2018-08','ATLAS-SUSY-2016-07','ATLAS-SUSY-2016-26','CMS-SUS-13-012','CMS-SUS-21-007','CMS-SUS-16-036','ATLAS-SUSY-2019-09','CMS-EXO-13-006','ATLAS-SUSY-2015-06','ATLAS-SUSY-2018-12','CMS-SUS-17-006','ATLAS-SUSY-2017-01','ATLAS-SUSY-2016-19','ATLAS-SUSY-2016-32','CMS-SUS-12-024','ATLAS-SUSY-2018-32','ATLAS-SUSY-2016-17','CMS-SUS-13-004','ATLAS-SUSY-2018-16','ATLAS-SUSY-2018-06','CMS-PAS-SUS-13-015','CMS-SUS-16-039','CMS-PAS-SUS-13-016','CMS-SUS-17-003'}
    return anaid in mets
    #print ( f"isSelectedCustomFunction {anaid} {txns}" )
    #for tx in txns:
    #    pass
    #return False

class Plotter ( LoggerBase ):
    """ the meta statistics plotter, see eg https://smodels.github.io/validation/300/significances.png
    """

    def __init__ ( self, args : Dict ):
        """
        :param filename: filename of dictionary
        :param filtervalue: filter out signal regions with expectedBG < filtervalue
        :param comment: an optional comment, to write in the plot
        :param likelihood: form of likelihood: "gauss", "gauss+poisson", or
        "lognormal+poisson"
        "gauss" or "g" means only a Gaussian for everything
        "gauss+poisson" or "gp" means Gauss * Poisson
        "lognormal+poisson" or "lp" means Lognormal * Poisson
        :param topologies: if not None, then filter for these topologies (e.g. T2tt)
        :param unscale: unscale, i.e. use the fudged bgError also for computing likelihoods
        :param signalmodel: use the signal+bg model for computing likelihoods
        :param filtersigma: filter out signal regions with 
        expectedBG/bgErr < filtersigma
        :param collaboration: select a specific collaboration
        :param fakes: add fakes to the plot
        :param analyses: if not None, then filter for these analyses
                         (e.g. CMS-SUS-16-039-ma5)
        :param disclaimer: add a disclaimer, "do not circulate"
        :param ulAlso: show UL results, also
        :param title: a title
        :param pvalues: if true then plot p-values, if false plot significances
        """
        super ( Plotter, self ).__init__ ( 0 )
        self.defaults()
        for a,value in args.items():
            if a=="options":
                self.options.update ( value )
                continue
            if a not in [ "topologies", "analyses", "select_topologies" ]:
                setattr ( self, a, value )
        if self.nbins == None:
            if "nbins" in self.options:
                self.nbins = self.options["nbins"]
            else:
                if not self.pvalues:
                    self.nbins = 29
                    # self.nbins = 13
                else:
                    self.nbins = 30
                    # self.nbins = 10
        self.origtopos = args["topologies"]
        self.printMsgs = { "pass": 0 }
        if self.origtopos == None:
            self.origtopos = "all"
        self.description = None
        if args['topologies'] != None:
                if args["topologies"].endswith ( ".py" ):
                    print ( f"[plotDBDict] you supplied {args['topologies']} as topologies. Did you supply the validation file instead?" )
                args["topologies"], descr = namesForSetsOfTopologies ( args['topologies'] )
                self.description = descr
                if args["topologies"].startswith("^"):
                    self.description = f"all but {descr} searches"
        if "select_topologies" in args and args['select_topologies'] not in  [ None, [] ]:
                args["select_topologies"], descr = namesForSetsOfTopologies ( args['select_topologies'] )
        if 'select_collaboration' in args:
            collaboration = args['select_collaboration'].upper()
            if collaboration in [ "", "*" ]:
                collaboration = "ALL"
            if not collaboration in [ "CMS", "ATLAS", "ALL" ]:
                print ( "[plotDBDict] error: collaboration must be either CMS, ATLAS, or ALL." )
                sys.exit(-1)
            self.collaboration = collaboration
        if "likelihood" in args:
            abbreviations = { "g": "gauss", "gp": "gauss+poisson", "lp": "lognormal+poisson" }
            likelihood = args['likelihood']
            if likelihood in abbreviations:
                likelihood = abbreviations[likelihood]
            if likelihood not in [ "gauss", "gauss+poisson", "lognormal+poisson" ]:
                self.pprint ( "error, likelihood is to be one of: gauss, gauss+poisson, lognormal+poisson" )
                sys.exit()
            self.likelihood = likelihood ## False: gauss, True: lognormal
        self.verbose = 1
        if "verbose" in args:
            v = args["verbose"]
            if v == True:
                v = 2
            self.verbose = v
        if "ulalso" in args:
            self.useAlsoULMaps = args['ulalso']
        topologies = args['topologies']
        if topologies not in [ None, "" ]:
            topos = topologies.split(",")
            for t in topos:
                if t.startswith ( "^" ):
                    self.negative_topologies.append ( t[1:] )
                else:
                    self.topologies.append ( t )
        if "select_topologies" in args:
            select_topologies = args['select_topologies']
            if select_topologies not in  [ None, "", [] ]:
                topos = select_topologies.split(",")
                for t in topos:
                    self.select_topologies.append ( t )
        if len(self.select_topologies)>0 and self.ignore_sqrts == False:
            self.pprint ( f"select_topologies is on, we will set ignore_sqrts to True" )
            self.ignore_sqrts = True

        if "analyses" in args and args["analyses"] not in [ None ]:
            analyses = args['analyses']
            if analyses not in [ None, "" ]:
                analyses = analyses.split(",")
                for a in analyses:
                    if a.startswith ( "^" ):
                        self.negativeanalyses.append ( a[1:] )
                    else:
                        self.analyses.append ( a )
        self.filenames = []
        if "comment" in args:
            comment = args['comment']
            if comment in [ "None", "", "none" ]:
                comment = None
            self.comment = comment
        if "dictfile" in args:
            pathname = args['dictfile']
            if type(pathname) in [ str ]:
                pathname = pathname.split(",")
            for pname in pathname:
                if os.path.isdir ( pname ):
                    pname = f"{pname}/*dict"
                    self.filenames += glob.glob ( pname )
                else:
                    self.filenames.append ( pname )
        else:
            self.pprint ( "we need dictfile in args" )
        if "outfile" in args:
            self.outfile = self.determineOutFile ( args["outfile"] )
        self.meta = {}
        self.data = {}
        self.read()
        if self.Zmax is None:
            self.determineZmax()
        self.plot( )
        self.logExecution ( )

    def showPng ( self, filename : os.PathLike ):
        if ("show" in self.options and self.options["show"]==True) or self.show == True:
            from smodels_utils.plotting.mpkitty import timg
            timg ( filename )

    def logExecution ( self ):
        """ log the call of the executable """
        cmd = ""
        prev=""
        for i,a in enumerate(sys.argv):
            if i > 0:
                cmd += " "
            if "select_t" in prev or "--ti" in prev or "--op" in prev or prev in [ "-O", "-T" ]:
                a = f'"{a}"'
            cmd += a
            prev = a
        cmd += "\n"
        with open ( "plotDBDict.log", "at" ) as f:
            import time
            f.write ( f"\n# {time.asctime()}\n" )
            f.write ( cmd )
            f.close()

    def getBins ( self, nbins : Union[None,int] = None ):
        """ get the bin edges """
        if nbins == None:
            nbins = self.nbins
        step = 1/nbins
        bins = np.arange ( 0., 1+1e-7, step )
        if not self.pvalues:
            step = (2*self.Zmax)/nbins
            bins = np.arange ( -self.Zmax, self.Zmax+1e-7, step )
        return step, bins

    def defaults ( self ):
        self.nbins = None # 10 for p-values, 13 for significances
        self.fudge = 1.
        self.Zmax = None
        self.use_custom_function= False
        self.before = None
        self.nosuperseded = False # yes superseded
        self.use_aggregated = False # yes add non-aggregated
        self.nofastlim = False # yes fastlim
        self.show = False
        self.pvalues = False # if False, then p-values if true then significances
        self.skippedAgg = set() # log all aggregated analyses that have been skipped
        self.origtopos = "all"
        self.collaboration = "ALL"
        self.likelihood = "gauss+poisson"
        self.useAlsoULMaps = False
        self.analyses = []
        self.comment = None
        self.topologies = []
        self.select_topologies = []
        self.ignore_sqrts = False
        self.negative_topologies = []
        self.negativeanalyses = []
        self.outfile = "not_specified.png"
        self.title = None
        self.options = { "alwayslegend": False }
        self.yrange = None
        self.unscale = False
        self.signalmodel = False
        self.fakes = False
        self.sqrts = [8,13,13.6]
        self.filter = 0.
        self.filtersigma = 0.
        self.disclaimer = False
        self.warnedBefore = False # warned that --before is set

    def determineZmax ( self ):
        """ obtain self.Zmax from data """
        Zmax = 0.
        for dictfile,filecontent in self.data.items():
            for anaid,values in filecontent.items():
                if not "orig_Z_fudged" in values:
                    # if no fudged version in values, then use original version
                    # self.pprint ( "no orig_Z_fudged in dictionaries, did you forget the '-C' flag when calling expResModifier.py?" )
                    # sys.exit(-1)
                    if values["orig_Z"] > Zmax and np.isfinite ( values["orig_Z"] ):
                        Zmax = values["orig_Z"]
                else:
                    if values["orig_Z_fudged"] > Zmax and np.isfinite ( values["orig_Z_fudged"] ):
                        Zmax = values["orig_Z_fudged"]
        self.Zmax = np.ceil ( Zmax * 4. ) / 4.
        # self.pprint ( f"Zmax is {Zmax:.2f}" )

    def selectedCollaboration( self, anaid : str ) -> bool:
        """ does anaid pass the collaboration selection? """
        if self.collaboration in [ "ALL", "all", "*" ]:
            return True
        if self.collaboration in anaid:
            return True
        return False

    def filterByTime ( self, D : Dict ) -> bool:
        """ filter by time, let everything before self.before pass! 
        :returns: true means pass the filter
        """
        if self.before == None:
            return True
        if not "timestamp" in D:
            if self.warnedBefore == False:
                self.pprint ( f"'--before' is set to {self.before} but no timestamps in dict file" )
                self.warnedBefore = True
            return True
        from datetime import datetime as dt
        deadline = dt.strptime ( self.before, "%Y/%m/%d")
        current = dt.strptime ( D["timestamp"], "%Y/%m/%d")
        return deadline >= current

    def selectedSqrts( self, id ):
        """ select for sqrt-s """
        from smodels_utils.helper.various import getSqrts
        s = getSqrts ( id )
        if s in self.sqrts:
            return True
        return False

    def filterSigma ( self, ratio : float ) -> bool:
        """ filter according to self.filtersigma
        :param ratio: sqrt(expectedBG) / bgError
        :returns: true: if we wish to keep the value, else false
        """
        if self.filternegativesigma != None:
            return ratio <= self.filternegativesigma
        if self.filtersigma >= 0.:
            return ratio >= self.filtersigma
        return ratio <= self.filtersigma

    def read ( self ):
        """ read in content of self.filenames """
        from multiverse.expResModifier import readDatabaseDictFile
        for fname in self.filenames:
            ret = readDatabaseDictFile ( fname )
            self.meta.update (  ret["meta"] )
            newdata = {}
            for i,v in ret["data"].items():
                if self.nosuperseded and 'superseded' in v and v['superseded']:
                    continue
                if self.nofastlim and 'fastlim' in v and v['fastlim']:
                    continue
                if not self.selectedCollaboration ( i ):
                    continue
                if not self.selectedSqrts ( i ):
                    continue
                if not self.filterByTime ( v ):
                    continue
                if "expectedBG" in v and v["expectedBG"]>=self.filter and \
                        self.filterSigma ( math.sqrt(v["expectedBG"]) / v["bgError"] ):
#                        v["expectedBG"]/v["bgError"]>=self.filtersigma:
                    newdata[i]=v
                else:
                    if ":ul" in i:
                        if self.useAlsoULMaps:
                            txname = i [ i.rfind(":")+1: ]
                            v["txns"] = txname
                            newdata[i]=v
                        else:
                            if self.verbose > 2:
                                self.pprint ( f"[plotDBDict] removing {i} (is an UL)" )
                    else:
                        eBG,bgerr=None,None
                        if "expectedBG" in v:
                            eBG = v["expectedBG"]
                            bgerr = v["bgError"]
            self.data[ret["basename"]] = newdata

    def getSqrts100 ( self, anaid : str, lumi : Union[int,float] ) -> str:
        """ get the sqrts of anaid plus > 100 fb^-1 lumi, as string 
        :param anaid: e.g. CMS-SUS-20-004
        :param lumi: lumi as number in 1/fb, e.g. 136
        :returns: e.g. '13_gt'
        """
        from smodels_utils.helper.various import getSqrts
        sqrts = getSqrts ( anaid )
        if sqrts < 10:
            return str(sqrts)
        if lumi>100:
            return f"{sqrts}_gt"
        return f"{sqrts}_lt"

    def countSRs ( self ):
        """ count the number of signal regions for each analysis,
            for later reweighting """
        self.srCounts = {}
        for filename in self.filenames:
            selfbase = os.path.basename ( filename ).replace(".dict","")
            for label,v in self.data[selfbase].items():
                p1 = label.find(":")
                anaid = label[:p1]
                sr = label[p1+1:]
                if not anaid in self.srCounts:
                    self.srCounts[anaid]=set()
                self.srCounts[anaid].add ( sr )

    def isSelected ( self, txns : List ) -> bool:
        for tx in txns:
            if tx in self.select_topologies:
                return True
        return False

    def compute ( self ):
        """ compute the p-values """
        empty = {"8":[], "13_lt":[], "13_gt":[] }
        P,Pfake,weights, weightsfake = cp ( empty ), cp ( empty ), cp ( empty ), cp ( empty )
        self.countSRs()
        hasComplained = False
        nSkipped = 0
        for filename in self.filenames:
            selfbase = os.path.basename ( filename )
            dname = selfbase.replace(".dict","")
            data = self.data [ dname ]
            skipped = []
            self.nanas = set()
            hasEffMaps = set()
            for k,v in data.items():
                if ":ul" in k:
                    continue
                p1 = k.find(":")
                anaid = k[:p1]
                hasEffMaps.add ( anaid )
            for k,v in data.items():
                p1 = k.find(":")
                anaid = k[:p1]
                passesAnas = False
                if len(self.analyses)==0 and len(self.negativeanalyses)==0:
                    passesAnas=True
                for ana in self.analyses:
                    if fnmatch.fnmatch ( anaid, ana ):
                        passesAnas=True
                        break
                if len(self.negativeanalyses) != 0:
                    passesAnas=True
                    for ana in self.negativeanalyses:
                        # print ( f"@@b does {anaid} match {ana}: {fnmatch.fnmatch ( anaid, ana )}" )
                        if fnmatch.fnmatch ( anaid, ana ):
                            passesAnas=False
                            break
                if not passesAnas:
                    if not anaid in skipped:
                        if nSkipped < 4:
                            self.pprint ( f"skipping {anaid} per request" )
                        if nSkipped == 4:
                            self.pprint ( f"(quenching more msgs re skipping)" )
                        nSkipped += 1
                    skipped.append ( anaid )
                    continue
                w = 1. / len(self.srCounts[anaid]) / len(self.filenames)
                txns = []
                if "txns" in v:
                    txns = v["txns"] # .split(",")
                passesTx=False
                if len(self.topologies)==0 and len(self.negative_topologies)==0:
                    passesTx=True
                for tx in self.topologies:
                    if tx in txns:
                        passesTx=True
                        break
                if len(self.negative_topologies) != 0:
                    passesTx=True
                    for tx in self.negative_topologies:
                        if tx in txns:
                            passesTx=False
                            break
                if not passesTx:
                    self.printMsgs["pass"]+=1
                    if self.printMsgs["pass"]<4:
                        self.pprint ( f"skipping {k}: does not pass txname filter" )
                    if self.printMsgs["pass"]==4:
                        self.pprint ( f"(quenching 'does not pass filter' msgs)" )
                    continue

                sqrts = self.getSqrts100 ( k, v["lumi"] )
                if self.ignore_sqrts:
                    sqrts = "13_gt"  # if we ignore sqrts, we treat all as 13_gt
                if self.isSelected ( txns ):
                    sqrts = "8"
                if self.use_custom_function:
                    if isSelectedCustomFunction ( anaid, txns ):
                        sqrts = "8"
                if ":ul" in k:
                    if self.useAlsoULMaps and anaid in hasEffMaps:
                        self.pprint ( f"[plotDBDict] skipping {anaid}:ul: has effmaps." )
                    if self.useAlsoULMaps and not anaid in hasEffMaps:
                        # lets take the upper limit results with us
                        p = scipy.stats.norm.cdf( v["x"] )
                        P[sqrts].append (p )
                        w = 1. / len(self.filenames)
                        weights[sqrts].append ( w )
                else:
                    obs = v["origN"]
                    # obs = v["newObs"]
                    fakeobs = float("nan")
                    if "newObs" in v:
                        fakeobs = v["newObs"]
                    vexp = v["expectedBG"]
                    fudge = 1.
                    if "fudge" in v:
                        fudge = v["fudge"]
                        self.fudge = fudge
                    bgErr = v["bgError"]/fudge
                    if self.unscale:
                        bgErr = v["bgError"]
                    if vexp < self.filter:
                        continue
                    if not self.filterSigma ( math.sqrt(vexp) / bgErr ):
                    # if vexp / bgErr < self.filtersigma:
                        continue
                    sigN = None
                    if "sigN" in v:
                        sigN = v["sigN"]
                    # bgErr = v["bgError"]# /v["fudge"]
                    if abs(v["fudge"]-1)<1e-3:
                        v["orig_p_fudged"]=v["orig_p"]
                    if "orig_p_fudged" in v and self.likelihood == "gauss+poisson":
                        p = v["orig_p_fudged"]
                    else:
                        if not hasComplained:
                            self.pprint ( "computing the p-values -- this might take a while, so consider doing this at expResModifier.py" )
                            hasComplained = True
                        lognormal = False
                        if self.likelihood == "lognormal+poisson":
                            lognormal = True
                        p = computeP ( obs, vexp, bgErr )
                    if self.use_aggregated:
                        if anaid+"-agg" in hasEffMaps:
                            if not anaid in self.skippedAgg:
                                self.skippedAgg.add ( anaid )
                                self.pprint ( f"skipping {len(self.srCounts[anaid])} SRs in {anaid}: we also have aggregated results for this analysis" )
                            continue
                    else:
                        if "-agg" in anaid:
                            nonaggid = anaid.replace("-agg","")
                            checkIfNonAgg = nonaggid in hasEffMaps
                            if checkIfNonAgg:
                                if not nonaggid in self.skippedAgg:
                                    self.skippedAgg.add ( nonaggid )
                                    self.pprint ( f"skipping {len(self.srCounts[anaid])} SRs in {anaid}: we also have non-aggregated results for this analysis" )
                                continue
                    if self.verbose > 5:
                        nds = sum ( [len(x) for x in  P.values() ] )
                        print ( f"[plotDBDict] adding #{nds+1}: {k}" )
                    P[sqrts].append( p )
                    weights[sqrts].append ( w )

                    pfake = float("nan")
                    if "new_p" in v:
                        pfake = v["new_p"]
                    else:
                        if not math.isnan ( fakeobs):
                            pfake = computeP ( fakeobs, vexp, bgErr,
                                               lognormal = lognormal )
                    if not math.isnan ( pfake):
                        Pfake[sqrts].append( pfake )
                        weightsfake[sqrts].append ( w )
                self.nanas.add ( removeAnaIdSuffices ( anaid ) )
        for s in P.keys():
            P[s]=np.array(P[s])
            Pfake[s]=np.array(Pfake[s])
            weights[s]=np.array(weights[s])
            weightsfake[s]=np.array(weightsfake[s])
        return P,Pfake,weights,weightsfake

    def discussPs ( self, P, Pfake, weights, weightsfake ):
        Ptot = np.concatenate ( [ P["8"], P["13_lt"], P["13_gt"] ] )
        Pfaketot = np.concatenate ( [ Pfake["8"], Pfake["13_lt"], Pfake["13_gt"] ] )
        self.pprint ( f"real Ps: {len(Ptot)} entries at {np.mean(Ptot):.3f} +/- {np.std(Ptot):.2f}" )
        self.pprint ( f"fake Ps: {len(Pfaketot)} entries at {np.mean(Pfaketot):.3f} +/- {np.std(Pfaketot):.2f}" )
        for i in [ "8", "13_lt", "13_gt" ]:
            w, v = self.computeWeightedMean ( P[i], weights[i] )
            n = len(P[i])
            if n > 0:
                self.pprint ( f"real Ps, {i}: {n} entries at {w:.3f} +/- {v:.2f}" )

    def computeWeightedMean ( self, ps, ws ):
        """ weighted average of p values
        :param ps: array of p values
        :param ws: array of weights
        """
        if len(ps)==0:
            return 0., 0.
        Pi = ps*ws
        wtot = sum(ws)
        central = float ( np.sum(Pi) / wtot )
        # var = np.sum ( ws*ws*ps ) / wtot**2
        var = math.sqrt ( 1. / ( 12. * len(Pi) ) )
        return central, var

    def determineOutFile ( self, outfile = None ):
        """ determine the actual output file name, i.e.
            plug in for the @@FILTER@@ placeholders """
        if outfile is None:
            if "outfile" in self.options:
                outfile = self.options["outfile"]
        if outfile is None:
            return "tmp.png"
        """
        origt = self.origtopos.replace(" ","").replace(",","_")
        flt = f"_{origt}{'_not'.join(self.negative_topologies)}"
        flt += "_".join(self.analyses) # 
        if len(self.negativeanalyses)>0:
            flt += "_not"
            flt += "_not".join(self.negativeanalyses )
        flt = flt.replace("*","star").replace("?","questionmark")
        """
        if self.description == None:
            flt = ""
        else:
            flt = self.description.replace(" ","_")
        outfile = outfile.replace("@@FILTER@@", flt )
        return outfile

    def toSignificance ( self, p ):
        """ translate a p-value to a significane, i.e. compute Phi(p)^-1 """
        if type(p) == dict:
            ret = {}
            for k,v in p.items():
                ret[k] = self.toSignificance ( v )
            return ret
        if type(p) in [ list, tuple, np.array, np.ndarray ]:
            ret = []
            for k in p:
                ret.append ( self.toSignificance ( k ) )
            if type(p) == tuple:
                ret = tuple(ret)
            if type(p) in [ np.array, np.ndarray ]:
                ret = np.array ( ret )
            return ret
        if type(p) in [ float, np.float32, np.float64 ]:
            if p == 0.:
                return -10 # big number
            Z = - scipy.stats.norm.ppf ( p )
            if not np.isfinite ( Z ):
                newZ = np.sign ( Z ) * 6
                self.pprint ( f"significance for p={p} was {Z} will cap to {int(newZ)}!" )
                Z = newZ
            #if abs(Z) > 2.5:
            #    self.pprint ( "@@2 Z", Z, "p", p )
            return Z
        self.pprint ( f"cannot compute significance for {p} {type(p)}" )
        sys.exit()
        return None
        # scipy.stats.norm.ppf

    def interactive ( self, container ):
        import IPython
        IPython.embed( colors = "neutral" )

    def getTitle ( self ):
        """ determine the plot title """
        dbname = os.path.basename ( self.meta["database"] )
        title = f"SModelS database v{dbname}"
        if "+" in dbname: ## too long
            title = f"v{dbname}"
        # title = f"$p$-values, SModelS database v{dbname}"
        fudge = 1.
        if "fudge" in self.meta:
            fudge = self.meta["fudge"]
        #if abs ( fudge - 1. ) > 1e-3:
        #    title += f", fudge={fudge:.2f}"
        selecting = "selecting "
        if self.description != None:
            self.pprint ( f"we selected {','.join(self.topologies)}" )
            title += f", {self.description}"
        if len (self.topologies )>0 and self.description == None:
            stopos = ""
            for i,t in enumerate(self.topologies):
                if "+" in t and not "+off" in t:
                    self.pprint ( f"WARNING: topology {t} has a + sign, did you mean to instead have a comma ','?" )
                stopos += prettyDescriptions.prettyTxname( t, False, "latex" )
                if i < len(self.topologies)-1:
                    stopos += ";"
            title += f", {selecting}{stopos}"
            selecting = ""
        if len (self.negative_topologies )>0 and self.description == None:
            stopos = ""
            for i,t in enumerate(self.negative_topologies):
                stopos += f"^{prettyDescriptions.prettyTxname(t, False, 'latex')}"
                if i < len(self.topologies)-1:
                    stopos += ";"
            title += f", {selecting}{stopos}"
            selecting = ""
        if len ( self.topologies ) + len ( self.negative_topologies ) == 0:
            title += f", all topologies"
        if len ( self.analyses ) > 0:
            title += f", {selecting}"
            for a in self.analyses:
                title += f" {a}"
            title = title.replace("  ", " " )
        if len ( self.negativeanalyses ) > 0:
            for a in self.negativeanalyses:
                title += f", {selecting}^{a}"
                selecting = ""
        if len ( self.analyses ) + len ( self.negativeanalyses ) == 0:
            title += ", all analyses"
        if self.unscale:
            title += f" (unscaling)"
        if self.signalmodel:
            title += f" (signalmodel)"
        if self.title != None:
            title = self.title
        self.title = title
        return title

    def getBinNr ( self, bins : np.array, x : float ) -> int:
        """ given a histogram with edges at bins,
        find index of entry <x> """
        for i,b in enumerate(bins):
            if x < b:
                return i-1
        return len(bins)-2 # to the right of the last one
        #ret=int(x*len(bins)) ## find the bin of the max
        #self.pprint ( "bins", bins, "binnr", ret )
        #return ret

    def plot( self ):
        """ plot the p-values / significances """
        P,Pfake,weights,weightsfake=self.compute ( )
        if not self.pvalues:
            P,Pfake=self.toSignificance((P,Pfake))
        weighted = False
        if "weighted" in self.options:
            weighted = self.options["weighted"]
        if not "database" in self.meta:
            if "orig_dbpath" in self.meta:
                self.meta["database"]=self.meta["orig_dbpath"]
            else:
                self.pprint ( "error: database not defined in meta. did you pick up any dict files at all?" )
                sys.exit()
        title = self.getTitle()

        fig, ax = plt.subplots()
        x = [ P["8"], P["13_lt"], P["13_gt"] ]
        step, bins = self.getBins()

        avgp8,varp8 =self.computeWeightedMean ( P["8"], weights["8"] )
        bin8=self.getBinNr ( bins, avgp8 ) ## find the bin of the max
        avgp13lt, var13lt = self.computeWeightedMean( P["13_lt"], weights["13_lt"] )
        avgp13gt, var13gt = self.computeWeightedMean( P["13_gt"], weights["13_gt"] )
        bin13lt=self.getBinNr ( bins, avgp13lt )
        bin13gt=self.getBinNr ( bins, avgp13gt )
        nm1 = 1. / len(self.filenames)
        wlist = [ [nm1]*len(weights["8"]), [nm1]*len(weights["13_lt"]), [nm1]*len(weights["13_gt"]) ]
        if weighted:
            wlist = [ weights["8"], weights["13_lt"], weights["13_gt"] ]
        nontrivial = [ len(x)>0 for x in wlist ]
        savgp8 = f"{avgp8:.2f}".lstrip('0').replace("-0","-")
        savgp13l = f"{avgp13lt:.2f}".lstrip('0').replace("-0","-")
        savgp13g = f"{avgp13gt:.2f}".lstrip('0').replace("-0","-")
        labels = [ "8 TeV", "13 TeV, $\\mathcal{L}<78/fb$", "13 TeV, full $\\mathcal{L}$" ]
        plotAverages = False
        if "plot_averages" in self.options:
            plotAverages = self.options["plot_averages"]
        if plotAverages:
            labels = [ f"8 TeV [{savgp8}]", f"13 TeV, $\\mathcal{{L}}<100/fb$ [{savgp13l}]", f"13 TeV, $\\mathcal{{L}}>100/fb$ [{savgp13g}]" ]
        if self.ignore_sqrts:
            rest_text = "all other searches"
            labels = [ self.select_text, "--", rest_text ]
        nLegendEntries=0
        for c,l in enumerate(labels):
            if not nontrivial[c]:
                labels[c]=""
            else:
                nLegendEntries+=1
        # colors = [ "tab:green", "tab:blue", "cyan" ]
        # colors = [ "tab:green", "tab:blue", "lightblue" ]
        colors = [ "tab:green", "lightblue", "tab:blue" ]
        for i in [0,1,2]:
            if f"color{i}" in self.options:
                colors[i] = self.options[f"color{i}"]

        H1 = plt.hist ( x, weights = wlist, bins=bins, histtype="bar",
                   label= labels, color= colors, stacked=True )
        if "yrange" in self.options and self.options["yrange"]!=None:
            ax = plt.gca()
            ax.set_ylim(self.options["yrange"])
        mx = max ( H1[0][2] ) ## highest y-value, like at all
        # eps = .2
        eps = mx / 50.
        l8 = 0. + eps
        h8 = H1[0][0][bin8] - eps
        h13lt = H1[0][1][bin13lt] - eps
        l13lt = H1[0][0][bin13lt] + eps
        if l13lt > h13lt:
            l13lt, h13lt = h13lt, l13lt
        h13gt = H1[0][2][bin13gt] - eps
        l13gt = H1[0][1][bin13gt] + eps
        if l13gt > h13gt:
            l13gt, h13gt = h13gt, l13gt

        if plotAverages:
            if 8 in self.sqrts and ( avgp8 > 0. or not self.pvalues):
                l81 = plt.plot ( [ avgp8, avgp8 ], [l8, h8 ], color = "darkgreen", zorder=1, label = r"averages of $p$-values, $\bar{p}$", linewidth=2 )
                l82 = plt.plot ( [ avgp8+varp8, avgp8+varp8 ], [l8, h8 ], color = "darkgreen", zorder=1, linestyle="dotted", linewidth=1 )
                l83 = plt.plot ( [ avgp8-varp8, avgp8-varp8 ], [l8, h8 ], color = "darkgreen", zorder=1, linestyle="dotted", linewidth=1 )
            if 13 in self.sqrts and ( avgp13lt > 0. or not self.pvalues ):
                l13l = plt.plot ( [ avgp13lt, avgp13lt ], [ l13lt, h13lt ], color = "darkblue", zorder=1, linewidth=2 )
                l13l2 = plt.plot ( [ avgp13lt+var13lt, avgp13lt+var13lt ], [ l13lt, h13lt ], color = "darkblue", zorder=1, linestyle="dotted", linewidth=1 )
                l13l3 = plt.plot ( [ avgp13lt-var13lt, avgp13lt-var13lt ], [ l13lt, h13lt ], color = "darkblue", zorder=1, linestyle="dotted", linewidth=1 )

            if 13 in self.sqrts and ( avgp13gt > 0. or not self.pvalues ):
                l13gt1 = plt.plot ( [ avgp13gt, avgp13gt ], [ l13gt, h13gt ], color = "darkblue", zorder=1, linewidth=2 )
                l13gt2 = plt.plot ( [ avgp13gt+var13gt, avgp13gt+var13gt ], [ l13gt, h13gt ], color = "darkblue", zorder=1, linestyle="dotted", linewidth=1 )
                l13gt3 = plt.plot ( [ avgp13gt-var13gt, avgp13gt-var13gt ], [ l13gt, h13gt ], color = "darkblue", zorder=1, linestyle="dotted", linewidth=1 )
            if self.fakes:
                fweights = np.concatenate ( [ [nm1]*len(weights["8"]), [nm1]*len(weights["13_lt"]), [nm1]*len(weights["13_gt"]) ] )
                if weighted:
                    fweights = np.concatenate ( [ weightsfake["8"], weightsfake["13_lt"], weightsfake["13_gt"] ] )
            # fweights = [ [ nm1 ]*len(Pfake[8]), [ nm1 ]*len(Pfake[13]) ]
                linewidth = 3
                if not self.pvalues:
                    linewidth = 2
                H2 = plt.hist ( np.concatenate ( [ Pfake["8"], Pfake["13_lt"],
                        Pfake["13_gt"] ] ), weights = fweights, bins=bins,
                        stacked=True, zorder=9, label="synthetic SM-only data",
                        color=["red" ], linewidth=linewidth, histtype="step" )
        self.discussPs ( P, Pfake, weights, weightsfake )
        loc, bbox_to_anchor = "best", None
        _, stdnmx = list (self.getBins ( 100 ) )
        nmcolor = "red"
        nmcolor = "black"
        if self.pvalues:
            loc = "upper right"
            bbox_to_anchor = (1.12,1.02)
            if self.options["draw_reference"]:
                Ptot = float(sum(np.concatenate ( [ Pfake["8"], Pfake["13_lt"],
                                          Pfake["13_gt"] ] )) )
                ex = Ptot / self.nbins
                plt.plot ( [0,1], [ex,ex], c=nmcolor, linestyle="dotted",
                           label="SM hypothesis" )

        else:
            if self.options["draw_reference"]:
                scale = 1. / 0.39894 * .75
                stdnmy = [ scipy.stats.norm.pdf(x)*mx * scale for x in stdnmx ]
                plt.plot ( stdnmx, stdnmy, c=nmcolor, linestyle="dotted",
                           label="standard normal" )
        if nLegendEntries > 1 or self.options["alwayslegend"]:
            legend = plt.legend( loc = loc, facecolor=(1, 1, 1, 0.2),
                    bbox_to_anchor = bbox_to_anchor )
            legend.set_zorder ( 10 )
        if self.likelihood == "lognormal+poisson":
            title += " (lognormal)"
        if self.likelihood == "gauss":
            title += " (simple)"
        if self.collaboration != "ALL":
            title += f", {self.collaboration} only"
        title = title.replace("\\n","\n")
        title = title.replace("<<newline>>","\n")
        plt.title  ( title )
        if self.pvalues:
            plt.plot ( [ .5, .5 ], [ -.003, .2 ], c="tab:grey", linewidth=1,
                       linestyle="-" )
        xlabel  = "$p$-values"
        ylabel = "# SRs"
        if not self.pvalues:
            xlabel = "significances"
        if weighted:
            ylabel = "#analyses (weighted)"
        if "ylabel" in self.options:
            ylabel = self.options["ylabel"]
        if "xlabel" in self.options:
            xlabel = self.options["xlabel"]
        plt.xlabel ( xlabel )
        plt.ylabel ( ylabel )
        Ptot = np.concatenate ( [ P["8"], P["13_lt"], P["13_gt"] ] )
        nAnas = len ( self.nanas )
        if False:
            nanas = list ( self.nanas )
            nanas.sort()
            self.pprint ( nanas )
        nSRs = int ( len(Ptot) / len(self.filenames ) )
        plotStats = True
        if "plotStats" in self.options:
            plotStats = self.options["plotStats"]
        if plotStats:
            plt.text ( .67, -.12, f"this plot contains {nSRs} SRs from {nAnas} analyses", transform=ax.transAxes, c="grey", fontsize=7 )
        if abs ( self.fudge - 1. ) > 1e-5:
            plt.text ( -.1, -.12, f"fudge={self.fudge:.2f}", transform=ax.transAxes,
                       c="black", fontsize=10 )

        # plt.ylabel ( "# Signal Regions" )
        self.pprint ( f"plotting {GREEN}{self.outfile}{RESET}" )
        if self.comment != None:
            plt.text ( .65, -.11, self.comment, transform=ax.transAxes,
                       style="italic" )
        if self.disclaimer not in [ False, None, "None", "False" ]:
            plt.text ( .3, .3, self.disclaimer, transform=ax.transAxes,
                       rotation=35, c="#ff3333", fontsize=20 )
        from installation import version as protomodels_version
        metadata = { "protomodels_version": protomodels_version() }
        plt.kittyPlot ( self.outfile, self.show, metadata )

        plt.clf()
        plt.close()

def getArgs( cmdline = None ):
    import argparse
    argparser = argparse.ArgumentParser(description="meta statistics plotter, i.e. the thing that plots pDatabase.png")
    argparser.add_argument ( '-d', '--dictfile', nargs='*',
            help='input dictionary file(s) or directory, as generated eg via "expResModifier.py -d <smodels-database> -C" [./dicts/]',
            type=str, default='./dicts/' )
    argparser.add_argument ( '-o', '--outfile', nargs='?',
            help='output file [./pDatabase@@FILTER@@.png]',
            type=str, default='./pDatabase@@FILTER@@.png' )
    argparser.add_argument ( '-c', '--comment', nargs='?',
            help='an optional comment, to put in the plot [None]',
            type=str, default=None )
    argparser.add_argument ( '-u', '--unscale',
            help='unscale, i.e. use the fudged bgError also for computing likelihoods', action='store_true' )
    argparser.add_argument ( '--nosuperseded',
            help='ignore entries in the .dict file that are marked as superseded',
            action='store_true' )
    argparser.add_argument ( '--use_aggregated',
            help='instead of the non-aggregated signal regions, add the aggregated ones of a given result',
            action='store_true' )
    argparser.add_argument ( '--use_custom_function',
            help='select based on isSelectedCustomFunction',
            action='store_true' )
    argparser.add_argument ( '--nofastlim',
            help='ignore entries in the .dict file that are marked as fastlim',
            action='store_true' )
    argparser.add_argument ( '-w', '--weighted',
            help='weighted plot, i.e. each analysis (not each SR) counts equally', action='store_true' )
    argparser.add_argument ( '-F', '--fakes',
            help='add the fakes to the plot', action='store_true' )
    argparser.add_argument ( '--draw_reference',
            help='draw the reference distribution', action='store_true' )
    argparser.add_argument ( '-p', '--pvalues',
            help='plot p-values, not significances', action='store_true' )
    argparser.add_argument ( '-b', '--before',
            help='plot only entries before a certain date, like 2017/2/27',
            type=str, default=None )
    argparser.add_argument ( '-v', '--verbose',
            help='verbosity level [1]',
            type=int, default=1 )
    argparser.add_argument ( '-S', '--signalmodel',
            help='use the signal+bg model for computing likelihoods', action='store_true' )
    argparser.add_argument ( '-l', '--likelihood', nargs='?',
            help='likelihood: gauss (g), gauss+poisson (gp), or lognormal+poisson (lp) [gauss+poisson]',
            type=str, default="gauss+poisson" )
    argparser.add_argument ( '-t', '--topologies', nargs='?',
            help='filter for certain topologies, e.g. T1, T2tt. Comma separated. The signal region must have a map for any one of the given topologies. "^" before the name acts as negation [None]',
            type=str, default=None )
    argparser.add_argument ( '--select_topologies', nargs='?',
            help='filter for certain topologies to hilight, e.g. T1, T2tt. Comma separated. The signal region must have a map for any one of the given topologies. [None]',
            type=str, default=None )
    argparser.add_argument ( '--select_text', nargs='?',
            help='the text to go in the legend for the selected topos [selected]',
            type=str, default="selected" )
    argparser.add_argument ( '--sqrts', nargs='*',
            help='sqrtses [8,13,13.6]', type=float, default=[8,13,13.6] )
    argparser.add_argument ( '-a', '--analyses', nargs='?',
            help='filter for certain analyses, e.g. CMS-SUS-16-039-*. Unix-type wildcards. Comma separated. "^" before the name acts as negation [None]',
            type=str, default=None )
    argparser.add_argument ( '-f', '--filter', nargs='?',
            help='filter out signal regions with expectedBG<x [x=0.]',
            type=float, default=0. )
    argparser.add_argument ( '-s', '--filtersigma', nargs='?',
            help='filter out (remove) systematics dominated signal regions with sqrt(expectedBG)/bgErr<x. [x=0.]',
            type=float, default=0. )
    argparser.add_argument ( '-ns', '--filternegativesigma', nargs='?',
            help='filter out (remove) statistics dominated signal regions with sqrt(expectedBG)/bgErr>x. [x=None]',
            type=float, default=None )
    argparser.add_argument ( '-C', '--select_collaboration', nargs='?',
            help='select a specific collaboration CMS, ATLAS, all [all]',
            type=str, default="all" )
    argparser.add_argument ( '-T', '--title', nargs='?',
            help='supply an alternative title [None]',
            type=str, default=None )
    argparser.add_argument ( '-D', '--disclaimer',
            type=str, default=None )
    argparser.add_argument ( '--ignore_sqrts',
            help='plot results from all runs with the same color', action='store_true' )
    argparser.add_argument ( '-O', '--options',
            help='dictionary of options, given as string {try nbins, xlabel, ylabel, plotStats, plot_averages, weighted, yrange, color0, color1, color2} [None]',
            type=str, default=None )
    argparser.add_argument ( '-U', '--ulalso',
            help='upper limit results also (but also if not eff maps exist for a given analysis)', action='store_true' )
    argparser.add_argument ( '--list_abbreviations',
            help='list all abbreviations of topology names', action='store_true' )
    argparser.add_argument ( '--show',
            help='show plot', action='store_true' )
    argparser.add_argument ( '--Zmax',
            help='maximum Z signifances to plot (|Z|) [None]',
            type=float, default=None )
    if type(cmdline) in [ str ]:
        cmdline = cmdline.split()
        if "plotDBDict.py" in cmdline[0]:
            cmdline = cmdline[1:]

    args=argparser.parse_args( cmdline )
    if args.list_abbreviations:
        print ( )
        print ( f"{RED}Defined abbreviations:{RESET}" )
        print ( f"{RED}======================{RESET}" )
        shorts, descriptions = namesForSetsOfTopologies ( "list" )
        for short,topos in shorts.items():
            topos, descr = namesForSetsOfTopologies ( short )
            print ( f"{GREEN}{short}{RESET}: {topos}" )
            print ( f"         {YELLOW}''{descr}''{RESET}" )
        print ( )
        sys.exit()
        
    if type(args.options) == str:
        args.options = eval ( args.options )
    if args.options is None:
        args.options = {}
    for k,v in args.__dict__.items():
        args.options[k]=v
    return args.__dict__

def main():
    args = getArgs()
    plotter = Plotter ( args )

def runNotebook( cmdline, options = {} ):
    """ meant to be run from with a jupyter notebook
    :param cmdline: the command line arguments, e.g "-d ./db310.dict  -r"
    :param options: additional options
    :returns: plotter object
    """
    args = getArgs( cmdline )
    plotter = Plotter ( args )

    ret = plotter.plot( args.outfile )
    return plotter

if __name__ == "__main__":
    main()
