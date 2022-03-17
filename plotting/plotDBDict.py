#!/usr/bin/env python3

""" plot the meta statistics of database.dict """

from smodels_utils.plotting import mpkitty as plt
from smodels_utils.helper import prettyDescriptions
#import matplotlib
#matplotlib.use('agg')
#from matplotlib import pyplot as plt
import numpy as np
import os, glob, sys, math
sys.path.insert(0,"../")
from ptools.helpers import computeP
from copy import deepcopy as cp
import scipy.stats
import matplotlib.mlab as mlab

class Plotter:
    def __init__ ( self, pathname, filtervalue: float, comment, likelihood: str,
                   topologies, unscale, signalmodel, filtersigma: float,
                   collaboration : str, doFakes : bool, analyses : str,
                   disclaimer : bool ):
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
        :param filtersigma: filter out signal regions with expectedBG/bgErr < filtersigma
        :param collaboration: select a specific collaboration
        :param doFakes: add fakes to the plot
        :param analyses: if not None, then filter for these analyses 
                         (e.g. CMS-SUS-16-039-ma5)
        :param disclaimer: add a disclaimer, "do not circulate"
        """
        collaboration = collaboration.upper()
        if collaboration in [ "", "*" ]:
            collaboration = "ALL"
        if not collaboration in [ "CMS", "ATLAS", "ALL" ]:
            print ( "[plotDBDict] error: collaboration must be either CMS, ATLAS, or ALL." )
            sys.exit(-1)
        self.collaboration = collaboration
        abbreviations = { "g": "gauss", "gp": "gauss+poisson", "lp": "lognormal+poisson" }
        if likelihood in abbreviations:
            likelihood = abbreviations[likelihood]
        if likelihood not in [ "gauss", "gauss+poisson", "lognormal+poisson" ]:
            print ( "error, likelihood is to be one of: gauss, gauss+poisson, lognormal+poisson" )
            sys.exit()
        self.likelihood = likelihood ## False: gauss, True: lognormal
        self.doFakes = doFakes
        self.unscale = unscale
        self.signalmodel = signalmodel
        self.topologies = []
        self.analyses = []
        self.disclaimer = disclaimer
        self.negativetopos = []
        self.negativeanalyses = []
        self.filtersigma = filtersigma
        self.verbose = 0
        self.useAlsoULMaps = False
        if topologies not in [ None, "" ]:
            topos = topologies.split(",")
            for t in topos:
                if t.startswith ( "^" ):
                    self.negativetopos.append ( t[1:] )
                else:   
                    self.topologies.append ( t )
        if analyses not in [ None, "" ]:
            analyses = analyses.split(",")
            for a in analyses:
                if a.startswith ( "^" ):
                    self.negativeanalyses.append ( a[1:] )
                else:   
                    self.analyses.append ( a )
        self.filenames = []
        if comment in [ "None", "", "none" ]:
            comment = None
        self.comment = comment
        for pname in pathname:
            if os.path.isdir ( pname ):
                pname = pname + "/db*dict"
            self.filenames += glob.glob ( pname )
        self.filter = filtervalue
        self.meta = {}
        self.data = {}
        self.read()

    def selectedCollaboration( self, anaid ):
        """ does anaid pass the collaboration selection? """
        if self.collaboration == "ALL":
            return True
        if self.collaboration in anaid:
            return True
        return False

    def read ( self ):
        """ read in content of filename """
        for fname in self.filenames:
            with open( fname,"rt") as f:
                tmp=f.readlines()
            lines = []
            for line in tmp:
                if line.startswith("#"):
                    continue
                lines.append ( line )
            basename = os.path.basename ( fname ).replace(".dict","")
            self.meta.update (  eval(lines[0]) )
            nan=float("nan")
            data = eval("\n".join(lines[1:]))
            newdata = {}
            for i,v in data.items():
                if not self.selectedCollaboration ( i ):
                    continue
                if "expectedBG" in v and v["expectedBG"]>=self.filter and \
                        v["expectedBG"]/v["bgError"]>=self.filtersigma:
                    newdata[i]=v
                else:
                    if ":ul" in i:
                        if self.useAlsoULMaps:
                            txname = i [ i.rfind(":")+1: ]
                            v["txns"] = txname
                            v["lumi"] = 136.1
                            newdata[i]=v
                        else:
                            print ( f"[plotDBDict] removing {basename}:{i} (is an UL)" )
                    else:
                        eBG,bgerr=None,None
                        if "expectedBG" in v:
                            eBG = v["expectedBG"]
                            bgerr = v["bgError"]
                        #print ( f"[plotDBDict] removing {basename}:{i} (eBG is {eBG}+-{bgerr})" )
            # print ( f"[plotDBDict] keeping {len(newdata)}/{len(data)} for {basename}" )
            self.data[basename] = newdata

    def getSqrts ( self, anaid ):
        """ get the sqrts of anaid """
        ret = 13
        t = anaid.replace("CMS-","").replace("ATLAS-","").replace("SUSY-","").\
                  replace("SUS-","").replace("PAS-","").replace("EXO-","").replace("CONF-","")
        t = t[:t.find("-")]
        t = int(t) % 2000
        if t < 15:
            ret = 8
        return ret

    def getSqrts100 ( self, anaid, lumi ):
        """ get the sqrts of anaid plus > 100 fb^-1 lumi, as string """
        ret = 13
        t = anaid.replace("CMS-","").replace("ATLAS-","").replace("SUSY-","").\
                  replace("SUS-","").replace("PAS-","").replace("EXO-","").replace("CONF-","")
        t = t[:t.find("-")]
        t = int(t) % 2000
        if t < 15:
            ret = "8"
        else:
            ret = "13"
            if lumi>100:
                ret += "_gt"
            else:
                ret += "_lt"
        return ret

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

    def compute ( self ):
        """ compute the p-values """
        empty = {"8":[], "13_lt":[], "13_gt":[] }
        P,Pfake,weights, weightsfake = cp ( empty ), cp ( empty ), cp ( empty ), cp ( empty )
        self.countSRs()
        hasComplained = False
        for filename in self.filenames:
            selfbase = os.path.basename ( filename )
            data = self.data [ selfbase.replace(".dict","") ]
            skipped = []
            self.nanas = set()
            for k,v in data.items():
                if "ul" in k:
                    print ( "k", k, v )
                p1 = k.find(":")
                anaid = k[:p1]
                passesAnas = False 
                if len(self.analyses)==0 and len(self.negativeanalyses)==0:
                    passesAnas=True
                for ana in self.analyses:
                    if ana in anaid:
                        passesAnas=True
                        break
                if len(self.negativeanalyses) != 0:
                    passesAnas=True
                    for ana in self.negativeanalyses:
                        if ana in anaid:
                            passesAnas=False
                            break
                if not passesAnas:
                    if not anaid in skipped:
                        print ( f"[plotDBDict] skipping {anaid} per request" )
                    skipped.append ( anaid )
                    continue 
                w = 1. / len(self.srCounts[anaid]) / len(self.filenames)
                txns = []
                if "txns" in v:
                    txns = v["txns"].split(",")
                passesTx=False
                if len(self.topologies)==0 and len(self.negativetopos)==0:
                    passesTx=True
                for tx in self.topologies:
                    if tx in txns:
                        passesTx=True
                        break
                if len(self.negativetopos) != 0:
                    passesTx=True
                    for tx in self.negativetopos:
                        if tx in txns:
                            passesTx=False
                            break
                if not passesTx:
                    if self.verbose > 1:
                        print ( f"[plotDBDict] skipping {k}: does not pass Tx filter" )
                    continue

                sqrts = self.getSqrts100 ( k, v["lumi"] )
                if ":ul" in k:
                    if self.useAlsoULMaps:
                        # lets take the upper limit results with us
                        p = scipy.stats.norm.cdf( v["x"] )
                        P[sqrts].append (p )
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
                    bgErr = v["bgError"]/fudge
                    if self.unscale:
                        bgErr = v["bgError"]
                    if vexp < self.filter:
                        continue
                    if vexp / bgErr < self.filtersigma:
                        continue
                    sigN = None
                    if "sigN" in v:
                        sigN = v["sigN"]
                    # bgErr = v["bgError"]# /v["fudge"]
                    if "orig_p" in v and self.likelihood == "gauss+poisson":
                        p = v["orig_p"]
                    else:
                        if not hasComplained:
                            print ( "computing the p-values -- this might take a while, so consider doing this at expResModifier.py" )
                            hasComplained = True
                        lognormal = False
                        if self.likelihood == "lognormal+poissohn":
                            lognormal = True
                        p = computeP ( obs, vexp, bgErr )
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
                self.nanas.add ( anaid )
        for s in P.keys():
            P[s]=np.array(P[s])
            Pfake[s]=np.array(Pfake[s])
            weights[s]=np.array(weights[s])
            weightsfake[s]=np.array(weightsfake[s])
        return P,Pfake,weights,weightsfake

    def discussPs ( self, P, Pfake, weights, weightsfake ):
        Ptot = np.concatenate ( [ P["8"], P["13_lt"], P["13_gt"] ] )
        Pfaketot = np.concatenate ( [ Pfake["8"], Pfake["13_lt"], Pfake["13_gt"] ] )
        print ( "[plotDBDict] real Ps: %d entries at %.3f +/- %.2f" % 
                ( len(Ptot), np.mean(Ptot), np.std(Ptot)  ) )
        print ( "[plotDBDict] fake Ps: %d entries at %.3f +/- %.2f" % 
                ( len(Pfaketot), np.mean(Pfaketot), np.std(Pfaketot) ) )
        for i in [ "8", "13_lt", "13_gt" ]:
            w, v = self.computeWeightedMean ( P[i], weights[i] )
            n = len(P[i])
            if n > 0:
                print ( "[plotDBDict] real Ps, %s: %d entries at %.3f +/- %.2f" % 
                        ( i, n, w, v ) )

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
        

    def plot( self, outfile ):
        """ plot the p-values """
        flt = "_".join(self.topologies)+"_^".join(self.negativetopos)
        flt += "_".join(self.analyses)+"_^".join(self.negativeanalyses )
        outfile = outfile.replace("@@FILTER@@", flt )
        P,Pfake,weights,weightsfake=self.compute ( )
        if not "database" in self.meta:
            print ( "error: database not defined in meta. did you pick up any dict files at all?" )
            sys.exit()
        dbname = os.path.basename ( self.meta["database"] )
        title = f"SModelS database v{dbname}"
        # title = f"$p$-values, SModelS database v{dbname}"
        fudge = 1.
        if "fudge" in self.meta:
            fudge = self.meta["fudge"]
        if abs ( fudge - 1. ) > 1e-3:
            title += ", fudge=%.2f" % fudge
        selecting = "selecting "
        if len (self.topologies )>0:
            stopos = ""
            for i,t in enumerate(self.topologies):
                stopos += prettyDescriptions.prettyTxname( t, "latex", False )
                if i < len(self.topologies)-1:
                    stopos += ";"
            title += f", {selecting}{stopos}"
            selecting = ""
        if len (self.negativetopos )>0:
            stopos = ""
            for i,t in enumerate(self.negativetopos):
                stopos += "^"+prettyDescriptions.prettyTxname( t, "latex", False )
                if i < len(self.topologies)-1:
                    stopos += ";"
            title += f", {selecting}{stopos}"
            selecting = ""
        if len ( self.topologies ) + len ( self.negativetopos ) == 0:
            title += f", all topologies"
        if len ( self.analyses ) > 0:
            for a in self.analyses:
                title += f", {selecting}{a}"
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
        nbins = 10 ## change the number of bins
        fig, ax = plt.subplots()
        x = [ P["8"], P["13_lt"], P["13_gt"] ]
        avgp8,varp8 =self.computeWeightedMean ( P["8"], weights["8"] )
        bin8=int(avgp8*nbins)
        avgp13lt, var13lt = self.computeWeightedMean( P["13_lt"], weights["13_lt"] ) 
        avgp13gt, var13gt = self.computeWeightedMean( P["13_gt"], weights["13_gt"] )
        bin13lt=int(avgp13lt*nbins)
        bin13gt=int(avgp13gt*nbins)
        nm1 = 1. / len(self.filenames)
        wlist = [ weights["8"], weights["13_lt"], weights["13_gt"] ]
        nontrivial = [ len(x)>0 for x in wlist ]
        bins = np.arange ( 0., 1+1e-7, 1/nbins )
        # labels = [ "real, 8 TeV", "real, 13 TeV", "real, 13 TeV, > 100 / fb" ]
        savgp8 = ( "%.2f" % avgp8 ).lstrip('0')
        savgp13l = ( "%.2f" % avgp13lt ).lstrip('0')
        savgp13g = ( "%.2f" % avgp13gt ).lstrip('0')
        labels = [ "8 TeV [%s]" % savgp8, "13 TeV, $\\mathcal{L}<100/fb$ [%s]" % savgp13l, "13 TeV, $\\mathcal{L}>100/fb$ [%s]" % savgp13g ]
        nLegendEntries=0
        for c,l in enumerate(labels):
            if not nontrivial[c]:
                labels[c]=""
            else:
                nLegendEntries+=1
        colors = [ "tab:green", "tab:blue", "cyan" ]
        H1 = plt.hist ( x, weights = wlist, bins=bins, histtype="bar",
                   label= labels, color= colors, stacked=True )
        mx = max ( H1[0][2] )
        #eps = .2
        eps = mx / 50.
        l8 = 0 + eps
        h8 = H1[0][0][bin8] - eps
        h13lt = H1[0][1][bin13lt] - eps
        l13lt = H1[0][0][bin13lt] + eps
        if l13lt > h13lt:
            l13lt, h13lt = h13lt, l13lt
        h13gt = H1[0][2][bin13gt] - eps
        l13gt = H1[0][1][bin13gt] + eps
        if l13gt > h13gt:
            l13gt, h13gt = h13gt, l13gt
        if avgp8 > 0.:
            l81 = plt.plot ( [ avgp8, avgp8 ], [l8, h8 ], color = "darkgreen", zorder=1, label = r"averages of $p$-values, $\bar{p}$", linewidth=2 )
            l82 = plt.plot ( [ avgp8+varp8, avgp8+varp8 ], [l8, h8 ], color = "darkgreen", zorder=1, linestyle="dotted", linewidth=1 )
            l83 = plt.plot ( [ avgp8-varp8, avgp8-varp8 ], [l8, h8 ], color = "darkgreen", zorder=1, linestyle="dotted", linewidth=1 )
        if avgp13lt > 0.:
            l13l = plt.plot ( [ avgp13lt, avgp13lt ], [ l13lt, h13lt ], color = "darkblue", zorder=1, linewidth=2 )
            l13l2 = plt.plot ( [ avgp13lt+var13lt, avgp13lt+var13lt ], [ l13lt, h13lt ], color = "darkblue", zorder=1, linestyle="dotted", linewidth=1 )
            l13l3 = plt.plot ( [ avgp13lt-var13lt, avgp13lt-var13lt ], [ l13lt, h13lt ], color = "darkblue", zorder=1, linestyle="dotted", linewidth=1 )

        if avgp13gt > 0.:
            l13gt1 = plt.plot ( [ avgp13gt, avgp13gt ], [ l13gt, h13gt ], color = "darkcyan", zorder=1, linewidth=2 )
            l13gt2 = plt.plot ( [ avgp13gt+var13gt, avgp13gt+var13gt ], [ l13gt, h13gt ], color = "darkcyan", zorder=1, linestyle="dotted", linewidth=1 )
            l13gt3 = plt.plot ( [ avgp13gt-var13gt, avgp13gt-var13gt ], [ l13gt, h13gt ], color = "darkcyan", zorder=1, linestyle="dotted", linewidth=1 )
        #fweights = [ nm1 ]*len(Pfake["8"]+Pfake["13_lt"]+Pfake["13_gt"])
        if self.doFakes:
            fweights = np.concatenate ( [ weightsfake["8"], weightsfake["13_lt"], weightsfake["13_gt"] ] )
        # fweights = [ [ nm1 ]*len(Pfake[8]), [ nm1 ]*len(Pfake[13]) ]
            H2 = plt.hist ( np.concatenate ( [ Pfake["8"], Pfake["13_lt"], Pfake["13_gt"] ] ), weights = fweights,
                        bins=bins, stacked=True, zorder=9,
                        label="fake", color=["red" ], linewidth=3, histtype="step" )
        self.discussPs ( P, Pfake, weights, weightsfake )
        # loc = "lower center"
        loc = "best"
        if nLegendEntries > 1:
            legend = plt.legend( loc = loc, facecolor=(1, 1, 1, 0.1) )
        if self.likelihood == "lognormal+poisson":
            title += " (lognormal)"
        if self.likelihood == "gauss":
            title += " (simple)"
        if self.collaboration != "ALL":
            title += f" {self.collaboration} only"
        plt.title  ( title )
        plt.plot ( [ .5, .5 ], [ -.003, .2 ], c="tab:grey", linewidth=1, 
                   linestyle="-" )
        plt.xlabel ( "$p$-values" )
        plt.ylabel ( "# analyses (weighted)" )
        Ptot = np.concatenate ( [ P["8"], P["13_lt"], P["13_gt"] ] )
        nAnas = len ( self.nanas )
        nSRs = len(Ptot)
        plt.text ( .69, -.12, f"this plot contains {nSRs} SRs from {nAnas} analyses", transform=ax.transAxes, c="black", fontsize=7 )
        # plt.ylabel ( "# Signal Regions" )
        print ( f"[plotDBDict] plotting {outfile}"  )
        if self.comment != None:
            plt.text ( .65, -.11, self.comment, transform=ax.transAxes, 
                       style="italic" )
        if self.disclaimer:
            plt.text ( .3, .3, "do not circulate!", transform=ax.transAxes,
                       rotation=35, c="#ff3333", fontsize=20 )
        plt.savefig ( outfile )
        if hasattr ( plt, "options" ) and plt.options["hasKittyBackend"]:
            plt.show()
        plt.clf()

def main():
    import argparse
    argparser = argparse.ArgumentParser(description="meta statistics plotter, i.e. the thing that plots pDatabase.png")
    argparser.add_argument ( '-d', '--dictfile', nargs='*',
            help='input dictionary file(s) or directory, as generated eg via "expResModifier.py -d <smodels-database> -C" [../data/database/]',
            type=str, default='../data/database/' )
    argparser.add_argument ( '-o', '--outfile', nargs='?',
            help='output file [./pDatabase@@FILTER@@.png]',
            type=str, default='./pDatabase@@FILTER@@.png' )
    argparser.add_argument ( '-c', '--comment', nargs='?',
            help='an optional comment, to put in the plot [None]',
            type=str, default=None )
    argparser.add_argument ( '-u', '--unscale', 
            help='unscale, i.e. use the fudged bgError also for computing likelihoods', action='store_true' )
    argparser.add_argument ( '-F', '--fakes', 
            help='add the fakes to the plot', action='store_true' )
    argparser.add_argument ( '-S', '--signalmodel', 
            help='use the signal+bg model for computing likelihoods', action='store_true' )
    argparser.add_argument ( '-l', '--likelihood', nargs='?',
            help='likelihood: gauss (g), gauss+poisson (gp), or lognormal+poisson (lp) [gauss+poisson]',
            type=str, default="gauss+poisson" )
    argparser.add_argument ( '-t', '--topologies', nargs='?',
            help='filter for certain topologies, e.g. T1, T2tt. Comma separated. The signal region must have a map for any one of the given topologies. "^" before the name acts as negation [None]',
            type=str, default=None )
    argparser.add_argument ( '-a', '--analyses', nargs='?',
            help='filter for certain analyses, e.g. CMS-SUS-16-039-ma5. Comma separated. "^" before the name acts as negation [None]',
            type=str, default=None )
    argparser.add_argument ( '-f', '--filter', nargs='?',
            help='filter out signal regions with expectedBG<x [x=0.]',
            type=float, default=0. )
    argparser.add_argument ( '-s', '--filtersigma', nargs='?',
            help='filter out signal regions with expectedBG/bgErr<x [x=0.]',
            type=float, default=0. )
    argparser.add_argument ( '-C', '--select_collaboration', nargs='?',
            help='select a specific collaboration CMS, ATLAS, all [all]',
            type=str, default="all" )
    argparser.add_argument ( '-D', '--disclaimer', 
            help='add a disclaimer', action='store_true' )
    args=argparser.parse_args()
    if args.topologies != None and args.topologies.endswith ( ".py" ):
        print ( f"[plotDBDict] you supplied {args.topologies} as topologies. Did you supply the validation file instead?" )
    plotter = Plotter ( args.dictfile, args.filter, args.comment, args.likelihood, 
                        args.topologies, args.unscale, args.signalmodel,
                        args.filtersigma, args.select_collaboration,
                        args.fakes, args.analyses, args.disclaimer )
    plotter.plot( args.outfile )

if __name__ == "__main__":
    main()
