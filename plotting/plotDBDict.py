#!/usr/bin/env python3

""" plot the meta statistics of database.dict """

from smodels_utils.plotting import mpkitty as plt
from smodels_utils.helper import prettyDescriptions

import numpy as np
import os, glob, sys, math
from copy import deepcopy as cp
import scipy.stats
import matplotlib.mlab as mlab
from typing import Union

sys.path.insert(0,"../")
from ptools.helpers import computeP
from ptools.moreHelpers import namesForSetsOfTopologies

class Plotter:

    def roughviz_template( self, data, labels, values, plot_svg, **kwargs):
        """ the template for the roughviz plot, we will overwrite the
            original one with this """
        from jinja2 import Template
        import random
        import string
        template = Template(data.decode("utf-8"))
        id_name = ''.join(random.choice(string.ascii_lowercase) for i in range(10))
        output = template.render(id_name = id_name,
                                 labels = labels,
                                 values = values,
                                 kwargs = kwargs)
        # print ( f"kwargs {kwargs}" )
        if(plot_svg):
            svg_id = "svg"+id_name
            script = """
            <style>
            div.output_area img, div.output_area svg{
            height: 100%!important;
            }
            </style>
            <script>
                var e = document.getElementById('"""+id_name+"""');
                var divCheckingInterval = setInterval(function(){
                if(e.getElementsByTagName('svg').length){
                    clearInterval(divCheckingInterval);
                    e.getElementsByTagName('svg')[0].setAttribute("id", '"""+svg_id+"""');
                    var svgElement = document.getElementById('"""+svg_id+"""');
                    var svgString = new XMLSerializer().serializeToString(svgElement);
                    var decoded = unescape(encodeURIComponent(svgString));
                    var base64 = btoa(decoded);
                    var imgSource = `data:image/svg+xml;base64,${base64}`;
                    e.innerHTML = "<img id='svgplot'>";
                    document.getElementById('svgplot').src = imgSource;       
                }}, 1);
            </script>
            """
            self.output = output
            self.script = script
            #self.display()
            self.saveRoughViz()
        else:
            self.output = output
            #self.display()

    def saveRoughViz ( self ):
        outfile = self.determineOutFile ( self.outfile )
        # htmlname = "plot.html"
        htmlname = outfile.replace(".png",".html" )
        f = open ( htmlname, "wt" )
        f.write ( self.output )
        if hasattr ( self, "script" ):
            f.write ( self.script )
        f.close()
        # self.pprint ( f"{htmlname} created, can be opened with ''xdg-open {htmlname}''" )
        cwd = os.getcwd()
        #pngname = "plot.png"
        pngname = outfile
        from shutil import which
        import subprocess
        cutycapt = which ( "cutycapt", path=f"/usr/bin:{os.environ['PATH']}" )
        if cutycapt is None:
            self.pprint ( f"cutycapt is not installed. maybe perform sudo apt install cutycapt" )
        else:
            cmd = f"cutycapt --zoom-factor=2. --delay=1000 --url=file://{cwd}/{htmlname} --out={pngname}"
            o = subprocess.getoutput ( cmd )
            self.pprint ( f"{pngname} created" )
            self.addLegendToRough ( pngname )
            self.showPng( pngname )
        if False:
            cmd = f"xdg-open {htmlname}"
            subprocess.getoutput ( cmd )
        if True and os.path.exists (htmlname ):
            # self.pprint ( f"removing {htmlname}" )
            os.unlink ( htmlname )
        
    def showPng ( self, filename ):
        if ("show" in self.options and self.options["show"]==True) or self.show == True:
            from smodels_utils.plotting.mpkitty import timg
            timg ( filename )

    def display ( self ):
        """ show html """
        from IPython.display import display, HTML
        display(HTML(self.output))
        if hasattr ( self, "script" ):
            display(HTML(self.script))

    def getBins ( self, nbins = None ):
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
        self.Zmax = 3.25
        self.before = None
        self.show = False
        self.pvalues = False # if False, then p-values if true then significances
        self.origtopos = "all"
        self.collaboration = "ALL"
        self.likelihood = "gauss+poisson"
        self.useAlsoULMaps = False
        self.analyses = []
        self.comment = None
        self.topologies = []
        self.negativetopos = []
        self.negativeanalyses = []
        self.outfile = self.determineOutFile ( "./pDatabase@@FILTER@@.png" )
        self.title = None
        self.roughviz = False
        self.options = { "alwayslegend": False }
        self.yrange = None
        self.unscale = False
        self.signalmodel = False
        self.fakes = False
        self.sqrts = [8,13,13.6]
        self.filter = 0.
        self.filtersigma = 0.
        self.disclaimer = False

    def __init__ ( self, args ):
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
        :param fakes: add fakes to the plot
        :param analyses: if not None, then filter for these analyses
                         (e.g. CMS-SUS-16-039-ma5)
        :param disclaimer: add a disclaimer, "do not circulate"
        :param ulAlso: show UL results, also
        :param title: a title
        :param pvalues: if true then plot p-values, if false plot significances
        """
        self.defaults()
        for a,value in args.items():
            if a=="options":
                self.options.update ( value )
                continue
            if a not in [ "topologies", "analyses" ]:
                setattr ( self, a, value )
        if self.nbins == None:
            if not self.pvalues:
                self.nbins = 13
            else:
                self.nbins = 10
        self.origtopos = args["topologies"]
        if self.origtopos == None:
            self.origtopos = "all"
        self.description = None
        if args['topologies'] != None:
                if args["topologies"].endswith ( ".py" ):
                    print ( f"[plotDBDict] you supplied {args['topologies']} as topologies. Did you supply the validation file instead?" )
                args["topologies"], descr = namesForSetsOfTopologies ( args['topologies'] )
                self.description = descr
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
                print ( "error, likelihood is to be one of: gauss, gauss+poisson, lognormal+poisson" )
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
                    self.negativetopos.append ( t[1:] )
                else:
                    self.topologies.append ( t )

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
                    pname = pname + "/db*dict"
                    self.filenames += glob.glob ( pname )
                else:
                    self.filenames.append ( pname )
        else:
            print ( "we need dictfile in args" )
        if "outfile" in args:
            self.outfile = self.determineOutFile ( args["outfile"] )
        self.meta = {}
        self.data = {}
        self.read()
        if self.roughviz:
            self.rough( )
        else:
            self.plot( )

    def pprint ( self, args, verbose = 0 ):
        # x = " ".join(map(str,args))
        x = args
        if verbose > self.verbose:
            print ( f"[plotDBDict] {x}" )

    def selectedCollaboration( self, anaid ):
        """ does anaid pass the collaboration selection? """
        if self.collaboration in [ "ALL", "all", "*" ]:
            return True
        if self.collaboration in anaid:
            return True
        return False

    def filterByTime ( self, D ):
        """ filter by time, let everything before self.before pass! """
        if self.before == None:
            return True
        if not "timestamp" in D:
            return True
        from datetime import datetime as dt
        deadline = dt.strptime ( self.before, "%Y/%m/%d")
        current = dt.strptime ( D["timestamp"], "%Y/%m/%d")
        # print ( "compare", self.before,"and",D["timestamp"], deadline >= current )
        return deadline >= current

    def selectedSqrts( self, id ):
        """ select for sqrt-s """
        from smodels_utils.helper.various import getSqrts
        s = getSqrts ( id )
        if s in self.sqrts:
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
                if not self.selectedSqrts ( i ):
                    continue
                if not self.filterByTime ( v ):
                    continue
                if "expectedBG" in v and v["expectedBG"]>=self.filter and \
                        v["expectedBG"]/v["bgError"]>=self.filtersigma:
                    newdata[i]=v
                else:
                    if ":ul" in i:
                        if self.useAlsoULMaps:
                            txname = i [ i.rfind(":")+1: ]
                            v["txns"] = txname
                            newdata[i]=v
                        else:
                            self.pprint ( f"removing {basename}:{i} (is an UL)", verbose = 2 )
                    else:
                        eBG,bgerr=None,None
                        if "expectedBG" in v:
                            eBG = v["expectedBG"]
                            bgerr = v["bgError"]
                        #print ( f"[plotDBDict] removing {basename}:{i} (eBG is {eBG}+-{bgerr})" )
            # print ( f"[plotDBDict] keeping {len(newdata)}/{len(data)} for {basename}" )
            self.data[basename] = newdata
            # print ( f"found {fname} {len(tmp)} basename >>{basename}<< newdata {len(newdata)}" )

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
                        self.pprint ( f"skipping {anaid} per request" )
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
                    self.pprint ( f"skipping {k}: does not pass Tx filter", verbose = 1 )
                    continue

                sqrts = self.getSqrts100 ( k, v["lumi"] )
                if ":ul" in k:
                    if self.useAlsoULMaps and anaid in hasEffMaps:
                        print ( f"[plotDBDict] skipping {anaid}:ul: has effmaps." )
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
        self.pprint ( "real Ps: %d entries at %.3f +/- %.2f" %
                ( len(Ptot), np.mean(Ptot), np.std(Ptot)  ), verbose = 1 )
        self.pprint ( "fake Ps: %d entries at %.3f +/- %.2f" %
                ( len(Pfaketot), np.mean(Pfaketot), np.std(Pfaketot) ),
                verbose = 1 )
        for i in [ "8", "13_lt", "13_gt" ]:
            w, v = self.computeWeightedMean ( P[i], weights[i] )
            n = len(P[i])
            if n > 0:
                self.pprint ( "real Ps, %s: %d entries at %.3f +/- %.2f" %
                        ( i, n, w, v ), verbose = 1 )

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
        origt = self.origtopos.replace(" ","").replace(",","_")
        flt = "_"+origt+"_^".join(self.negativetopos)
        flt += "_".join(self.analyses)+"_^".join(self.negativeanalyses )
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
            return Z
        print ( "cannot compute significance for",p,type(p) )
        sys.exit()
        return None
        # scipy.stats.norm.ppf

    def rough ( self ):
        """ roughviz plot of the same data """
        outfile = self.determineOutFile ( self.outfile )
        debug = []
        P,Pfake,weights,weightsfake=self.compute ( )
        if not self.pvalues:
            P,Pfake=self.toSignificance((P,Pfake))
        if not "database" in self.meta:
            print ( "error: database not defined in meta. did you pick up any dict files at all?" )
            sys.exit()
        title = self.getTitle()
        weighted = False
        if "weighted" in self.options:
            weighted = self.options["weighted"]
        import roughviz
        roughviz.roughviz.generate_template = self.roughviz_template
        if hasattr ( roughviz, "charts" ):
            print ( "I think you installed py-roughviz, not roughviz" )
            sys.exit(-1)
        import pandas as pd
        if not "database" in self.meta:
            print ( "error: database not defined in meta. did you pick up any dict files at all?" )
            sys.exit()
        title = self.getTitle()

        step, bins = self.getBins()
        # print ( f"[plotDBDict bins are at {bins}" )

        (p8,x8) = np.histogram ( P["8"], bins )
        (p13lt,x13lt) = np.histogram ( P["13_lt"], bins )
        (p13gt,x13gt) = np.histogram ( P["13_gt"], bins )
        factor = 1.
        if weighted:
            factor = 100.
            (p8,x8) = np.histogram ( P["8"], bins, weights=weights["8"] )
            (p13lt,x13lt) = np.histogram ( P["13_lt"], bins, weights=weights["13_lt"] )
            (p13gt,x13gt) = np.histogram ( P["13_gt"], bins, weights=weights["13_gt"] )
        sbins = [ f"{x+step/2.:.2f}" for x in x8[:-1] ]
        if not self.pvalues:
            sbins = [ f"{x+step/2.:.1f}" for x in x8[:-1] ]
        p8l = [ factor*float(x) for x in p8 ]
        p13ltl = [ factor*float(x) for x in p13lt ]
        #p13ltl = [ float(x)+float(y) for x,y in zip(p13lt,p8) ]
        p13gtl = [ factor*float(x) for x in p13gt ]
        #p13gtl = [ float(x)+float(y) for x,y in zip(p13gt,p13ltl) ]
        d = { "labels": sbins, "8 TeV": p8l, "13 TeV, < 100/fb": p13ltl, "13 TeV, > 100/fb": p13gtl }
        df = pd.DataFrame ( data = d )
        if "tilde" in title:
            title = f"${title}$"
        columns = [ "8 TeV", "13 TeV, < 100/fb", "13 TeV, > 100/fb" ]
        yLabel = "# SRs"
        if weighted:
            yLabel = "# analyses (weighted, x 100)"
        if "ylabel" in self.options:
            yLabel = self.options["ylabel"]
        roughness = 6
        if "roughness" in self.options:
            roughness = self.options["roughness"]
        if "title" in self.options:
            if self.options["title"] in [ False, None ]:
                # title = None
                pass
            if type(self.options["title"]) == str:
                title = self.options["title"]
        xlabel = "p-values"
        if not self.pvalues:
            xlabel = "significances"
        if "ylabel" in self.options:
            yLabel = self.options["ylabel"]
        if "xlabel" in self.options:
            xlabel = self.options["xlabel"]
        bar = roughviz.stackedbar ( df["labels"], df[ columns],
                xLabel=xlabel, roughness = roughness,
                yLabel = yLabel, title = title,
                titleFontSize = 18, plot_svg = True, interactive = False,
                labelFontSize = 16, axisFontSize = 16, legend = "true" )
        # bar = roughviz.outputs
        # self.interactive( { "df": df, "bar": bar, "debug": debug }  )
        # self.addLegendToRough ( outfile )
        return bar, debug

    def addLegendToRough ( self, filename ):
        """ rough plot does not have legend, so we write it ourselves """
        from PIL import Image, ImageDraw, ImageFont
        import ptools
        path = ptools.__file__.replace("ptools/__init__.py","shared/")
        font = os.path.join ( os.path.abspath ( path ), "Gaegu-Regular.ttf" )
        img = Image.open ( filename )
        d1 = ImageDraw.Draw(img)
        myFont = ImageFont.truetype( font, 36 )
        labels = { "8": "8 TeV", "13lt": "13 TeV, low lumi", "13gt": "13 TeV, high lumi" }
        colors = { "8": (135, 207, 236), "13lt": (121, 202, 176), "13gt": (218, 194, 161 ) }
        ymin, dy = 80, 50
        ycoords = { "8": ymin+2*dy, "13lt": ymin+dy, "13gt": ymin }
        for l in labels:
            txt = labels[l]
            c = colors[l]
            y = ycoords[l]
            d1.text((1100, y), txt, fill = c,font=myFont)
        img.save ( filename )

    def interactive ( self, container ):
        import IPython
        IPython.embed( colors = "neutral" )

    def getTitle ( self ):
        """ determine the plot title """
        dbname = os.path.basename ( self.meta["database"] )
        title = f"SModelS database v{dbname}"
        # title = f"$p$-values, SModelS database v{dbname}"
        fudge = 1.
        if "fudge" in self.meta:
            fudge = self.meta["fudge"]
        if abs ( fudge - 1. ) > 1e-3:
            title += ", fudge=%.2f" % fudge
        selecting = "selecting "
        if self.description != None:
            self.pprint ( f"we selected {','.join(self.topologies)}", verbose = 1 )
            title += f", {self.description}"
            # title += f",selecting {self.origtopos}"
        if len (self.topologies )>0 and self.description == None:
            stopos = ""
            for i,t in enumerate(self.topologies):
                if "+" in t and not "+off" in t:
                    print ( f"[plotDBDict] WARNING: topology {t} has a + sign, did you mean to instead have a comma ','?" )
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
            title += f", {selecting}"
            for a in self.analyses:
                title += f", {a}"
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
        #print ( "bins", bins, "binnr", ret )
        #return ret

    def plot( self ):
        """ plot the p-values """
        P,Pfake,weights,weightsfake=self.compute ( )
        if not self.pvalues:
            P,Pfake=self.toSignificance((P,Pfake))
        weighted = False
        if "weighted" in self.options:
            weighted = self.options["weighted"]
        if not "database" in self.meta:
            print ( "error: database not defined in meta. did you pick up any dict files at all?" )
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
        wlist = [ [1.]*len(weights["8"]), [1.]*len(weights["13_lt"]), [1.]*len(weights["13_gt"]) ]
        if weighted:
            wlist = [ weights["8"], weights["13_lt"], weights["13_gt"] ]
        nontrivial = [ len(x)>0 for x in wlist ]
        # labels = [ "real, 8 TeV", "real, 13 TeV", "real, 13 TeV, > 100 / fb" ]
        savgp8 = ( "%.2f" % avgp8 ).lstrip('0')
        savgp13l = ( "%.2f" % avgp13lt ).lstrip('0')
        savgp13g = ( "%.2f" % avgp13gt ).lstrip('0')
        # labels = [ "8 TeV", "13 TeV, $\\mathcal{L}<100/fb$", "13 TeV, $\\mathcal{L}>100/fb$" ]
        #labels = [ "8 TeV", "13 TeV, $\\mathcal{L}<78/fb$", "13 TeV, full lumi" ]
        labels = [ "8 TeV", "13 TeV, $\\mathcal{L}<78/fb$", "13 TeV, full $\\mathcal{L}$" ]
        plotAverages = True
        if "plot_averages" in self.options:
            plotAverages = self.options["plot_averages"]
        if plotAverages:
            labels = [ "8 TeV [%s]" % savgp8, "13 TeV, $\\mathcal{L}<100/fb$ [%s]" % savgp13l, "13 TeV, $\\mathcal{L}>100/fb$ [%s]" % savgp13g ]
        nLegendEntries=0
        for c,l in enumerate(labels):
            if not nontrivial[c]:
                labels[c]=""
            else:
                nLegendEntries+=1
        # colors = [ "tab:green", "tab:blue", "cyan" ]
        colors = [ "tab:green", "tab:blue", "lightblue" ]
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
                fweights = np.concatenate ( [ weightsfake["8"], weightsfake["13_lt"], weightsfake["13_gt"] ] )
            # fweights = [ [ nm1 ]*len(Pfake[8]), [ nm1 ]*len(Pfake[13]) ]
                H2 = plt.hist ( np.concatenate ( [ Pfake["8"], Pfake["13_lt"], Pfake["13_gt"] ] ), weights = fweights,
                            bins=bins, stacked=True, zorder=9,
                            label="fake", color=["red" ], linewidth=3, histtype="step" )
        self.discussPs ( P, Pfake, weights, weightsfake )
        # loc = "lower center"
        loc = "best"
        _, stdnmx = list (self.getBins ( 100 ) )
        scale = 1. / 0.39894 * .75
        stdnmy = [ scipy.stats.norm.pdf(x)*mx * scale for x in stdnmx ]
        if not self.pvalues:
            plt.plot ( stdnmx, stdnmy, c="red", linestyle="dotted", label="standard normal" )
        if nLegendEntries > 1 or self.options["alwayslegend"]:
            legend = plt.legend( loc = loc, facecolor=(1, 1, 1, 0.1) )
        if self.likelihood == "lognormal+poisson":
            title += " (lognormal)"
        if self.likelihood == "gauss":
            title += " (simple)"
        if self.collaboration != "ALL":
            title += f" {self.collaboration} only"
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
        nSRs = len(Ptot)
        plotStats = True
        if "plotStats" in self.options:
            plotStats = self.options["plotStats"]
        if plotStats:
            plt.text ( .67, -.12, f"this plot contains {nSRs} SRs from {nAnas} analyses", transform=ax.transAxes, c="black", fontsize=7 )
        # plt.ylabel ( "# Signal Regions" )
        self.pprint ( f"plotting {self.outfile}", verbose = 5 )
        if self.comment != None:
            plt.text ( .65, -.11, self.comment, transform=ax.transAxes,
                       style="italic" )
        if self.disclaimer:
            plt.text ( .3, .3, "do not circulate!", transform=ax.transAxes,
                       rotation=35, c="#ff3333", fontsize=20 )
        plt.kittyPlot ( self.outfile, self.show )

        plt.clf()
        plt.close()

def getArgs( cmdline = None ):
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
    argparser.add_argument ( '-w', '--weighted',
            help='weighted plot, i.e. each analysis (not each SR) counts equally', action='store_true' )
    argparser.add_argument ( '-F', '--fakes',
            help='add the fakes to the plot', action='store_true' )
    argparser.add_argument ( '-p', '--pvalues',
            help='plot p-values, not significances', action='store_true' )
    argparser.add_argument ( '-b', '--before',
            help='plot only entries before a certain date, like 2017/2/27', 
            type=str, default=None )
    argparser.add_argument ( '-S', '--signalmodel',
            help='use the signal+bg model for computing likelihoods', action='store_true' )
    argparser.add_argument ( '-l', '--likelihood', nargs='?',
            help='likelihood: gauss (g), gauss+poisson (gp), or lognormal+poisson (lp) [gauss+poisson]',
            type=str, default="gauss+poisson" )
    argparser.add_argument ( '-t', '--topologies', nargs='?',
            help='filter for certain topologies, e.g. T1, T2tt. Comma separated. The signal region must have a map for any one of the given topologies. "^" before the name acts as negation [None]',
            type=str, default=None )
    argparser.add_argument ( '--sqrts', nargs='*',
            help='sqrtses [8,13,13.6]', type=float, default=[8,13,13.6] )
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
    argparser.add_argument ( '-T', '--title', nargs='?',
            help='supply an alternative title [None]',
            type=str, default=None )
    argparser.add_argument ( '-D', '--disclaimer',
            help='add a disclaimer', action='store_true' )
    argparser.add_argument ( '-O', '--options',
            help='dictionary of options, given as string {try xlabel, ylabel, plotStats, plot_averages, weighted, yrange} [None]',
            type=str, default=None )
    argparser.add_argument ( '-U', '--ulalso',
            help='upper limit results also (but also if not eff maps exist for a given analysis)', action='store_true' )
    argparser.add_argument ( '-r', '--roughviz',
            help='roughviz plot', action='store_true' )
    argparser.add_argument ( '--show',
            help='show plot', action='store_true' )
    argparser.add_argument ( '--Zmax',
            help='maximum Z signifances to plot (|Z|) [3.25]',
            type=float, default=3.25 )
    if type(cmdline) in [ str ]:
        cmdline = cmdline.split()
        if "plotDBDict.py" in cmdline[0]:
            cmdline = cmdline[1:]

    args=argparser.parse_args( cmdline )
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
    :param cmdline: the command line arguments, e.g "-d ./db222pre1.dict  -r"
    :param options: additional options
    :returns: plotter object
    """
    args = getArgs( cmdline )
    plotter = Plotter ( args )

    if args.roughviz:
        ret, _ = plotter.rough( args.outfile, options )
    else:
        ret = plotter.plot( args.outfile )
    return plotter

if __name__ == "__main__":
    main()
