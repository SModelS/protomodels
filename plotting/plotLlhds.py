#!/usr/bin/env python3

""" the plotting script for the llhd scans """

from smodels.base.physicsUnits import TeV
import pickle, sys, copy, subprocess, os, time, glob, math
from colorama import Fore as ansi
import IPython
import numpy as np
from protomodels.csetup import setup
setup()
import matplotlib
matplotlib.use("Agg")
from matplotlib import pyplot as plt
import matplotlib.patches as mpatches
if "plotting" in os.getcwd():
    sys.path.insert(0,"../")
from ptools.sparticleNames import SParticleNames
matplotlib.rcParams['hatch.linewidth'] = .5  # previous svg hatch linewidth
from protomodels.builder.loggerbase import LoggerBase
from typing import Dict, Tuple, Union

def findMin ( oldZ ):
    """ find the minimum in Z """
    idx = np.nanargmin ( oldZ ) 
    y = idx % oldZ.shape[1] 
    x = int ( ( idx - y ) / oldZ.shape[1] )
    m = oldZ[x][y]
    return x,y,m

def filterSmaller ( X, Y ):
    """ filter out whenever X < Y """
    Xs,Ys = [], []
    for irow,row in enumerate ( zip ( X, Y ) ):
        xt = []
        yt = []
        if irow % 3 == 1:
            continue
        if irow % 3 == 2:
            continue
        for icol,col in enumerate ( zip ( row[0], row[1] ) ):
            if col[0]>col[1]: ## all is good
                xt.append ( float(col[0]) )
                yt.append ( float(col[1]) )
            else:
                xt.append ( float("nan") )
                yt.append ( float("nan") )
        Xs.append ( xt )
        Ys.append ( yt )
    return np.array(Xs), np.array(Ys)
            
def isCloseToExisting ( minXY, existingPoints ):
    """ is the point at minXY close to any existing Points? """
    for ep in existingPoints:
        d2 = (minXY[0] - ep[0])**2 + (minXY[1] - ep[1])**2
        if d2 < 5.:
            return True
    return False

def getAlpha ( color ):
    """ different alpha for different colors """
    rets = { "red": .3, "gray": .4 }
    if color in rets:
        return rets[color]
    return .3

def getPidList( xvariable, rundir ):
    """ obtain the list of pids to produce plots for """
    if xvariable > 0:
        return [ xvariable ]
    pids = set()
    ## obtain pids from mp files
    # files = glob.glob ( "%s/mp*pcl" % rundir )
    files = glob.glob ( f"{rundir}/llhd*pcl" )
    for f in files:
        t = f.replace(rundir,"")
        t = t.replace("mp","")
        t = t.replace("llhd","")
        t = t.replace(".pcl","")
        t = t.replace(".pcl","")
        t = t.replace("1000022","")
        print ( f"[plotLlhds] @@2 adding {t} {type(t)}" )
        pids.add ( self.namer.pid(t) )
        # pids.add ( int(t) )
    pids = list ( pids )
    if len(pids)==0:
        print ( "[plotLlhds] could not find any llhd*pcl files. Perhaps you wish to perform ../moretools/fetchFromClip.py --llhds=" )
        sys.exit()
    print ( f"[plotLlhds] creating plots for pids: {', '.join ( map(str,pids) )}" )
    return pids

class LlhdPlot ( LoggerBase ):
    """ A simple class to make debugging the plots easier """
    def __init__ ( self, xvariable, yvariable, verbose, copy, max_anas, 
                   interactive, drawtimestamp, compress, rundir,
                   upload, dbpath ):
        """
        :param xvariable: pid for x axis, possibly a range of pids
        :param yvariable: pid for y axis
        :param verbose: verbosity (debug, info, warn, or error)
        :param copy: copy plot to ../../smodels.github.io/protomodels/latest
        :param max_anas: maximum number of analyses on summary plot
        :param interactive: prepare for an interactive session?
        :param drawtimestamp: if true, put a timestamp on plot
        :param compress: prepare for compression
        :param upload: upload directory, default is "latest"
        :param dbpath: path to database
        """
        super ( LlhdPlot, self ).__init__ ( 0 )
        self.namer = SParticleNames ( susy = False )
        xvariable, yvariable = self.namer.pid ( xvariable ), self.namer.pid ( yvariable )
        self.dbpath = dbpath
        self.usePrettyNames = False
        self.rundir = rundir
        self.upload = upload
        self.setup( xvariable, yvariable )
        self.DEBUG, self.INFO = 40, 30
        self.drawtimestamp = drawtimestamp
        self.max_anas = max_anas ## maximum number of analyses
        self.copy = copy
        self.rthreshold = 1.7
        self.interactive = interactive
        self.hiscorefile = "./hiscores.dict"
        if rundir != None:
            self.hiscorefile = f"{rundir}/hiscores.dict"
        self.setVerbosity ( verbose )
        self.compress = compress
        masspoints,mx,my,nevents,topo,timestamp = self.loadPickleFile( compress )
        self.masspoints = masspoints
        self.mx = mx
        self.my = my
        self.nevents = nevents
        self.topo = topo
        from ptools.moreHelpers import namesForSetsOfTopologies
        self.toponames = namesForSetsOfTopologies ( self.topo )[0]
        self.timestamp = timestamp
        self.massdict = {}
        self.rdict = {}
        if masspoints == None:
            return
        countCritics = {} ## count occurrences of analyses in critic
        # to determine which analyses dominate the critic
        for m in masspoints:
            masstuple = (m["mx"],m["my"])
            self.massdict[ masstuple ] = m["llhd"]
            self.rdict[ masstuple ] = m["critic"]
            for ana,r in m["critic"].items():
                if r>self.rthreshold:
                    if not ana in countCritics:
                        countCritics[ana]=0
                    countCritics[ana]+=1
        significances = {} ## get the Z's at the winner
        for name,eul in masspoints[0]["eul"].items():
            if not "oul" in masspoints[0]:
                continue
            if not name in masspoints[0]["oul"]:
                continue
            oul = masspoints[0]["oul"][name]
            if eul is None or oul is None:
                continue
            sigma_exp = eul / 1.96
            Z = ( oul - eul ) / sigma_exp
            if Z in significances:
                Z+=1e-6
            significances[name] = Z
        self.significances = significances

        self.criticsOccurences = countCritics
        self.printCritic()

    def printCritic ( self ):
        self.pprint ( "Leading critics:" )
        reverted = { v:k for k,v in self.criticsOccurences.items() }
        keys = list ( reverted.keys())
        keys.sort( reverse=True )
        for k in keys[:3]:
            v = reverted[k]
            self.pprint ( f"    {v}: {k} vetoes" )

    def getMostOutspokenCritic ( self ) -> str:
        """ given self.criticsOccurences, report the most
        important critic """
        ret,count="?",0
        for k,v in self.criticsOccurences.items():
            if v > count:
                ret=k
                count=v
        ret = ret.replace("(comb)","")
        return ret

    def integrateLlhds ( self, Z, RMAX ):
        """ compute the integral of the likelihood over all points """
        I = 0.
        for x,row in enumerate(Z):
            for y,nll in enumerate(row):
                if RMAX[x][y]>self.rthreshold:
                    continue
                if not np.isnan(nll):
                    I += np.exp ( - nll )
        return I

    def computeHPD ( self, Z, RMAX, alpha : float = .9, verbose : bool = True ):
        """ compute the regions of highest posterior density to the alpha quantile
        """
        newZ = copy.deepcopy ( Z )
        if alpha > .999999: # give all points with likelihoods
            for x,row in enumerate(newZ):
                for y,_ in enumerate(row):
                    if _ > 0.:
                        newZ[x][y]=1.
                    else:
                        newZ[x][y]=0.
            return newZ
        I = self.integrateLlhds ( Z, RMAX )
        S = 0.
        points = []
        n = 0
        oldZ = copy.deepcopy ( Z )
        for x,row in enumerate(newZ):
            for y,_ in enumerate(row):
                rmax = 0.
                if type(RMAX) != type(None):
                    rmax = RMAX[x][y]
                if rmax > self.rthreshold: ## kill the excluded areas
                    oldZ[x][y]= float("nan") # oldZ[x][y] # float("nan")
                n += 1
                newZ[x][y] = 0.
        ctr = 0
        if not np.all(np.isnan(oldZ)):
            while S < alpha: ## as long as we dont have enough area
                x,y,m = findMin(oldZ)
                ctr+= 1
                S += np.exp ( -m)/I ## add up
                oldZ[x][y]=float("nan") ## kill this one
                newZ[x][y]=1 +1./ctr
        if verbose:
            print ( f"{sum(sum(newZ))}/{n} points in {int(alpha*100)}% HPD" )
        return newZ

    def topoMatches ( self, topo ):
        """ does topo match self.topo """
        ret = topo in self.toponames
        return ret

    def setVerbosity ( self, verbose ):
        self.verbose = verbose
        if type(verbose)==str:
            verbose = verbose.lower()
            if "deb" in verbose:
                self.verbose = 40
                return
            if "inf" in verbose:
                self.verbose = 30
                return
            if "warn" in verbose:
                self.verbose = 20
                return
            if "err" in verbose:
                self.verbose = 10
                return
            self.pprint ( f"I dont understand verbosity ``{verbose}''. Setting to debug." )
            self.verbose = 40

    def getHash ( self, m1=None, m2=None ):
        """ get hash for point. if None, get hash for self.mx, self.my """
        if m1 == None:
            m1 = self.mx
        if m2 == None:
            m2 = self.my
        return int(1e3*m1) + int(1e0*m2)

    def getHighestLlhdFor ( self, ana : str, llhddict : Dict ) -> Dict:
        """ return llhds for ana

        :param ana: the analysis id. optionally a data type can be specificed, e.g.
        as :em. Alternatively, a signal region can be specified.
        :param llhddict: a dictionary of likelihoods of one of the llhddicts, 
        e.g. self.llhddicts[0]["llhd"]

        :returns: Dictionary with highest likelihood and name of signal region, 
        (both being None if nothing is found)
        """
        max_llhd,sr = None, None
        dType = "any"
        n_ana = ana
        if ":" in ana:
            n_ana,dType = ana.split(":")
            if dType == "(combined)":
                dType = "(comb)"
        mname = "?"
        for name,llhd in llhddict.items():
            tokens = name.split(":")
            if dType == "ul" and tokens[1] != "None":
                continue
            if dType == "em" and tokens[1] == "None":
                continue
            if n_ana != tokens[0]:
                continue
            if tokens[1] != None and dType not in [ "any", "ul", "None" ]:
                # if signal regions are given, they need to match
                if tokens[1] != dType:
                    continue
                self.debug ( f"found a match for {tokens[0]}, {tokens[1]}, l={llhd}" )
            if not self.topoMatches ( tokens[2] ):
                self.pprint ( f"topology {tokens[2]} does not match {self.topo}, will skip" )
                # continue
            if max_llhd == None or llhd > max_llhd:
                max_llhd = llhd
                sr = tokens[1]
                mname = name
        ret = { "llhd": max_llhd, "sr": sr }
        ret["Z"] = -30.
        if mname in self.significances:
            ret["Z"] = self.significances[mname]
        return ret

    def writeScriptFile ( self ):
            syv = "_"+self.namer.asciiName(self.yvariable)
            syv = syv.replace(",","").replace(" ","").replace("~","m")
            if syv == "_X1Z":
                syv = ""
            scriptfilename = f"llhdPlot_{self.namer.asciiName(self.xvariable)}{syv}.py"
            with open ( scriptfilename, "wt" ) as f:
                print ( f"[llhdScanner] created llhdPlotScript.py" )
                f.write ( "#!/usr/bin/env python3\n\n" )
                f.write ( "import sys\n" )
                f.write ( "interactive=False\n" )
                f.write ( "if '-i' in sys.argv:\n" )
                f.write ( "    interactive=True\n" )
                f.write ( "from plotting import plotLlhds\n" )
                f.write ( f"plot = plotLlhds.LlhdPlot ( xvariable={self.xvariable}, yvariable={self.yvariable}, verbose='{self.verbose}', copy={self.copy},\n" )
                f.write ( f"    max_anas={self.max_anas}, interactive=interactive, drawtimestamp={self.drawtimestamp}, compress={self.compress},\n" )
                f.write ( f"    rundir='{self.rundir}',\n" )
                f.write ( f"    upload='{self.upload}', dbpath='{self.dbpath}' )\n" )
                f.write ( f"plot.plot()\n" )
                f.write ( f"if '-s' in sys.argv:\n" )
                f.write ( f"    plot.show()\n" )
                f.close()
                os.chmod ( scriptfilename, 0o755 )

    def loadPickleFile ( self, returnAll=False ):
        """ load dictionary from picklefile 
        :param returnAll: return all likelihoods info
        """
        topo, timestamp = "?", "?"
        masspoints = None
        ctr = 0
        success = False
        while ctr < 15 and not success:
            with open ( self.picklefile, "rb" ) as f:
                try:
                    dic = pickle.load ( f )
                    masspoints = dic["masspoints"]
                    mx = dic["mxvariable"]
                    my = dic["myvariable"]
                    nevents = dic["nevents"]
                    topo = dic["topo"]
                    timestamp = dic["timestamp"]
                    success = True
                except Exception as e:
                    self.pprint ( f"Exception {e}, when reading {self.picklefile}")
                    ctr += 1
                    time.sleep (.1 * ctr )
                f.close()
        self.pprint ( f"loaded {len(masspoints)} masspoints." )
        if masspoints == None:
            self.pprint ( f"couldnt read llhds in {self.picklefile}" )
            return None,None,None,None,None,None
        if returnAll:
            return masspoints,mx,my,nevents,topo,timestamp
        llhds=[]
        mu = 1.
        def getMu1 ( L ):
            for k,v in L.items():
                if abs(k-mu)<1e-9:
                    return v
            print ( "couldnt find anything" )
            return None
        for point in masspoints:
            if not "llhd" in point:
                print ( f"point has no llhds? {point.keys()}" )
            if self.xvariable in [ 1000001, 1000002, 1000003, 1000004 ]:
                if point['mx']<310.:
                    print ( f"light squark mass wall, skipping mx {llhd['mx']} < 310 GeV" )
                    continue
            app = copy.deepcopy ( point )
            if type(app)==tuple:
                print ( f"FIXME why is this a tuple? {app[:2]}" )
                app= { "mx": point[0], "my": point[1], "llhd": getMu1(point[2]),
                       "critic": point[3] }
            else:
                app["llhd"] = getMu1(point["llhd"])
            llhds.append ( app )
        return llhds,mx,my,nevents,topo,timestamp

    def setup ( self, xvariable, yvariable ):
        """ setup rundir, picklefile path and hiscore file path """
        self.hiscorefile = self.rundir + "/hiscores.dict"
        if not os.path.exists ( self.hiscorefile ):
            self.pprint ( f"could not find hiscore file {self.hiscorefile}" )
 
        self.xvariable = xvariable
        self.yvariable = yvariable
        if type(self.xvariable) in [ tuple, list ]:
            xvariable = self.xvariable[0]
        self.picklefile = f"{self.rundir}/llhd{self.namer.asciiName(xvariable)}{self.namer.asciiName(self.yvariable).replace(',','').replace(' ','')}.pcl"
        if not os.path.exists ( self.picklefile ):
            llhdp = self.picklefile
            self.picklefile = f"{self.rundir}/mp{self.namer.asciiName(xvariable)}{self.namer.asciiName(self.yvariable)}.pcl" 
        if not os.path.exists ( self.picklefile ):
            self.pprint(f"could not find pickle files {llhdp} and {self.picklefile}")

    def describe ( self ):
        """ describe the situation """
        print ( f"{len ( self.masspoints)} masspoints obtained from {self.picklefile}, hiscore stored in {self.hiscorefile}" )
        print ( "Data members: plot.masspoints, plot.massdict, plot.timestamp, plot.mx, plot.my" )
        print ( "              plot.xvariable, plot.yvariable, plot.topo" )
        print ( "Function members: plot.findClosestPoint()" )


    def getLClosestTo ( self, L, mx=None, my=None ):
        """ get the L closest to your point """
        if mx == None:
            mx=self.mx
        if my == None:
            my=self.my
        def distance_ ( k, mx, my ):
            _x = int(math.floor(k/1000.))
            _y = int(math.floor(k % 1000 ) )
            ret= (mx - _x)**2 + (my - _y)**2
            return ret

        dmmin, vmin = float("inf"), 23.
        for k,v in L.items():
            dm = distance_ ( k, mx, my )
            if dm < dmmin and not np.isnan(v):
                dmmin = dmmin
                vmin = v
        return vmin

    def getPrettyName ( self, anaid ):
        """ get pretty name of ana id """
        if not self.usePrettyNames: ## set to true and we have the old analysis Ids
            anaid = anaid.replace("(combined)","(comb)" )
            return anaid
        if not hasattr ( self, "database" ):
            from smodels.experiment.databaseObj import Database
            #dbname = "./original.pcl" 
            # dbname = "official"
            self.database = Database ( self.dbpath )
        from smodels_utils.helper.prettyDescriptions import prettyTexAnalysisName
        if ":" in anaid:
            anaid = anaid[:anaid.find(":")]
        ers = self.database.getExpResults ( analysisIDs = [ anaid ] )
        for er in ers:
           if hasattr ( er.globalInfo, "prettyName" ):
              pn = er.globalInfo.prettyName
              sqrts = er.globalInfo.sqrts.asNumber(TeV)
              ret = prettyTexAnalysisName ( pn, sqrts, dropEtmiss = True,
                                        collaboration = True, anaid = er.globalInfo.id )
              # for the 2020 paper to be consistent
              ret = ret.replace( "+ top tag", "stop" )
              ret = ret.replace( "+ 4 (1 b-)jets", "multijet" )
              # ret += " -> " + anaid
              return ret
        # print ( "found no pretty name", ers[0].globalInfo )
        return anaid

    def plot ( self, ulSeparately : bool = True, xvariable : Union[None,int] = None, 
               dbpath : str = "official" ):
        """ a summary plot, overlaying all contributing analyses 

        :param ulSeparately: if true, then plot UL results on their own
        """
        if xvariable == None and type(self.xvariable) in [ list, tuple ]:
            for p in self.xvariable:
                self.plot ( ulSeparately, p )
            return
        if type(xvariable) in [ tuple, list ]:
            for p in xvariable:
                self.plot ( ulSeparately, p )
            return
        if xvariable == None:
            xvariable = self.xvariable
        self.pprint ( f"plotting likelihoods for {self.namer.asciiName(xvariable)}: {self.topo}" )
        resultsForPIDs = {}
        ## this is just to obtain the hiscore
        from plotting.plotHiscore  import HiscorePlotter
        plotter= HiscorePlotter()
        protomodel = plotter.obtain ( 0, self.hiscorefile, dbpath = dbpath )
        for tpred in protomodel.bestCombo:
            resultsForPIDs = plotter.getPIDsOfTPred ( tpred, resultsForPIDs, 
                                integrateSRs=False )
        stats = self.getAnaStats( integrateSRs=False )
        if stats == None:
            self.pprint ( "found no ana stats?" )
            return
        anas = list(stats.keys())
        if xvariable in resultsForPIDs:
            self.debug ( f"results for PIDs {', '.join(resultsForPIDs[xvariable])}" )
            anas = list ( resultsForPIDs[xvariable] )
        anas.sort()
        self.pprint ( f"summary plot: {', '.join ( anas )}" )
        colors = [ "red", "green", "blue", "orange", "cyan", "magenta", "grey", 
            "brown", "pink", "indigo", "olive", "orchid", "darkseagreen", "teal" ]
        xmin,xmax,ymin,ymax=9000,0,9000,0
        for m in self.masspoints:
            if m["mx"] < xmin:
                xmin = m["mx"]
            if m["mx"] > xmax:
                xmax = m["mx"]
            if m["my"] < ymin:
                ymin = m["my"]
            if m["my"] > ymax:
                ymax = m["my"]
        if abs(xmin-310.)<1e-5:
            xmin=330. ## cut off the left margin
        self.pprint ( f"plot ranges: x=[{xmin:.1f},{xmax:.1f}] y=[{ymin:.1f},{ymax:.1f}]" )
        handles = []
        existingPoints = []
        combL = {}
        rankthem = {}
        for ana in anas: ## loop over the analyses
            ret = self.getHighestLlhdFor ( ana, self.masspoints[0]["llhd"] )
            rankthem[ ret["Z"] ] = ana
        Zs = list ( rankthem.keys())
        Zs.sort( reverse = True )
        newanas = []
        for mZ in Zs:
            if mZ > -20.:
                newanas.append ( rankthem[mZ] )
        for ctr,ana in enumerate ( newanas ): ## loop over the analyses
            if ctr >= self.max_anas:
                self.pprint ( f"too many ({len(anas)} > {self.max_anas}) analyses." )
                for ana in anas[ctr:]:
                    self.pprint ( f"  - skipping {ana}" )
                break
            color = colors[ctr]
            x,y=set(),set()
            L, R = {}, {}
            minXY=( float("nan"),float("nan"), float("inf") )
            s="none"
            ## first, check for the hiscore point
            ret = self.getHighestLlhdFor ( ana, self.masspoints[0]["llhd"] )
            sr = ret["sr"]
            llhd = ret["llhd"]
            if llhd:
                s=f"{-np.log(llhd):.2f}"
            self.pprint ( f"{ana} @ hiscore m=({self.masspoints[0]['mx']:.1f},{self.masspoints[0]['my']:.1f}): nll_max={s}" )
            cresults = 0
            ## then, run on all other points
            for cm,masspoint in enumerate(self.masspoints[1:]):
                if cm % 100 == 0:
                    print ( ".", end="", flush=True )
                m1,m2,llhds,critic=masspoint["mx"],masspoint["my"],masspoint["llhd"],masspoint["critic"]
                rmax=float("nan")
                if len(critic)>0:
                    rmax=max(critic.values())
                if m2 > m1:
                    print ( f"m2,m1 mass inversion? {m1,m2}" )
                x.add ( m1 )
                y.add ( m2 )
                zt = float("nan")
                ret = self.getHighestLlhdFor ( ana, llhds )
                result = ret [ "llhd" ]
                sr = ret [ "sr" ]
                if result:
                    zt = - np.log( result )
                    cresults += 1
                    if zt < minXY[2] and rmax<=self.rthreshold:
                        minXY=(m1,m2,zt)
                h = self.getHash(m1,m2)
                L[h]=zt
                if not h in combL:
                    combL[h]=0.
                if np.isnan(zt):
                    combL[h] = combL[h] + 100.
                else:
                    combL[h] = combL[h] + zt
                R[h]=rmax
            print ()
            self.pprint ( f"{ana}: {cresults}/{len(self.masspoints)} results" )
            if cresults == 0:
                continue
                # return
            if False:
                x.add ( xmax*1.03 )
                x.add ( xmin*.93 )
                y.add ( ymax*1.01 )
                y.add ( ymin*.97 )
            x,y=list(x),list(y)
            x.sort(); y.sort()
            X, Y = np.meshgrid ( x, y )
            Z = float("nan")*X
            RMAX = float("nan")*X
            for irow,row in enumerate(Z):
                for icol,col in enumerate(row):
                    h = 0
                    if len(x)>= icol and len(y) >= irow:
                        h = self.getHash(list(x)[icol],list(y)[irow])
                    if h in L:
                        Z[irow,icol]=L[h]
                    if h in R:
                        RMAX[irow,icol]=R[h]
            if self.interactive:
                self.RMAX = RMAX
                # self.ZCOMB = ZCOMB
                self.Z = Z
                self.L = L
                self.R = R
                self.X = X
                self.Y = Y
            hldZ100 = self.computeHPD ( Z, None, 1., False )
            cont100 = plt.contour ( X, Y, hldZ100, levels=[0.25], colors = [ color ], linestyles = [ "dotted" ], zorder=-1 )
            #hldZ95 = self.computeHPD ( Z, .95, False )
            #cont95 = plt.contour ( X, Y, hldZ95, levels=[0.5], colors = [ color ], linestyles = [ "dashed" ] )
            #plt.clabel ( cont95, fmt="95%.0s" )
            hldZ50 = self.computeHPD ( Z, RMAX, .68, False )
            cont50c = plt.contour ( X, Y, hldZ50, levels=[1.0], colors = [ color ], zorder=-1 )
            cont50 = plt.contourf ( X, Y, hldZ50, levels=[1.,10.], colors = [ color, color ], alpha=getAlpha( color ), zorder=-1 )
            plt.clabel ( cont50c, fmt="68%.0s" )
            if hasattr ( cont50, "axes" ):
                ax = cont50.axes
            else:
                ax = cont50.ax
            while isCloseToExisting ( minXY, existingPoints ):
                if type( self.yvariable ) == tuple:
                    minXY = ( minXY[0]*1.2, minXY[1]*1.2, minXY[2] )
                else:
                    minXY = ( minXY[0]+8., minXY[1]+8., minXY[2] )
            a = ax.scatter( [ minXY[0] ], [ minXY[1] ], marker="*", s=180, color="black", zorder=20 )
            anan = ana.replace(":None",":UL") # + " (%.2f)" % (minXY[2])
            label = self.getPrettyName ( ana )
            # print ( f"@@5 ana {ana} label {label} dType {dType}" )
            a = ax.scatter( [ minXY[0] ], [ minXY[1] ], marker="*", s=110, 
                    color=color, label=label, alpha=1., zorder=20 )
            existingPoints.append ( minXY )
            handles.append ( a )
        ZCOMB = float("nan")*X
        for irow,row in enumerate(Z):
            for icol,col in enumerate(row):
                h = 0
                if len(x)> icol and len(y) > irow:
                    h = self.getHash(list(x)[icol],list(y)[irow])
                if h in combL and not np.isnan(combL[h]):
                    ZCOMB[irow,icol]=combL[h]
                    if combL[h]==0.:
                        ZCOMB[irow,icol]=float("nan")
        self.ZCOMB = ZCOMB
        contRMAX = plt.contour ( X, Y, RMAX, levels=[self.rthreshold], colors = [ "gray" ], zorder=-1 )
        contRMAXf = plt.contourf ( X, Y, RMAX, levels=[self.rthreshold,float("inf")], colors = [ "gray" ], hatches = ['////'], alpha=getAlpha( "gray" ), zorder=-1 )
        hldZcomb68 = self.computeHPD ( ZCOMB, RMAX, .68, False  )
        contZCOMB = plt.contour ( X, Y, hldZcomb68, levels=[.25], colors = [ "black" ], zorder=-1 )

        # ax.scatter( [ minXY[0] ], [ minXY[1] ], marker="s", s=110, color="gray", label="excluded", alpha=.3, zorder=20 )
        print()
        # self.pprint ( "timestamp:", self.timestamp, self.topo, max(x) )
        dx,dy = max(x)-min(x),max(y)-min(y)
        if self.drawtimestamp:
            plt.text( max(x)-.37*dx,min(y)-.11*dy,self.timestamp, c="gray" )
        ### the altitude of the alpha quantile is l(nuhat) - .5 chi^2_(1-alpha);ndf
        ### so for alpha=0.05%, ndf=1 the dl is .5 * 3.841 = 1.9207
        ### for ndf=2 the dl is ln(alpha) = .5 * 5.99146 = 2.995732
        ### folien slide 317
        if hasattr ( cont50, "axes" ):
            ax = cont50.axes
        else:
            ax = cont50.ax
        Xs,Ys = filterSmaller ( X, Y )
        h = self.getHash()
        s=""
        ax.scatter( [ self.mx ], [ self.my ], marker="*", s=200, color="white", zorder=20 )
        c = ax.scatter( [ self.mx ], [ self.my ], marker="*", s=160, color="black", 
                      label="proto-model%s" % s, zorder=20 )
        handles.append ( c )
        if sr == None:
            sr = "UL"
        if "electroweakinos_offshell" in self.topo:
            self.error ( "FIXME fix the names of the topo sets! electroweakinos!!" )
            self.topo = self.topo.replace("electroweakinos_offshell","electroweakinos" )
            self.error ( "FIXME fix the names of the topo sets! electroweakinos!!" )
            
        plt.title ( f"HPD regions, {self.namer.texName(xvariable, addSign=False, addDollars=True)} [{self.topo}]", fontsize=14 )
        plt.xlabel ( f"m({self.namer.texName(xvariable,addSign=False, addDollars=True)}) [GeV]", fontsize=14 )
        var, postfix = "m", " [GeV]"
        if type ( self.yvariable ) == tuple:
            var, postfix = "ssm", ""
        plt.ylabel ( f"{var}({self.namer.texName(self.yvariable, addSign=False, addDollars=True)}){postfix}" )
        hasCritic = np.any ( RMAX > self.rthreshold )
        if hasCritic:
            circ1 = mpatches.Patch( facecolor="gray",alpha=getAlpha("gray"),hatch=r'////',label=f'excluded by critic (r>{self.rthreshold}):\n{self.getMostOutspokenCritic()} et al', edgecolor="black" )
            handles.append ( circ1 )
        legend = ax.legend( handles=handles, loc="best", fontsize=12 )
        figname = f"{self.rundir}/llhd{self.namer.asciiName(xvariable)}_{self.namer.asciiName(self.yvariable).replace(',','').replace(' ','')}.png"
        self.pprint ( f"saving to {figname}" )
        plt.savefig ( figname )
        if self.interactive:
            self.axes = ax
            self.plt = plt
            import IPython
            IPython.embed( colors = "neutral" )
        plt.close()
        if self.copy:
            self.copyFile ( figname )
        self.figname = figname ## store filename

    def copyFile ( self, filename ):
        """ copy filename to smodels.github.io/protomodels/<upload>/ """
        dest = os.path.expanduser ( "~/git/smodels.github.io" )
        cmd = f"cp {filename} {dest}/protomodels/{self.upload}/"
        o = subprocess.getoutput ( cmd )
        self.pprint ( f"{cmd}: {o}" )


    def getAnaStats ( self, integrateSRs : bool = True, integrateTopos : bool = True,
                      integrateDataType : bool =True  ) -> Dict:
        """ given the likelihood dictionaries D, get
            stats of which analysis occurs how often 
        :param integrateTopos: sum over all topologies
        :param integrateSRs: sum over all signal regions
        :param integrateDataType: ignore data type
        """
        anas = {}
        if self.masspoints == None:
            return None
        for masspoint in self.masspoints:
            # print ( "masspoint", masspoint )
            m1,m2,llhds=masspoint["mx"],masspoint["my"],masspoint["llhd"]
            critic = masspoint["critic"]
            for k,v in llhds.items():
                tokens = k.split(":")
                if not integrateTopos and self.topo not in tokens[2]:
                    continue
                dType = ":em"
                if tokens[1] in [ "None", None ]:
                    dType = ":ul"
                name = tokens[0]
                if not integrateDataType:
                    name = name + dType
                if not integrateTopos:
                    name = tokens[0]+tokens[1]
                if not name in anas.keys():
                    anas[name]=0
                anas[name]=anas[name]+1
        return anas

    def listAnalyses( self ):
        """
        :param verbose: verbosity: debug, info, warn, or error
        """
        stats = self.getAnaStats( integrateDataType=False )
        print ( f"{len(self.masspoints):6d} masspoints with {self.topo}" )
        for k,v in stats.items():
            print ( f"{v:6d}: {k}" )

    def compress ( self ):
        """ produce a pcl file with only a fraction of the points. 
            good for testing and development """
        backupfile = self.picklefile.replace(".pcl",".bu.pcl")
        subprocess.getoutput ( f"cp {self.picklefile} {backupfile}" )
        newfile = self.picklefile.replace(".pcl",".comp.pcl")
        mx,my=set(),set()
        for m in self.masspoints:
            mx.add ( m[0] )
            my.add ( m[1] )
        mx=list(mx)
        my=list(my)

        with open ( newfile, "wb" ) as f:
            mps = []
            for i,m in enumerate(self.masspoints):
                if mx.index (m[0] ) % 2 == 0 and \
                   my.index (m[1] ) % 2 == 0:
                # if i % 5 == 0:
                    mps.append ( m )
            pickle.dump ( mps, f )
            pickle.dump ( self.mx, f )
            pickle.dump ( self.my, f )
            pickle.dump ( self.nevents, f )
            pickle.dump ( self.topo, f )
            pickle.dump ( self.timestamp, f )
            f.close()

    def pprintPointNr ( self, idx : int ):
        """ print some relevant info about mass point # idx """
        if idx >= len(self.masspoints):
            self.pprint ( f"asked to print nonexistant point {idx}, we only have {len(self.masspoints)}" )
            return
        point = self.masspoints[idx]
        print ( f"point #{idx}: m({point['mx']:.2f},{point['my']:.3f})" )
        print ( f"critic:" )
        print ( f"=======" )
        for ctr,(k,v) in enumerate(point["critic"].items()):
            print ( f"    {k:14s}: {v:.2f}" )
            if ctr > 3:
                break
        print ( f"llhds:" )
        print ( f"======" )
        for ctr,(k,v) in enumerate(point["llhd"].items()):
            print ( f"    {k:14s}: {v:.2g}" )
            if ctr > 1:
                break

    def listExcludedPoints ( self, revert : bool = False ):
        """ print to stdout the list of excluded points.
        :param revert: if true, list non-excluded points
        """
        se = "non-excluded" if revert else "excluded"
        filteredpoints = []
        for idx,point in enumerate ( self.masspoints ):
            rmax = max(point["critic"].values())
            point["idx"] = idx
            point["rmax"] = rmax
            if rmax >= self.rthreshold and not revert:
                filteredpoints.append ( point )
            if rmax <= self.rthreshold and revert:
                filteredpoints.append ( point )
        filteredpoints.sort ( key = lambda x : 1e6*x["mx"]+x["my"] )
        print ( f"List of {len(filteredpoints)} {se} points (r>{self.rthreshold}):" )
        for point in filteredpoints:
            print ( f"    #{point['idx']}: m=({point['mx']:.2f},{point['my']:.2f}) r={point['rmax']:.2f}" )


    def findClosestPoint ( self, m1 : Union[None,float]=None, 
            m2 : Union[None,float]=None, nll : bool = False ) -> Dict:
        """ find the mass point closest to m1, m2. If not specified, 
            return the hiscore point.
        :param m1: if None, use best fit point coord
        :param m2: if None, use best fit point coord
        :param nll: if True, report nlls, else report likelihoods.

        :returns: dictionary with coordinates and llhd
        """
        if m1 == None:
            m1 = self.mx
        if m2 == None:
            m2 = self.my
        dm,point, idx = float("inf"),None, 0

        def distance ( m ):
            return (m["mx"]-m1)**2 + (m["my"]-m2)**2

        for idx_,m in enumerate ( self.masspoints ):
            tmp = distance(m)
            if tmp < dm:
                dm = tmp
                point = m
                idx = idx_
        ret = { "mx": point["mx"], "my": point["my"], "dm": np.sqrt(dm), "idx": idx }
        llhdname = "llhd"
        if not nll:
            return ret
        llhdname = "nll"
        # asked for NLLs
        D = {}
        for k,v in point[2].items():
            D[k]=-np.log(v)
        ret[llhdname] = D
        return ret

    def interact ( self ):
        import IPython
        varis = "plot.describe()"
        print ( f"{ansi.GREEN}[plot] interactive session. Try: {varis}{ansi.RESET}" )
        IPython.embed( using=False )

    def show( self ):
        cmd = f"see {self.figname}"
        o = subprocess.getoutput ( cmd )
        print ( o )
        # from smodels_utils.plotting.mpkitty import timg
        # timg ( self.figname )


if __name__ == "__main__":
    import argparse
    argparser = argparse.ArgumentParser(
            description='plot likelihoods scans')
    argparser.add_argument ( '-v', '--verbose',
            help='verbosity: debug, info, warn, or error [warn]',
            type=str, default="warn" )
    argparser.add_argument ( '-1', '--xvariable',
            help='xvariable, if 0 then search for llhd*pcl files [0]',
            type=int, default=0 )
    argparser.add_argument ( '-M', '--max_anas',
            help='maximum number of analyses [4]',
            type=int, default=4 )
    argparser.add_argument ( '-N', '--notimestamp',
            help='dont put a timestamp on it',
            action="store_true" )
    argparser.add_argument ( '-2', '--yvariable',
            help='yvariable [1000022]',
            type=int, default=1000022 )
    argparser.add_argument ( '-a', '--analyses',
            help="analyses, comma separated. '*' means all analyses [*]",
            type=str, default="*" )
    argparser.add_argument ( '-l', '--list_analyses',
            help='list all analyses for these pids',
            action="store_true" )
    argparser.add_argument ( '-C', '--compress',
            help='compress the pickle file so that things work on a laptop',
            action="store_true" )
    argparser.add_argument ( '-c', '--copy',
            help='copy plots to ~/git/smodels.github.io/protomodels/<upload>',
            action="store_true" )
    argparser.add_argument ( '-u', '--upload',
            help='choose upload directory [latest]',
            type=str, default="latest" )
    argparser.add_argument ( '-R', '--rundir',
            help='override the default rundir [None]',
            type=str, default=None )
    argparser.add_argument ( '-I', '--interactive',
            help='interactive mode',
            action="store_true" )
    argparser.add_argument ( '-s', '--show',
            help='show plot at the end',
            action="store_true" )
    args = argparser.parse_args()
    drawtimestamp = not args.notimestamp

    rundir = gsetup( args.rundir )
    pids = getPidList ( args.xvariable, rundir )

    if args.interactive and len(pids)>1:
        print ( "[plotLlhds] interactive mode plus several plots. interactive is only for one plot." )
        args.interactive = False

    for xvariable in pids:
        plot = LlhdPlot ( xvariable, args.yvariable, args.verbose, args.copy, args.max_anas,
                          args.interactive, drawtimestamp, args.compress, rundir,
                          args.upload )

        if args.list_analyses:
            plot.listAnalyses()

        if args.compress:
            plot.compress()
            sys.exit()

        filename = plot.plot()

        if args.show:
            from smodels_utils.plotting.mpkitty import timg
            timg ( filename )

        if args.interactive:
            plot.interact()

