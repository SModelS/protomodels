#!/usr/bin/env python3

""" the plotting script for the llhd scans """

from smodels.tools.physicsUnits import TeV
import pickle, sys, copy, subprocess, os, colorama, time, glob, math
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
    # print ( f"{colorama.Fore.RED}[plotLlhds] plotting tool is meant to be called from protomodels top directory, not 'plotting'{colorama.Fore.RESET}" )
from ptools.sparticleNames import SParticleNames
matplotlib.rcParams['hatch.linewidth'] = .5  # previous svg hatch linewidth

def integrateLlhds ( Z, RMAX, rthreshold ):
    """ compute the integral of the likelihood over all points """
    I = 0.
    for x,row in enumerate(Z):
        for y,nll in enumerate(row):
            if RMAX[x][y]>rthreshold:
                continue
            if not np.isnan(nll):
                I += np.exp ( - nll )
    return I

def findMin ( oldZ ):
    """ find the minimum in Z """
    idx = np.nanargmin ( oldZ ) 
    y = idx % oldZ.shape[1] 
    x = int ( ( idx - y ) / oldZ.shape[1] )
    m = oldZ[x][y]
    return x,y,m

def computeHPD ( Z, RMAX, alpha = .9, verbose = True, rthreshold=1.7 ):
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
    I = integrateLlhds ( Z, RMAX, rthreshold )
    S = 0.
    points = []
    n = 0
    oldZ = copy.deepcopy ( Z )
    for x,row in enumerate(newZ):
        for y,_ in enumerate(row):
            rmax = 0.
            if type(RMAX) != type(None):
                rmax = RMAX[x][y]
            if rmax > rthreshold: ## kill the excluded areas
                oldZ[x][y]= float("nan") # oldZ[x][y] # float("nan")
            n += 1
            newZ[x][y] = 0.
    ctr = 0
    while S < alpha: ## as long as we dont have enough area
        x,y,m = findMin(oldZ)
        ctr+= 1
        S += np.exp ( -m)/I ## add up
        oldZ[x][y]=float("nan") ## kill this one
        newZ[x][y]=1 +1./ctr
    if verbose:
        print ( "%d/%d points in %d%s HPD" % ( sum(sum(newZ)), n, int(alpha*100), "%" ) )
    return newZ

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

def getPidList( pid1, rundir ):
    """ obtain the list of pids to produce plots for """
    if pid1 > 0:
        return [ pid1 ]
    pids = set()
    ## obtain pids from mp files
    # files = glob.glob ( "%s/mp*pcl" % rundir )
    files = glob.glob ( "%s/llhd*pcl" % rundir )
    for f in files:
        t = f.replace(rundir,"")
        t = t.replace("mp","")
        t = t.replace("llhd","")
        t = t.replace(".pcl","")
        t = t.replace(".pcl","")
        t = t.replace("1000022","")
        pids.add ( int(t) )
    pids = list ( pids )
    if len(pids)==0:
        print ( "[plotLlhds] could not find any llhd*pcl files. Perhaps you wish to perform ../moretools/fetchFromClip.py --llhds=" )
        sys.exit()
    print ( "[plotLlhds] creating plots for pids: %s" % ", ".join ( map(str,pids) ) )
    return pids

class LlhdPlot:
    """ A simple class to make debugging the plots easier """
    def __init__ ( self, pid1, pid2, verbose, copy, max_anas, 
                   interactive, drawtimestamp, compress, rundir,
                   upload, dbpath ):
        """
        :param pid1: pid for x axis, possibly a range of pids
        :param pid2: pid for y axis
        :param verbose: verbosity (debug, info, warn, or error)
        :param copy: copy plot to ../../smodels.github.io/protomodels/latest
        :param max_anas: maximum number of analyses on summary plot
        :param interactive: prepare for an interactive session?
        :param drawtimestamp: if true, put a timestamp on plot
        :param compress: prepare for compression
        :param upload: upload directory, default is "latest"
        :param dbpath: path to database
        """
        self.dbpath = dbpath
        self.rundir = rundir
        self.upload = upload
        self.setup( pid1, pid2 )
        self.DEBUG, self.INFO = 40, 30
        self.drawtimestamp = drawtimestamp
        self.max_anas = max_anas ## maximum number of analyses
        self.copy = copy
        self.rthreshold = 1.7
        self.interactive = interactive
        self.hiscorefile = "./hiscore.hi"
        if rundir != None:
            self.hiscorefile = f"{rundir}/hiscore.hi"
        self.setVerbosity ( verbose )
        masspoints,mx,my,nevents,topo,timestamp = self.loadPickleFile( compress )
        self.masspoints = masspoints
        self.mx = mx
        self.my = my
        self.nevents = nevents
        self.topo = topo
        self.timestamp = timestamp
        self.massdict = {}
        self.rdict = {}
        if masspoints == None:
            return
        for m in masspoints:
            self.massdict[ (m[0],m[1]) ] = m[2]
            if len(m)>3:
                self.rdict[ (m[0],m[1]) ] = m[3]

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
            self.pprint ( "I dont understand verbosity ``%s''. Setting to debug." % verbose )
            self.verbose = 40

    def getHash ( self, m1=None, m2=None ):
        """ get hash for point. if None, get hash for self.mx, self.my """
        if m1 == None:
            m1 = self.mx
        if m2 == None:
            m2 = self.my
        return int(1e3*m1) + int(1e0*m2)

    def getResultFor ( self, ana, masspoint ):
        """ return result for ana/topo pair 
        :param ana: the analysis id. optionally a data type can be specificed, e.g.
                    as :em. Alternatively, a signal region can be specified.
        :param masspoint: a point from self.masspoints
        :returns: results for this analysis (possibly data type, possibly signal region) 
                  and topology
        """
        #self.pprint ( "asking for %s" % ana )
        ret,sr = None, None
        dType = "any"
        if ":" in ana:
            ana,dType = ana.split(":")
        for k,v in masspoint.items():
            tokens = k.split(":")
            if dType == "ul" and tokens[1] != "None":
                continue
            if dType == "em" and tokens[1] == "None":
                continue
            if ana != tokens[0]:
                continue
            # self.pprint ( "asking for %s, %s %s" % ( tokens[0], tokens[1], dType ) )
            if tokens[1] != None and dType not in [ "any", "ul", "None" ]:
                # if signal regions are given, they need to match
                if tokens[1] != dType:
                    continue
                self.debug ( "found a match for", tokens[0], tokens[1], v )
            if self.topo not in tokens[2]:
                continue
            if ret == None or v > ret:
                ret = v
                sr = tokens[1]
        return ret,sr

    def loadPickleFile ( self, returnAll=False ):
        """ load dictionary from picklefile 
        :param returnAll: return all likelihoods info
        """
        topo, timestamp = "?", "?"
        allhds = None
        with open ( self.picklefile, "rb" ) as f:
            try:
                allhds = pickle.load ( f )
                mx = pickle.load ( f )
                my = pickle.load ( f )
                nevents = pickle.load ( f )
                topo = pickle.load ( f )
                timestamp = pickle.load ( f )
            except EOFError as e:
                print ( "[plotLlhds] EOF error %s, when reading %s" % \
                        ( e, self.picklefile ) )
            f.close()
        if allhds == None:
            print ( "couldnt read llhds in %s" % self.picklefile )
            return None,None,None,None,None,None
        if returnAll:
            return allhds,mx,my,nevents,topo,timestamp
        llhds=[]
        mu = 1.
        def getMu1 ( L ):
            for k,v in L.items():
                if abs(k-mu)<1e-9:
                    return v
            print ( "couldnt find anything" )
            return None
        for llhd in allhds:
            if self.pid1 in [ 1000001, 1000002, 1000003, 1000004 ]:
                if llhd[0]<310.:
                    print ( "light squark mass wall, skipping mx %d < 310 GeV" % llhd[0] )
                    continue
            if len(llhd)==4:
                llhds.append ( (llhd[0],llhd[1],getMu1(llhd[2]),llhd[3]) )
            else:
                llhds.append ( (llhd[0],llhd[1],getMu1(llhd[2]),[0.,0.,0.]) )
        return llhds,mx,my,nevents,topo,timestamp

    def pprint ( self, *args ):
        print ( "[plotLlhds] %s" % " ".join(map(str,args)) )  

    def debug ( self, *args ):
        if self.verbose >= self.DEBUG:
            print ( "[plotLlhds] %s" % " ".join(map(str,args)) )  

    def setup ( self, pid1, pid2 ):
        """ setup rundir, picklefile path and hiscore file path """
        self.hiscorefile = self.rundir + "/hiscore.hi"
        if not os.path.exists ( self.hiscorefile ):
            self.pprint ( "could not find hiscore file %s" % self.hiscorefile )
 
        self.pid1 = pid1
        self.pid2 = pid2
        if type(self.pid1) in [ tuple, list ]:
            pid1 = self.pid1[0]
        self.picklefile = "%s/llhd%d%d.pcl" % ( self.rundir, pid1, self.pid2 )
        if not os.path.exists ( self.picklefile ):
            llhdp = self.picklefile
            self.picklefile = "%s/mp%d%d.pcl" % ( self.rundir, pid1, self.pid2 )
        if not os.path.exists ( self.picklefile ):
            self.pprint ( "could not find pickle files %s and %s" % \
                          ( llhdp, self.picklefile ) )

    def describe ( self ):
        """ describe the situation """
        print ( "%d masspoints obtained from %s, hiscore stored in %s" % \
                ( len ( self.masspoints), self.picklefile, self.hiscorefile ) )
        print ( "Data members: plot.masspoints, plot.massdict, plot.timestamp, plot.mx, plot.my" )
        print ( "              plot.pid1, plot.pid2, plot.topo" )
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
        if False: ## set to true and we have the old analysis Ids
            return anaid
        if not hasattr ( self, "database" ):
            from smodels.experiment.databaseObj import Database
            #dbname = "./original.pcl" 
            #dbname = "/home/walten/git/smodels-database"
            #dbname = "/scratch-cbe/users/wolfgan.waltenberger/rundir/db31.pcl"
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

    def plot ( self, ulSeparately=True, pid1=None, dbpath = "official" ):
        """ a summary plot, overlaying all contributing analyses 
        :param ulSeparately: if true, then plot UL results on their own
        """
        if pid1 == None and type(self.pid1) in [ list, tuple ]:
            for p in self.pid1:
                self.plot ( ulSeparately, p )
            return
        if type(pid1) in [ tuple, list ]:
            for p in pid1:
                self.plot ( ulSeparately, p )
            return
        if pid1 == None:
            pid1 = self.pid1
        self.pprint ( "plotting summary for %s, %s" % ( pid1, self.topo ) )
        resultsForPIDs = {}
        from plotting.plotHiscore  import HiscorePlotter
        plotter= HiscorePlotter()
        protomodel = plotter.obtain ( 0, self.hiscorefile, dbpath = dbpath )
        for tpred in protomodel.bestCombo:
            resultsForPIDs = plotter.getPIDsOfTPred ( tpred, resultsForPIDs, integrateSRs=False )
        stats = self.getAnaStats( integrateSRs=False )
        if stats == None:
            self.pprint ( "found no ana stats?" )
            return
        anas = list(stats.keys())
        if pid1 in resultsForPIDs:
            self.debug ( "results for PIDs %s" % ", ".join ( resultsForPIDs[pid1] ) )
            anas = list ( resultsForPIDs[pid1] )
        anas.sort()
        self.pprint ( "summary plot: %s" % ", ".join ( anas ) )
        # print ( stats.keys() )
        colors = [ "red", "green", "blue", "orange", "cyan", "magenta", "grey", "brown",
                   "pink", "indigo", "olive", "orchid", "darkseagreen", "teal" ]
        xmin,xmax,ymin,ymax=9000,0,9000,0
        for m in self.masspoints:
            if m[0] < xmin:
                xmin = m[0]
            if m[0] > xmax:
                xmax = m[0]
            if m[1] < ymin:
                ymin = m[1]
            if m[1] > ymax:
                ymax = m[1]
        if abs(xmin-310.)<1e-5:
            xmin=330. ## cut off the left margin
        print ( "[plotLlhds] range x [%d,%d] y [%d,%d]" % ( xmin, xmax, ymin, ymax ) )
        handles = []
        existingPoints = []
        combL = {}
        namer = SParticleNames ( susy = False )
        for ctr,ana in enumerate ( anas ): ## loop over the analyses
            if ctr >= self.max_anas:
                self.pprint ( "too many (%d > %d) analyses." % (len(anas),self.max_anas) )
                for ana in anas[ctr:]:
                    self.pprint ( "  - skipping %s" % ana )
                break
            color = colors[ctr]
            x,y=set(),set()
            L, R = {}, {}
            minXY=( 0.,0., float("inf") )
            s=""
            r,sr = self.getResultFor ( ana, self.masspoints[0][2] )
            if r:
                s="(%.2f)" % (-np.log(r))
            print ( "[plotLlhds] result for", ana,"is", s )
            cresults = 0
            for cm,masspoint in enumerate(self.masspoints[1:]):
                #if cm % 10 != 0:
                #    continue
                if cm % 1000 == 0:
                    print ( ".", end="", flush=True )
                m1,m2,llhds,robs=masspoint[0],masspoint[1],masspoint[2],masspoint[3]
                rmax=float("nan")
                if len(robs)>0:
                    rmax=robs[0]
                if m2 > m1:
                    print ( "m2,m1 mass inversion?",m1,m2 )
                x.add ( m1 )
                y.add ( m2 )
                zt = float("nan")
                result,sr = self.getResultFor ( ana, llhds )
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
            # print ( "\n[plotLlhds] min(xy) for %s is at m=(%d/%d): %.2f(%.2g)" % ( ana, minXY[0], minXY[1], minXY[2], np.exp(-minXY[2] ) ) )
            if cresults == 0:
                print ( f"[plotLlhds] warning: found no results for {masspoint}. skip" )
                continue
                # return
            x.add ( xmax*1.03 )
            x.add ( xmin*.93 )
            y.add ( ymax+50. )
            y.add ( 0. )
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
            hldZ100 = computeHPD ( Z, None, 1., False, rthreshold=self.rthreshold )
            cont100 = plt.contour ( X, Y, hldZ100, levels=[0.25], colors = [ color ], linestyles = [ "dotted" ], zorder=10 )
            #hldZ95 = computeHPD ( Z, .95, False )
            #cont95 = plt.contour ( X, Y, hldZ95, levels=[0.5], colors = [ color ], linestyles = [ "dashed" ] )
            #plt.clabel ( cont95, fmt="95%.0s" )
            hldZ50 = computeHPD ( Z, RMAX, .68, False, rthreshold=self.rthreshold )
            cont50c = plt.contour ( X, Y, hldZ50, levels=[1.0], colors = [ color ], zorder=10 )
            cont50 = plt.contourf ( X, Y, hldZ50, levels=[1.,10.], colors = [ color, color ], alpha=getAlpha( color ), zorder=10 )
            plt.clabel ( cont50c, fmt="68%.0s" )
            if hasattr ( cont50, "axes" ):
                ax = cont50.axes
            else:
                ax = cont50.ax
            while isCloseToExisting ( minXY, existingPoints ):
                minXY = ( minXY[0]+8., minXY[1]+8., minXY[2] )
            a = ax.scatter( [ minXY[0] ], [ minXY[1] ], marker="*", s=180, color="black", zorder=20 )
            anan = ana.replace(":None",":UL") # + " (%.2f)" % (minXY[2])
            label = self.getPrettyName ( ana )
            a = ax.scatter( [ minXY[0] ], [ minXY[1] ], marker="*", s=110, color=color, 
                            label=label, alpha=1., zorder=20 )
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
        contRMAX = plt.contour ( X, Y, RMAX, levels=[self.rthreshold], colors = [ "gray" ], zorder=10 )
        contRMAXf = plt.contourf ( X, Y, RMAX, levels=[self.rthreshold,float("inf")], colors = [ "gray" ], hatches = ['////'], alpha=getAlpha( "gray" ), zorder=10 )
        hldZcomb68 = computeHPD ( ZCOMB, RMAX, .68, False, rthreshold=self.rthreshold )
        contZCOMB = plt.contour ( X, Y, hldZcomb68, levels=[.25], colors = [ "black" ], zorder=10 )

        # ax.scatter( [ minXY[0] ], [ minXY[1] ], marker="s", s=110, color="gray", label="excluded", alpha=.3, zorder=20 )
        print()
        self.pprint ( "timestamp:", self.timestamp, self.topo, max(x) )
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
        # Xs,Ys=X,Y
        Xs,Ys = filterSmaller ( X, Y )
        h = self.getHash()
        # print ( "hash is", h )
        #s=" (??)"
        #if h in L:
        #    s=" (%.2f)" % L[h]
        #s=" (%.2f)" % self.getLClosestTo ( L )
        s=""
        ax.scatter( [ self.mx ], [ self.my ], marker="*", s=200, color="white", zorder=20 )
        c = ax.scatter( [ self.mx ], [ self.my ], marker="*", s=160, color="black", 
                      label="proto-model%s" % s, zorder=20 )
        handles.append ( c )
        if sr == None:
            sr = "UL"
        # plt.title ( "HPD regions, %s [%s]" % ( namer.texName(pid1, addSign=False, addDollars=True), self.topo ), fontsize=14 )
        plt.xlabel ( "m(%s) [GeV]" % namer.texName(pid1,addSign=False, addDollars=True), fontsize=14 )
        plt.ylabel ( "m(%s) [GeV]" % namer.texName(self.pid2, addSign=False, addDollars=True), fontsize=14 )
        circ1 = mpatches.Patch( facecolor="gray",alpha=getAlpha("gray"),hatch=r'////',label='excluded by critic (r>1.7)', edgecolor="black" )
        handles.append ( circ1 )
        plt.legend( handles=handles, loc="upper left", fontsize=12 )
        figname = "%s/llhd%d.png" % ( self.rundir, pid1 )
        self.pprint ( "saving to %s" % figname )
        plt.savefig ( figname )
        if self.interactive:
            self.axes = ax
            self.plt = plt
        plt.close()
        if self.copy:
            self.copyFile ( figname )
        return

    def copyFile ( self, filename ):
        """ copy filename to smodels.github.io/protomodels/<upload>/ """
        dest = os.path.expanduser ( "~/git/smodels.github.io" )
        cmd = "cp %s %s/protomodels/%s/" % ( filename, dest, self.upload )
        o = subprocess.getoutput ( cmd )
        self.pprint ( "%s: %s" % ( cmd, o ) )


    def getAnaStats ( self, integrateSRs=True, integrateTopos=True,
                      integrateDataType=True  ):
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
            m1,m2,llhds=masspoint[0],masspoint[1],masspoint[2]
            if len(masspoint)>3:
                robs = masspoint[3]
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
        print ( "%6d masspoints with %s" % ( len(self.masspoints), self.topo ) )
        for k,v in stats.items():
            print ( "%6d: %s" % ( v, k ) )

    def compress ( self ):
        """ produce a pcl file with only a fraction of the points. 
            good for testing and development """
        backupfile = self.picklefile.replace(".pcl",".bu.pcl")
        subprocess.getoutput ( "cp %s %s" % ( self.picklefile, backupfile ))
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

    def findClosestPoint ( self, m1=None, m2=None, nll=False ):
        """ find the mass point closest to m1, m2. If not specified, 
            return the hiscore point.
        :param nll: if True, report nlls, else report likelihoods.
        """
        if m1 == None:
            m1 = self.mx
        if m2 == None:
            m2 = self.my
        dm,point = float("inf"),None
        def distance ( m ):
            return (m[0]-m1)**2 + (m[1]-m2)**2

        for m in self.masspoints:
            tmp = distance(m)
            if tmp < dm:
                dm = tmp
                point = m
        if not nll:
            return point
        # asked for NLLs
        D = {}
        for k,v in point[2].items():
            D[k]=-np.log(v)
        return ( point[0], point[1], D )

    def interact ( self ):
        import IPython
        varis = "plot.describe()"
        print ( "%s[plot] interactive session. Try: %s%s" % \
                ( colorama.Fore.GREEN, varis, colorama.Fore.RESET ) )
        IPython.embed( using=False )

if __name__ == "__main__":
    import argparse
    argparser = argparse.ArgumentParser(
            description='plot likelihoods scans')
    argparser.add_argument ( '-v', '--verbose',
            help='verbosity: debug, info, warn, or error [warn]',
            type=str, default="warn" )
    argparser.add_argument ( '-1', '--pid1',
            help='pid1, if 0 then search for llhd*pcl files [0]',
            type=int, default=0 )
    argparser.add_argument ( '-M', '--max_anas',
            help='maximum number of analyses [4]',
            type=int, default=4 )
    argparser.add_argument ( '-N', '--notimestamp',
            help='dont put a timestamp on it',
            action="store_true" )
    argparser.add_argument ( '-2', '--pid2',
            help='pid2 [1000022]',
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
    args = argparser.parse_args()
    drawtimestamp = not args.notimestamp

    rundir = gsetup( args.rundir )
    pids = getPidList ( args.pid1, rundir )

    if args.interactive and len(pids)>1:
        print ( "[plotLlhds] interactive mode plus several plots. interactive is only for one plot." )
        args.interactive = False

    for pid1 in pids:
        plot = LlhdPlot ( pid1, args.pid2, args.verbose, args.copy, args.max_anas,
                          args.interactive, drawtimestamp, args.compress, rundir,
                          args.upload )

        if args.list_analyses:
            plot.listAnalyses()

        if args.compress:
            plot.compress()
            sys.exit()

        plot.plot()

        if args.interactive:
            plot.interact()

