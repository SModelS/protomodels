#!/usr/bin/env python3

""" script used to produce the likelihood scans """

__all__ = [ "LlhdScanner" ]

import os, sys, multiprocessing, time, numpy, subprocess, copy, glob
from protomodels.csetup import setup
setup()
from smodels.tools.wrapperBase import WrapperBase
WrapperBase.defaulttempdir="./" ## keep the temps in our folder
from smodels.base.physicsUnits import fb
from smodels.base.runtime import nCPUs
from smodels.matching.theoryPrediction import TheoryPrediction
from tester.combiner import Combiner
from tester.predictor import Predictor
from plotting import plotLlhds
from typing import Dict, Tuple, Union, List
from ptools.sparticleNames import SParticleNames
from builder.loggerbase import LoggerBase

namer = SParticleNames ( False )

def findPids ( rundir ):
    """ search for llhd*pcl files, report the corresponding pids.
    :returns: set of pids
    """
    ret = set()
    files = glob.glob("llhd*pcl")
    files += glob.glob( f"{rundir}/llhd*pcl" )
    for f in files:
        p = f.find("llhd")
        s = f[p+4:]
        s = s.replace(".pcl","")
        s = s.replace("1000022","")
        ret.add ( int(s) )
    print ( f"[llhdScanner] pids are {ret}" )
    return ret

class LlhdThread ( LoggerBase ):
    """ one thread of the sweep """
    def __init__ ( self, threadnr: int, rundir: str,
                   protomodel, xvariable, yvariable, mxvariable, myvariable, nevents: int,
                   predictor ):
        super ( LlhdThread, self ).__init__ ( threadnr )
        self.rundir = setup( rundir )
        self.threadnr = threadnr
        self.M = copy.deepcopy ( protomodel )
        self.origmasses = copy.deepcopy ( self.M.masses )
        self.M.createNewSLHAFileName ( prefix=f"lthrd{threadnr}_{xvariable}" )
        self.xvariable = xvariable
        self.yvariable = yvariable
        self.mxvariable = mxvariable
        self.myvariable = myvariable
        self.nevents = nevents
        self.predictor = predictor

    def massesAreTied ( self, xvariable, yvariable ):
        """ are the masses of xvariable and yvariable tied originally? """
        if not xvariable in self.origmasses:
            return False
        if not yvariable in self.origmasses:
            return False
        dm = self.origmasses[xvariable] - self.origmasses[yvariable]
        if abs(dm)<1e-5:
            return True
        return False

    def getPredictions ( self, recycle_xsecs : bool = True ) -> Dict:
        """ get predictions, return likelihoods 

        :param recycle_xsecs: if true, then recycle the cross sections, dont
        recompute
        :returns: a diction with likelihoods ("llhd"), critics' responses ("critic"),
        observed ("oul") and expected ("eul") upper limits on mu.
        """
        self.M.createSLHAFile( )
        sigmacut=.02*fb
        if max(self.M.masses)>1600:
            sigmacut=.01*fb
        if max(self.M.masses)>1800:
            sigmacut=.003*fb
        if max(self.M.masses)>2000:
            sigmacut=.001*fb
        ## first get rmax
        if hasattr ( self.predictor, "predictions" ):
            del self.predictor.predictions
        worked = self.predictor.predict ( self.M, keep_predictions = True )
        
        # print ( "@@5 worked", worked, len(self.predictor.predictions) )

        ## now get the likelihoods
        llhds={}
        ## start with the SM likelihood
        llhds[0.] = self.getLikelihoods ( self.predictor.predictions, mu=0. )
        ## get for the others FIXME should adapt to ssm?
        for mu in numpy.arange(.4,1.8,.05):
            llhds[float(mu)] = self.getLikelihoods ( self.predictor.predictions, mu=mu )
        ouls = self.getLimits ( self.predictor.predictions, False )
        euls = self.getLimits ( self.predictor.predictions, True )
        del self.predictor.predictions
        self.M.delCurrentSLHA()
        critics={}
        for critic in self.M.critic_description.split(","):
            tokens = critic.split(":")
            if len(tokens)>1:
                critics[tokens[0]]=float(tokens[1])

        return { "llhd": llhds, "critic": critics, "oul": ouls, "eul": euls }

    def getLimits ( self, predictions : List[TheoryPrediction], 
                    expected : bool ) -> Dict:
        """ get the limits for all predictions 

        :param expected: if true, get expected limits on mu, else observed
        """
        limits = {}
        for tp in predictions:
            txname = ','.join ( set( [ i.txName for i in tp.txnames ] ) )
            dId = tp.dataId()
            if dId == "(combined)":
                dId = "(comb)"
            name = f"{tp.analysisId()}:{dId}:{txname}"
            limits[ name ] = tp.getUpperLimitOnMu ( expected = expected )
        return limits

    def getLikelihoods ( self, predictions, mu = 1. ) -> Dict:
        """ return dictionary with the likelihoods per analysis """
        llhds= {}
        for tp in predictions:
            txname = ','.join ( set( [ i.txName for i in tp.txnames ] ) )
            dId = tp.dataId()
            if dId == "(combined)":
                dId = "(comb)"
            name = f"{tp.analysisId()}:{dId}:{txname}"
            llhds[ name ] = tp.likelihood ( mu )
        return llhds

    def clean ( self ):
        """ clean up after the run """
        cmd = f"rm {self.M.currentSLHA}"
        subprocess.getoutput ( cmd )

    def setMass ( self, pid : int, mass : float ):
        """ set mass of <pid> to <mass> """
        partners = [ ( 1000023, 1000024 ) ]
        self.M.masses[pid]=mass
        for pair in partners:
            if not pid in pair:
                continue
            for p in pair:
                if p in self.M.masses and self.massesAreTied ( p, pid ):
                    self.M.masses[p]=mass

    def run ( self, rxvariable, ryvariable ):
        """ run for the points given """
        oldmasses = {}
        masspoints=[]
        nxvariables = len(rxvariable)
        for i1,m1 in enumerate(rxvariable):
            self.pprint ( f"now starting with {i1}/{nxvariables}" )
            self.setMass ( self.xvariable, m1 )
            self.M.masses[self.yvariable]=self.myvariable ## reset LSP mass
            for k,v in oldmasses.items():
                self.pprint ( f"WARNING: setting mass of {k} back to {v}" )
                self.M.masses[k]=v
            oldmasses={}
            self.M.delXSecs() ## make sure we compute
            xsecs = self.M.getXsecs()
            xsectot = 0.*fb
            for xsec in xsecs[0]:
                xsectot += xsec.value
            if xsectot.asNumber ( fb ) < 1e-10:
                self.pprint ( "WARNING no xsec??" )
            for i2,m2 in enumerate(ryvariable):
                if m2 > m1: ## we assume yvariable to be the daughter
                    continue
                self.M.masses[self.yvariable]=m2
                for pid_,m_ in self.M.masses.items():
                    if pid_ != self.yvariable and m_ < m2: ## make sure LSP remains the LSP
                        self.pprint ( f"WARNING: have to raise {pid_} from {m_} to {m2+1.}" )
                        oldmasses[pid_]=m_
                        self.M.masses[pid_]=m2 + 1.
                point = self.getPredictions ( False )
                llhds = point["llhd"]
                nllhds,nnonzeroes=0,0
                for mu,llhd in llhds.items():
                    nllhds+=len(llhd)

                self.pprint ( f"{i1}/{nxvariables}: m({namer.asciiName(self.xvariable)})={m1:.1f}, m2({namer.asciiName(self.yvariable)})={m2:.1f}, {len(llhds)} mu's, {nllhds} llhds." )
                point["mx"] = m1 
                point["my"] = m2
                masspoints.append ( point )
        return masspoints

def runThread ( threadid: int, rundir: str, M, xvariable, yvariable, mxvariable,
                myvariable, nevents: int, rxvariable, ryvariable, predictor, return_dict ):
    """ the method needed for parallelization to work """
    thread = LlhdThread ( threadid, rundir, M, xvariable, yvariable, mxvariable, myvariable, nevents,
                          predictor )
    newpoints = thread.run ( rxvariable, ryvariable )
    if return_dict != None:
        return_dict[threadid]=newpoints
    thread.clean()
    return newpoints

class LlhdScanner ( LoggerBase ):
    """ class that encapsulates a likelihood sweep """
    def __init__ ( self, protomodel, xvariable, yvariable, nproc, rundir : str,
                   dbpath : str = "official", select : str = "all",
                   do_srcombine : bool = True, skip_production : bool = False, 
                   dry_run : bool = False ):
        """
        :param rundir: the rundir
        :param dbpath: the database path
        :param skip_production: if possible, skip production, go to plotting
        :param dry_run: dont actually perform the actions
        """
        super ( LlhdScanner, self ).__init__ ( 0 )
        self.dry_run = dry_run
        self.rundirarg = rundir
        self.rundir = setup( rundir )
        self.M = protomodel
        self.xvariable = xvariable
        self.yvariable = yvariable
        self.nproc = nproc
        self.skip_production = skip_production
        self.predictor = Predictor ( 0, dbpath=dbpath, 
                select=select, do_srcombine = do_srcombine )
        self.pprint ( f"self.predictor = Predictor ( 0, dbpath='{dbpath}', select='{select}', do_srcombine = {do_srcombine} )" )

    def describeRange ( self, r ):
        """ describe range r in a string """
        if len(r)==0:
            return ""
        if len(r)==1:
            return f"{r[0]:.2f}"
        if len(r)==2:
            return f"{r[0]:.2f},{r[1]:.2f}"
        return f"{r[0]:.2f},{r[1]:.2f} ... {r[-1]:.2f} -> {len(r)} points"

    def getMassPoints ( self, rxvariable, ryvariable ):
        """ run for the given mass points
        :param rxvariable: list of masses for xvariable
        :param ryvariable: list of masses for yvariable
        :returns: masspoints
        """
        if self.dry_run:
            self.pprint ( f"dry_run. stopping here" )
            sys.exit()
            self.pprint ( f"dry_run. would run for xvariable={rxvariable}" )
            self.pprint ( f"yvariable={ryvariable}" )
            sys.exit()
        if self.nproc == 1:
            return runThread ( 0, self.rundir, self.M, self.xvariable, self.yvariable, \
                               self.mxvariable, self.myvariable, self.nevents, rxvariable, ryvariable,
                               self.predictor, None )
        chunkedRxvariable = [ list(rxvariable[i::self.nproc]) for i in range(self.nproc) ]
        processes = []
        manager = multiprocessing.Manager()
        return_dict=manager.dict()
        # print ( "chunked", chunkedRxvariable )
        for ctr,chunk in enumerate(chunkedRxvariable):
            self.M.walkerid = 2000+ctr
            p = multiprocessing.Process ( target = runThread, args = ( ctr, self.rundir, self.M, self.xvariable, self.yvariable, self.mxvariable, self.myvariable, self.nevents, chunk, ryvariable, self.predictor, return_dict ) )
            p.start()
            processes.append ( p )

        for p in processes:
            p.join()
        masspoints = []
        hasStored=set()
        for k,v in return_dict.items():
            for mp in v:
                key=(mp["mx"],mp["my"])
                if key in hasStored:
                    continue
                hasStored.add ( key )
                masspoints.append ( mp )
        return masspoints

    def scanLikelihoodFor ( self, range1 : Dict, range2 : Dict,
                            nevents : int, topo : str, output : str ):
        """ plot the likelihoods as a function of xvariable and yvariable

        :param range1: dictionary for range1 with min, max, dm
        :param range2: dictionary for range1 with min, max, dm
        :param output: prefix for output file [mp]
        """
        self.nevents = nevents
        xvariable = self.xvariable
        yvariable = self.yvariable
        if yvariable != self.M.LSP:
            print ( f"[llhdScanner] we currently assume yvariable to be the LSP, but it is {yvariable}" )
        picklefile = f"{output}{namer.asciiName(xvariable)}{namer.asciiName(yvariable)}.pcl"
        if os.path.exists ( picklefile ) and self.skip_production:
            print ( f"[llhdScanner] we were asked to skip production: {picklefile} exists." )
            return
        import numpy
        c = Combiner()
        anaIds = c.getAnaIdsWithPids ( self.M.bestCombo, [ xvariable, yvariable ] )
        ## mass range for xvariable
        self.mxvariable = self.M.masses[xvariable]
        self.myvariable = self.M.masses[yvariable]
        
        rxvariable = numpy.arange ( range1["min"], range1["max"]+1e-8, range1["dm"] )
        ryvariable = numpy.arange ( range2["min"], range2["max"]+1e-8, range2["dm"] )
        print ( f"[llhdScanner] range for {namer.asciiName(xvariable)}: {self.describeRange( rxvariable )}" )
        print ( f"[llhdScanner] range for {namer.asciiName(yvariable)}: {self.describeRange( ryvariable )}" )
        print ( f"[llhdScanner] total {len(rxvariable)*len(ryvariable)} points, {nevents} events for {topo}" )
        self.M.createNewSLHAFileName ( prefix=f"llhd{xvariable}" )
        #self.M.initializePredictor()
        self.predictor.filterForTopos ( topo )
        self.M.walkerid = 2000

        thread0 = LlhdThread ( 0, self.rundir, self.M, self.xvariable, self.yvariable, \
                               self.mxvariable, self.myvariable, self.nevents, self.predictor )
        point = thread0.getPredictions ( False )
        llhds = point["llhd"]
        critics = point["critic"]
        thread0.clean()
        self.pprint ( f"protomodel point: m1({namer.asciiName(self.xvariable)})={self.mxvariable:.2f}, m2({namer.asciiName(self.yvariable)})={self.myvariable:.2f}, {len(llhds)} llhds" )
        point [ "mx" ] = self.mxvariable
        point [ "my" ] = self.myvariable
        masspoints = [ point ]

        if True:
            ## freeze out all other particles
            for pid_,m_ in self.M.masses.items():
                if pid_ not in [ self.xvariable, self.yvariable ]:
                    self.M.masses[pid_]=1e6

        newpoints = self.getMassPoints ( rxvariable, ryvariable )
        masspoints += newpoints
        import pickle
        if os.path.exists ( picklefile ):
            subprocess.getoutput ( f"cp {picklefile} {picklefile}.old" )
        self.pprint ( f"now saving to {picklefile}" )
        f=open( picklefile ,"wb" )
        mydict = { "masspoints": masspoints, "mxvariable": self.mxvariable,
                   "myvariable": self.myvariable, "nevents": nevents, "topo": topo,
                   "timestamp": time.asctime(), "xvariable": xvariable, "yvariable": yvariable,
                   "model": self.M.dict() }
        pickle.dump ( mydict, f )
        f.close()
        self.M.delCurrentSLHA()

    def overrideWithDefaults ( self, args ):
        topo = { 1000005: "T2bb",1000006: "T2tt", 2000006: "T2tt", 1000021: "T1", \
                 1000023: "electroweakinos_offshell", 
                 1000024: "electroweakinos_offshell,T1",
                 1000001: "T2",  1000002: "T2", 1000003: "T2", 1000004: "T2" }
        ### make the LSP scan depend on the mother
        if args.topo == None:
            args.topo = topo[args.xvariable]
        self.mxvariable = self.M.masses[self.xvariable]
        self.myvariable = self.M.masses[self.yvariable]
        nbinsx, nbinsy = 20, 20 # how many bins do we want per dimension
        if args.min1 == None:
            args.min1 = self.mxvariable*.6
        if args.max1 == None:
            args.max1 = self.mxvariable*1.7
        if args.deltam1 == None:
            args.deltam1 = ( args.max1 - args.min1 ) / nbinsx
        if args.min2 == None:
            args.min2 = max ( self.myvariable*.6 - 10., 1. )
        if args.max2 == None:
            args.max2 = self.myvariable*1.9 + 10.
        if args.deltam2 == None:
            args.deltam2 = ( args.max2 - args.min2 ) / nbinsy
        return args

def main ():
    import argparse
    argparser = argparse.ArgumentParser(
            description='perform likelhood scans')
    argparser.add_argument ( '-n', '--number',
            help='which hiscore to plot [0]',
            type=int, default=0 )
    argparser.add_argument ( '-x', '--xvariable',
            help='variable for x axis, e.g. 1000006 or "Xt" [Xt]',
            type=str, default='Xt' )
    argparser.add_argument ( '-y', '--yvariable',
            help='variable for y axis, e.g. 1000022 or "(Xt,Xt)", [X1Z]',
            type=str, default="X1Z" )
    argparser.add_argument ( '-P', '--nproc',
            help='number of process to run in parallel. zero is autodetect. Negative numbers are added to autodetect [0]',
            type=int, default=0 )
    argparser.add_argument ( '-m1', '--min1',
            help='minimum mass of xvariable [None]',
            type=float, default=None )
    argparser.add_argument ( '-M1', '--max1',
            help='maximum mass of xvariable [None]',
            type=float, default=None )
    argparser.add_argument ( '-d1', '--deltam1',
            help='delta m of xvariable [None]',
            type=float, default=None )
    argparser.add_argument ( '-m2', '--min2',
            help='minimum mass of yvariable [None]',
            type=float, default=None )
    argparser.add_argument ( '-M2', '--max2',
            help='maximum mass of yvariable [None]',
            type=float, default=None )
    argparser.add_argument ( '-d2', '--deltam2',
            help='delta m of yvariable [None]',
            type=float, default=None )
    argparser.add_argument ( '-t', '--topo',
            help='topology [None]',
            type=str, default=None )
    argparser.add_argument ( '-R', '--rundir',
            help='override the default rundir [None]',
            type=str, default=None )
    argparser.add_argument ( '-e', '--nevents',
            help='number of events [50000]',
            type=int, default=50000 )
    argparser.add_argument ( '-H', '--hiscores',
            help='hiscore file to draw from [<rundir>/hiscores.dict]',
            type=str, default="default" )
    argparser.add_argument ( '-D', '--draw',
            help='also perform the plotting, ie call plotLlhds',
            action='store_true' )
    argparser.add_argument ( '--dry_run',
            help='just tell us what you would be doing, dont actually do it',
            action='store_true' )
    argparser.add_argument ( '-v', '--verbosity',
            help='verbosity -- debug, info, warn, err [info]',
            type=str, default="info" )
    argparser.add_argument ( '-o', '--output',
            help="prefix for output file [llhd]",
            type=str, default="llhd" )
    argparser.add_argument ( '-u', '--uploadTo',
            help="where do we upload to, on smodels.github.io [latest]",
            type=str, default="latest" )
    argparser.add_argument ( '-s', '--select',
            help="what do we select for [all]",
            type=str, default="all" )
    argparser.add_argument ( '--dbpath',
            help="path to database [official]",
            type=str, default="official" )
    argparser.add_argument ( '-c', '--do_srcombine',
            help='do_srcombine', action='store_true' )
    argparser.add_argument ( '-S', '--skip_production',
            help='if possible, skip production', action='store_true' )
    args = argparser.parse_args()
    rundir = setup( args.rundir )
    nproc = args.nproc
    if nproc < 1:
        nproc = nCPUs() + nproc
    if args.hiscores == "default":
        args.hiscores = f"{rundir}/hiscores.dict"
    from ptools.hiscoreTools import fetchHiscoresObj
    hi = fetchHiscoresObj ( args.hiscores, None, args.dbpath )
    protomodel = hi.hiscores[0]
    #self.pprint ( f"fetched {protomodel} from {args.hiscores}" )

    xvariables = [ namer.pid ( args.xvariable ) ]
    if args.xvariable == 0:
        xvariables = findPids( rundir )
    for xvariable in xvariables:
        yvariable = namer.pid ( args.yvariable )
        scanner = LlhdScanner( protomodel, xvariable, yvariable, nproc, rundir, 
                dbpath = args.dbpath, select = args.select, 
                do_srcombine = args.do_srcombine, 
                skip_production = args.skip_production,
                dry_run = args.dry_run )
        args.xvariable = xvariable
        args = scanner.overrideWithDefaults ( args )
        range1 = { "min": args.min1, "max": args.max1, "dm": args.deltam1 }
        range2 = { "min": args.min2, "max": args.max2, "dm": args.deltam2 }
        scanner.scanLikelihoodFor ( range1, range2, args.nevents, args.topo, 
                args.output )
        if args.draw:
            verbose = args.verbosity
            copy = True
            max_anas = 5
            interactive = False
            drawtimestamp = True
            compress = False
            upload = args.uploadTo
            plot = plotLlhds.LlhdPlot ( xvariable, yvariable, verbose, copy, max_anas,
                  interactive, drawtimestamp, compress, rundir, upload, args.dbpath )
            syv = "_"+namer.asciiName(yvariable)
            if syv == "_X1Z":
                syv = ""
            scriptfilename = f"llhdPlot_{namer.asciiName(xvariable)}{syv}.py"
            with open ( scriptfilename, "wt" ) as f:
                print ( f"[llhdScanner] created llhdPlotScript.py" )
                f.write ( "#!/usr/bin/env python3\n\n" )
                f.write ( "import sys\n" )
                f.write ( "interactive=False\n" )
                f.write ( "if '-i' in sys.argv:\n" )
                f.write ( "    interactive=True\n" )
                f.write ( "from plotting import plotLlhds\n" )
                f.write ( f"plot = plotLlhds.LlhdPlot ( xvariable={xvariable}, yvariable={args.yvariable}, verbose='{verbose}', copy={copy},\n" )
                f.write ( f"    max_anas={max_anas}, interactive=interactive, drawtimestamp={drawtimestamp}, compress={compress},\n" )
                f.write ( f"    rundir='{rundir}',\n" )
                f.write ( f"    upload='{upload}', dbpath='{args.dbpath}' )\n" )
                f.write ( f"plot.plot()\n" )
                f.write ( f"if '-s' in sys.argv:\n" )
                f.write ( f"    plot.show()\n" )
                f.close()
                os.chmod ( scriptfilename, 0o755 )
            plot.plot()

if __name__ == "__main__":
    main()
