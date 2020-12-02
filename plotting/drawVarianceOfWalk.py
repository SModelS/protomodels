#!/usr/bin/env python3

""" a plotting script that draws the evolution of the Ks, and their variances, as 
    function of the step """

import glob, pickle, sys, os, subprocess, argparse
import numpy as np
import matplotlib.pyplot as plt


class VariancePlotter:
    def __init__ ( self ):
        self.filenames = "history*.list"

    def pprint ( self, *args ):
        print ( "[drawVarianceOfWalk] %s" % " ".join(map(str,args)) )  

    def fetchData( self ):
        cmd = "rm -r history*list"
        subprocess.getoutput ( cmd )
        dirs = [ "twice10", "twice1", "twice2" ]
        for d in dirs:
            cmd = f"mkdir -p {d}"
            subprocess.getoutput ( cmd )
            cmd = f"rm -r {d}/history*list"
            subprocess.getoutput ( cmd )
            self.pprint ( f"fetching data from clip, {d}" )
            cmd = f"scp clip-login-1:/scratch-cbe/users/wolfgan.waltenberger/rundir.{d}/history\*.list {d}/"
            subprocess.getoutput ( cmd )
            number = d.replace("twice","")
            number=int(number)
            cmd = f'mcp "{d}/history*.list" "history#1_{number:02d}.list"'
            subprocess.getoutput ( cmd )

    def getData ( self ):
        files = glob.glob ( self.filenames )
        Kvalues = {}
        for f in files:
            h = open ( f, "rt" )
            walkerid = f.replace( ".list", "" ).replace( "history", "" )
            walkerid = int ( walkerid )
            try:
                txt = h.read()
                txt = txt.replace("nan",'"nan"')
                # print ( "txt >>%s<< {{%s}}" % ( txt[-5:], txt[-1] ) )
                if txt[-2]==",":
                    txt = txt[:-2]+"]\n"
                if "history recording" in txt[-80:]:
                    txt = txt[:-10]+"\n]\n"
                D=eval( txt )
                Kvalues[walkerid]=D
            except Exception as e:
                print ( "[drawVarianceOfWalk] could not read %s: %s. skipping" % \
                        ( f, e ) )
            h.close()
        self.data = Kvalues

    def interact ( self ):
        import IPython
        IPython.embed( using=False )

    def loadData ( self ):
        if not os.path.exists ( "var.pcl" ):
            return None
        f = open  ( "var.pcl", "rb" )
        self.data = pickle.load ( f )
        f.close()

    def storeData ( self ):
        """ store it all in a pickle file """
        f = open  ( "var.pcl", "wb" )
        pickle.dump ( self.data, f )
        f.close()

    def getMaxSteps ( self, data ):
        """ find out the max step size """
        nmax = 0
        maxes = []
        for d,v in data.items():
            walkermax = -5
            for s in v:
                step = s["step"]
                if step > nmax:
                    nmax = step
                if step > walkermax:
                    walkermax = step
            maxes.append ( walkermax )
        ret = min(maxes)
        if ret > 1999:
            ret = 1999
        return ret

    def getKsPerStep ( self, var="K" ):
        """ from the data, get the Ks of all walkers as a function
            of the step """
        nsteps = self.getMaxSteps ( self.data )
        ret = {}
        for walker,walk in self.data.items():
            for wstep in walk:
                nstep = wstep["step"]
                K = wstep[var]
                if K == None:
                    continue
                if not nstep in ret:
                    ret[nstep]=[]
                ret[nstep].append ( K )
        return ret

    def getMassOfPidPerStep ( self, pid=1000006 ):
        """ from the data, get the Ks of all walkers as a function
            of the step """
        nsteps = self.getMaxSteps ( self.data )
        ret = {}
        for walker,walk in self.data.items():
            for wstep in walk:
                nstep = wstep["step"]
                K = None
                if pid in wstep["masses"]:
                    K = float(wstep["masses"][pid])
                if K == None:
                    continue
                if not nstep in ret:
                    ret[nstep]=[]
                ret[nstep].append ( K )
        return ret

    def draw ( self, var = "m1000006", drawMax=True, output = "@VAR.png" ):
        """ draw the evolution of var
        :param var: draw variable "var". If starts with "m", draw mass of pid
        :param drawMax: if True, also draw the max value per step
        """
        output = output.replace( "@VAR", var )
        if var in [ "K", "Z" ]:
            Ks = self.getKsPerStep() ## this is per step
        if var.startswith("m"):
            Ks = self.getMassOfPidPerStep ( int(var[1:]) )
            
        keys = list ( Ks.keys() )
        keys.sort()
        means, maxs = [], []
        lastMax = -5
        for step in keys:
            KsInStep = Ks[step] ## the K values in this step
            means.append ( np.mean ( KsInStep ) )
            thisMax =  np.max ( KsInStep )
            if thisMax > lastMax:
                lastMax = thisMax 
            maxs.append ( lastMax )

        nvalues = [ len(Ks[x]) for x in keys ]
        def minus ( x ):
            if x>50:
                return x-50
            return 0
        avgedmeans = [ np.mean(means[minus(x):x+50]) for x in range(len(keys)) ]
        # avgedmaxs = [ np.max(maxs[minus(x):x+50]) for x in range(len(keys)) ]
        # print ( "avgedmeans", avgedmeans )
        plt.plot ( keys, means, c="orange" )
        if drawMax:
            plt.plot ( keys, maxs, c="red" )
        plt.plot ( keys, avgedmeans, c="black" )
        plt.title ( "evolution of $K$, for %d walkers" % len(self.data) )
        plt.xlabel ( "step" )
        plt.ylabel ( var )
        # plt.plot ( keys, nvalues, c="orange" )
        # outputfile = "var.png"
        outputfile = output
        self.pprint ( f"print to {outputfile}" )
        plt.savefig ( outputfile )

if __name__ == "__main__":
    import argparse
    argparser = argparse.ArgumentParser(description="tool to plot how K develops over the steps")
    argparser.add_argument ( '-f', '--fetch',
            help='fetch the history files', action="store_true" )
    argparser.add_argument ( '-i', '--interact',
            help='enter interactive mode', action="store_true" )
    argparser.add_argument ( '-m', '--drawmax',
            help='draw max value', action="store_true" )
    argparser.add_argument ( '-v', '--variable',
            help='variable to draw, e.g. K, m1000006 [m1000006]',
            type=str, default='m1000006' )
    argparser.add_argument ( '-o', '--output',
            help='output file name, replace @VAR with variable [@VAR.png]',
            type=str, default='@VAR.png' )
    args=argparser.parse_args()

    plotter = VariancePlotter()
    if args.fetch: 
        plotter.fetchData()
        plotter.getData()
        plotter.storeData()
    else:
        plotter.loadData()
    plotter.draw( args.variable, args.drawmax, args.output )
    if args.interact:
        plotter.interact()
