#!/usr/bin/env python3

""" a second attempt at plotting likelihoods """

import os
from base.loggerbase import LoggerBase

class NLLPlotter ( LoggerBase ):
    """ our second generation 2d nll plotter """

    def __init__ ( self, args : dict ):
        super ( NLLPlotter, self ).__init__ ( "nll" )
        self.args = args
        if self.args["inputfile"]==None:
            self.findInputFile()
        self.readInputFile()
        self.plot()
        if self.args["interact"]:
            self.interact()

    def findInputFile ( self ):
        """ no -I argument was given, so find an input file.
        currently we simply return the first in the list of 
        all matching files """
        import glob
        files = glob.glob ( "nll*.dict" )
        if len(files)==0:
            files = glob.glob ( "nll*.pcl" )
        if len(files)==0:
            self.pprint ( "no input files found" )
            sys.exit()
        if len(files)>1:
            self.pprint ( f"found multiple input files, will plot for {files[0]}" )
        self.args["inputfile"]=files[0]

    def readInputFile ( self ):
        """ read in the content of args["inputfile"] """
        ifile = self.args["inputfile"]
        if not os.path.exists ( ifile ):
            self.error ( "inputfile {ifile} does not exist" )
            sys.exit()
        if ifile.endswith ( ".dict" ):
            with open ( ifile, "rt" ) as f:
                self.data = eval ( f.read() )
                return
        import pickle
        with open ( ifile, "rb" ) as f:
            self.data = pickle.load ( f )

    def plot ( self ):
        self.pprint ( f"plotting {self.args['inputfile']}" )

    def interact ( self ):
        import sys, IPython; IPython.embed( colors = "neutral" )


if __name__ == "__main__":
    import argparse
    argparser = argparse.ArgumentParser(
            description='plot likelihoods scans')
    argparser.add_argument ( '-i', '--inputfile',
            help='input file',
            type=str, default=None )
    argparser.add_argument ( '-I', '--interact',
            help='start interactive shell', action="store_true" )
    args = argparser.parse_args()
    plotter = NLLPlotter ( vars(args) )
