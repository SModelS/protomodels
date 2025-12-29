#!/usr/bin/env python3

""" a second attempt at plotting likelihoods """

import os
from base.loggerbase import LoggerBase
import matplotlib.pyplot as plt

class NLLPlotter ( LoggerBase ):
    """ our second generation 2d nll plotter """

    def __init__ ( self, args : dict ):
        super ( NLLPlotter, self ).__init__ ( "nll" )
        self.args = args
        if self.args["inputfile"]==None:
            self.findInputFile()
        self.readInputFile()
        self.plot()
        if self.args["show"]:
            self.show()
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

    def getCriticList ( self, critic_type : str = "ul" ) -> list:
        """
        get the list of critic outputs per point
        :param critic_type: one of: ul, llhd

        :returns: list of dictionaries, "mx", "my", "passes" as keys
        """
        assert critic_type in [ "ul", "llhd" ], \
             f"critic type should be one of: ul, llhd"
        ret = []
        for masspoint in self.data["masspoints"]:
            mx, my = masspoint["mx"], masspoint["my"]
            passes = masspoint["critic"][critic_type]["passes"]
            tmp = { "mx": mx, "my": my, "passes": passes }
            ret.append ( tmp )
        return ret

    def plot ( self ):
        critic_points = self.getCriticList()
        self.plotBooleanMap ( critic_points )
        self.pprint ( f"plotting {self.args['inputfile']} -> {self.outputfile}" )

    def plotBooleanMap ( self, points : list[dict] ):
        """ given a list of dicts, draw the contour
        :param points: a list of dictionaries, "mx", "my", "passes"
        """
        import numpy as np
        from scipy.interpolate import griddata

        # Extract arrays
        x = np.array([d["mx"] for d in points])
        y = np.array([d["my"] for d in points])
        z = np.array([d["passes"] for d in points], dtype=float)

        # Create interpolation grid
        xi = np.linspace(x.min(), x.max(), 200)
        yi = np.linspace(y.min(), y.max(), 200)
        Xi, Yi = np.meshgrid(xi, yi)

        # Interpolate boolean field
        Zi = griddata((x, y), z, (Xi, Yi), method="linear")

        # Draw contour where passes == True
        plt.contour(Xi, Yi, Zi, levels=[0.5])
        plt.scatter(x, y, c=z, cmap="coolwarm", s=30)
        plt.xlabel("mx")
        plt.ylabel("my")
        plt.title("Contour of passes == True")
        self.savefig()

    def savefig ( self ):
        """ save the figure to file """
        figname = self.args["outputfile"]
        if "@@I@@" in figname:
            from pathlib import Path
            ifile = Path ( self.args["inputfile"] ).stem
            figname = figname.replace("@@I@@",ifile)
        self.outputfile = figname
        # self.pprint ( f"saving to {figname}" )
        plt.savefig ( figname )

    def show ( self ):
        from smodels_utils.plotting.mpkitty import timg
        timg ( self.outputfile )

    def interact ( self ):
        import sys, IPython; IPython.embed( colors = "neutral" )


if __name__ == "__main__":
    import argparse
    argparser = argparse.ArgumentParser(
            description='plot likelihoods scans')
    argparser.add_argument ( '-i', '--inputfile',
            help='input file',
            type=str, default=None )
    argparser.add_argument ( '-o', '--outputfile',
            help='output file [@@I@@.png]',
            type=str, default="@@I@@.png" )
    argparser.add_argument ( '-I', '--interact',
            help='start interactive shell', action="store_true" )
    argparser.add_argument ( '-s', '--show',
            help='show image', action="store_true" )
    args = argparser.parse_args()
    plotter = NLLPlotter ( vars(args) )
