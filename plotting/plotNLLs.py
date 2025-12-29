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
        self.handles = []
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
        critic_points = self.getCriticList( critic_type = "ul" )
        options = { "label": "excluded by ul" }
        self.plotBooleanMap ( critic_points, options )
        critic_points = self.getCriticList( critic_type = "llhd" )
        options = { "c_area": "gray", "c_line": "dimgray", "hatches": "\\\\",
                    "label": "excluded by llhd"  }
        self.plotBooleanMap ( critic_points, options )

        # Existing scatter handles (from plt.scatter calls)
        # handles, labels = plt.gca().get_legend_handles_labels()

        # Add the area patch to the legend
        plt.legend( handles=self.handles, loc="best")
        plt.title ( "llhd based critic" )
        self.savefig()
        self.pprint ( f"plotting {self.args['inputfile']} -> {self.outputfile}" )

    def plotBooleanMap ( self, points : list[dict], options : dict  ):
        """ given a list of dicts, draw the contour
        :param points: a list of dictionaries, "mx", "my", "passes"
        """
        import numpy as np
        from scipy.interpolate import griddata
        defaults = { "c_area": "gray", "c_line": "dimgray",
            "scatter": False, "xlabel": "mx", "ylabel": "my", "hatches": "////",
            "label": "excluded" }
        opts = defaults
        opts.update ( options ) 

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
        plt.contourf(Xi, Yi, Zi, levels=[-0.1,0.5], colors=[ opts["c_area"] ],
                     alpha=0.8,hatches = [ opts["hatches"] ] )
        plt.contour(Xi, Yi, Zi, levels=[0.5], colors=[ opts["c_line"] ] )
        # plt.scatter(x, y, c=z, cmap="coolwarm", s=30)
        # passes == False → red
        if opts["scatter"]:
            mask_true = z == 1
            mask_false = z == 0
            plt.scatter( x[mask_false], y[mask_false],
                color="red", edgecolor="black", s=40 )

            # passes == True → green
            plt.scatter( x[mask_true], y[mask_true],
                color="green", edgecolor="black", s=40 )
        plt.xlabel( opts["xlabel"] )
        plt.ylabel( opts["ylabel"] )
        
        if opts["label"] not in [ None, "" ]:
            import matplotlib.patches as mpatches
            # Legend entry for hatched grey area (passes == False region)
            false_area_patch = mpatches.Patch(
                facecolor= opts[ "c_area" ],
                edgecolor= opts [ "c_line" ],
                hatch= opts[ "hatches" ],
                label= opts[ "label" ]
            )
            self.handles.append ( false_area_patch )

    def savefig ( self ):
        """ save the figure to file """
        figname = self.args["outputfile"]
        if "@@I@@" in figname:
            from pathlib import Path
            ifile = Path ( self.args["inputfile"] ).stem
            figname = figname.replace("@@I@@",ifile)
        self.outputfile = figname
        from smodels_utils.helper.various import pngMetaInfo
        metadata = pngMetaInfo()
        metadata["Prod Commandline"] = self.data["meta"]["cmdline"]
        # self.pprint ( f"saving to {figname}" )
        plt.savefig ( figname, metadata = metadata )

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
