#!/usr/bin/env python3

""" a second attempt at plotting likelihoods """

import os
from base.loggerbase import LoggerBase
import matplotlib.pyplot as plt
from scipy.interpolate import griddata
import numpy as np
from ptools.sparticleNames import SParticleNames

namer = SParticleNames ( False )

class NLLPlotter ( LoggerBase ):
    """ our second generation 2d nll plotter """

    def __init__ ( self, args : dict ):
        super ( NLLPlotter, self ).__init__ ( "nll" )
        self.args = args
        self.outputfile = None
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

    def getCriticList ( self, critic_type : str = "both" ) -> list:
        """
        get the list of critic outputs per point
        :param critic_type: one of: ul, llhd, both

        :returns: list of dictionaries, "mx", "my", "passes" as keys
        """
        assert critic_type in [ "ul", "llhd", "both" ], \
             f"critic type should be one of: ul, llhd, both"
        ret = []
        for masspoint in self.data["masspoints"]:
            mx, my = masspoint["mx"], masspoint["my"]
            if critic_type == "both":
                p_ul = masspoint["critic"]["ul"]["passes"]
                p_llhd = masspoint["critic"]["llhd"]["passes"]
                passes = p_ul and p_llhd
            else:
                passes = masspoint["critic"][critic_type]["passes"]
            tmp = { "mx": mx, "my": my, "passes": passes }
            ret.append ( tmp )
        return ret

    def getNLLList ( self, anaid : str = "combined", 
                     removeDisallowed : bool = True ) -> list:
        """
        get the list of likelihoods for anaid per point
        :param anaid: e.g. ATLAS-SUSY-2018-06:EM6, or combined
        :param removeDisallowed: remove points not allowed by critic

        :returns: list of dictionaries, "mx", "my", "nll" as keys
        """
        ret = [] 
        # nlls, llhds = [], []
        for masspoint in self.data["masspoints"]:
            if removeDisallowed and masspoint["critic"]["ul"]["passes"]==False\
                    or masspoint["critic"]["llhd"]["passes"]==False:
                continue
            mx, my = masspoint["mx"], masspoint["my"]
            anlls = masspoint["nll"]
            for ssm, anas in anlls.items():
                if abs(ssm-1.)<1e-5 and anaid in anas:
                    nll = anas[anaid]
                    # nlls.append ( nll )
                    llhd = float ( np.exp ( - nll ) )
                    # llhds.append ( llhd )
                    # nll = masspoint["nll"][1.0][anaid]
                    tmp = { "mx": mx, "my": my, "nll": nll, "llhd": llhd }
                    ret.append ( tmp )
                    break
        return ret

    def normalize ( self, points : list[dict], how : str = "max_llhd" ) -> list[dict]:
        """ normalize the likelihoods given in points, in various ways
        :param how: one of: max_llhd

        :returns: normalized points
        """
        ret = []
        nlls = [ x["nll"] for x in points ]
        llhds = [ x["llhd"] for x in points ]
        min_nll = min ( nlls )
        max_llhd = max ( llhds )
        for p in points:
            tmp = { "mx": p["mx"], "my": p["my"] }
            tmp["dnll"]=p["nll"]-min_nll
            tmp["llhd_rel"]=p["llhd"]/max_llhd
            ret.append ( tmp )
        return ret

    def plotLikelihoodMass ( self, points : list[dict], options : dict ):
        """ plot the likelihood mass given in nll_points """
        defaults = { "text": False, "colors": [ "red", "darkred" ],
                     "label": "probability mass" }
        opts = defaults
        opts.update ( options )
        points = self.normalize ( points, "max_llhd" )
        # ---- input data ----
        # example: data = [{"mx": ..., "my": ..., "llhd": ...}, ...]
        xs = np.array([d["mx"] for d in points])
        ys = np.array([d["my"] for d in points])
        ll = np.array([d["llhd_rel"] for d in points])

        # ---- make a regular grid ----
        nx, ny = 200, 200
        xi = np.linspace(xs.min(), xs.max(), nx)
        yi = np.linspace(ys.min(), ys.max(), ny)
        X, Y = np.meshgrid(xi, yi)

        Z = griddata((xs, ys), ll, (X, Y), method="linear")
        Z = np.nan_to_num(Z, nan=0.0)

        # ---- convert likelihood to probability ----
        Z = np.maximum(Z, 0)
        Z /= Z.sum()

        # ---- find contour levels for given probability mass ----
        def contour_level_for_mass(Z, mass):
            z_sorted = np.sort(Z.ravel())[::-1]
            cumsum = np.cumsum(z_sorted)
            idx = np.searchsorted(cumsum, mass)
            return z_sorted[idx]

        ## 1 and 2 sigma
        levels = [ float ( 1-np.exp(-.5) ), 
                   float ( 1-np.exp(-2) ) ]
        level_1s = contour_level_for_mass(Z, levels[0] )
        level_2s = contour_level_for_mass(Z, levels[1] )

        plt.contourf( X, Y, Z,
            levels=[level_2s, level_1s,Z.max()],
            colors=opts["colors"], alpha=0.15 )
        # ---- plot ----
        # plt.figure(figsize=(6, 5))
        cs = plt.contour(X, Y, Z, levels=[level_2s, level_1s],
                    colors=opts["colors"], linewidths=2 )
        # Label contours
        fmt = {
            level_1s: f"{int(levels[0]*100):d}%",
            level_2s: f"{int(levels[1]*100):d}%"
        }
        plt.clabel(cs, cs.levels, inline=True, fmt=fmt, fontsize=10)
        # plt.scatter(xs, ys, s=5, c="k", alpha=0.3)
        from matplotlib.lines import Line2D
        legend_elements = [
            Line2D([0], [0], color=opts["colors"][1], lw=2, label=opts["label"] ),
#            Line2D([0], [0], color=opts["colors"][0], lw=2, label=opts["label"] ),
        ]
        for legend_element in legend_elements:
            self.handles.append ( legend_element )

        plt.scatter ( xs, ys, s=1 )
        if opts["text"]:
            for d in points:
                plt.text ( d["mx"], d["my"], f"{d['llhd_rel']:.1g}" )

    def plot ( self ):
        """
        critic_points = self.getCriticList( critic_type = "ul" )
        options = { "label": "excluded by ul" }
        self.plotBooleanMap ( critic_points, options )
        """
        critic_points = self.getCriticList( critic_type = "both" )
        options = { "c_area": "gray", "c_line": "dimgray", "hatches": "\\\\",
                    "label": "excluded by critic"  }
        self.plotBooleanMap ( critic_points, options )

        rmCritic = True
        anaid = "CMS-SUS-20-004:(comb):TChiHH"
        nll_points = self.getNLLList( anaid = anaid,
               removeDisallowed = rmCritic )
        options = { "text": True, "label": anaid }
        self.plotLikelihoodMass ( nll_points, options )

        # Existing scatter handles (from plt.scatter calls)
        # handles, labels = plt.gca().get_legend_handles_labels()

        # Add the area patch to the legend
        plt.legend( handles=self.handles, loc="best")
        plt.xlabel( rf"m$\left({namer.texName(self.data['meta']['xvariable'])}\right)$ [GeV]" )
        plt.ylabel( rf"m$\left({namer.texName(self.data['meta']['yvariable'])}\right)$ [GeV]" )
        plt.title ( "probability mass" )
        self.savefig()
        self.pprint ( f"plotting {self.args['inputfile']} -> {self.outputfile}" )

    def plotBooleanMap ( self, points : list[dict], options : dict  ):
        """ given a list of dicts, draw the contour
        :param points: a list of dictionaries, "mx", "my", "passes"
        """
        import numpy as np
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
        if self.outputfile == None:
            return
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
