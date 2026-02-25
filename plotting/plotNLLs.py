#!/usr/bin/env python3

""" a second attempt at plotting likelihoods """

import os, copy, sys
from base.loggerbase import LoggerBase
import matplotlib.pyplot as plt
from scipy.interpolate import griddata
import numpy as np
from ptools.sparticleNames import SParticleNames
from typing import Union

namer = SParticleNames ( False )

class NLLPlotter ( LoggerBase ):
    """ our second generation 2d nll plotter """

    def __init__ ( self, args : dict ):
        super ( NLLPlotter, self ).__init__ ( "nll" )
        self.args = args
        options = eval ( args["options"] )
        self.setOptions ( options )
        self.outputfile = None
        self.handles = []
        if self.args["inputfile"]==None:
            self.findInputFile()

    def run ( self ):
        self.readInputFile()
        self.plot()
        if self.args["show"]:
            self.show()
        if self.args["interact"]:
            self.interact()

    def setOptions ( self, options : dict ):
        """ given the options string at command line, draw options """
        self.options = {}
        for name,value in options.items():
            if not ":" in name:
                self.options[ name ] = value
                continue
            tnames = name.split(":")
            if len(tnames)!=2:
                self.warn ( f"option {name} unknown" )
                continue
            if not tnames[0] in self.options:
                self.options[ tnames[0] ] = {}
            self.options[ tnames[0] ][ tnames[1] ] = value

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
            self.error ( f"inputfile {ifile} does not exist" )
            sys.exit()
        if ifile.endswith ( ".dict" ):
            with open ( ifile, "rt" ) as f:
                self.data = eval ( f.read() )
                return
        import pickle
        with open ( ifile, "rb" ) as f:
            self.data = pickle.load ( f )

    def getCriticList ( self, critic_type : str = "both" ) -> list[dict]:
        """
        get the list of critic outputs per point
        :param critic_type: one of: ul, llhd, both

        :returns: list of dictionaries: example:
        [{ "x": .. , "y": ..., "robs": ..., "rexp": ..., "passes": ... },...]
        """
        assert critic_type in [ "ul", "llhd", "both" ], \
             f"critic type should be one of: ul, llhd, both"
        ret = []
        selection = self.getOptions ( "selection" )
        for parameterpoint in self.data["parameterpoints"]:
            mx, my = parameterpoint["x"], parameterpoint["y"]
            if my > selection["ymax"]:
                continue
            if mx > selection["xmax"]:
                continue
            tmp = { "x": mx, "y": my }
            if critic_type == "both":
                passes = False
                if parameterpoint["critic"]!=None:
                    p_ul = parameterpoint["critic"]["ul"]["passes"]
                    p_llhd = parameterpoint["critic"]["llhd"]["passes"]
                    passes = p_ul and p_llhd
                    tmp["robs"] = parameterpoint["critic"]["llhd"]["robs"]
                    tmp["rexp"] = parameterpoint["critic"]["llhd"]["rexp"]
                tmp["passes"]=passes
            else:
                passes = parameterpoint["critic"][critic_type]["passes"]
                tmp["robs"] = parameterpoint["critic"][critic_type]["robs"]
                if "rexp" in parameterpoint["critic"][critic_type]:
                    tmp["rexp"] = parameterpoint["critic"][critic_type]["rexp"]
                tmp["passes"]=passes
            ret.append ( tmp )
        return ret

    def findAnaId ( self, anaid : str, anas : list[dict] ) -> Union[bool,str]:
        """ see if you can find anaid in anas
        :param anaid: e.g. CMS-EXO-20-004:(comb):TChiISR, CMS-EXO-20-004:(comb)
        :param anas: list of anaids with txnames
        :returns: the anaid that it identified, false if none found
        """
        if anaid in anas:
            return anaid
        for contender, nll in anas.items():
            if anaid in contender:
                return contender
        return False

    def getNLLList ( self, anaid : str = "joint",
                     removeDisallowed : bool = True ) -> list:
        """
        get the list of likelihoods for anaid per point
        :param anaid: e.g. CMS-EXO-20-004:(comb):TChiISR, CMS-EXO-20-004:(comb), or joint
        :param removeDisallowed: remove points not allowed by critic

        :returns: list of dictionaries, "x", "y", "nll" as keys
        """
        ret = []
        options = self.getOptions("selection")
        for parameterpoint in self.data["parameterpoints"]:
            if parameterpoint["critic"]==None:
                continue
            if removeDisallowed and parameterpoint["critic"]["ul"]["passes"]==False\
                    or parameterpoint["critic"]["llhd"]["passes"]==False:
                continue
            mx, my = parameterpoint["x"], parameterpoint["y"]
            if my > options["ymax"]:
                continue
            if mx > options["xmax"]:
                continue
            anlls = parameterpoint["nll"]
            for ssm, anas in anlls.items():
                if abs(ssm-1.)<1e-5:
                    myId = self.findAnaId ( anaid, anas )
                    if myId != False:
                        nll = anas[myId]
                        llhd = float ( np.exp ( - nll ) )
                        tmp = { "fullid": myId, "x": mx, "y": my, "nll": nll, "llhd": llhd }
                        ret.append ( tmp )
                        break
        self.pprint ( f"getNLLList: returning {len(ret)}/{len(self.data['parameterpoints'])} points for {anaid}" )
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
            tmp = { "x": p["x"], "y": p["y"] }
            tmp["dnll"]=p["nll"]-min_nll
            tmp["llhd_rel"]=p["llhd"]
            if max_llhd > 0.:
                tmp["llhd_rel"]=p["llhd"]/max_llhd
            ret.append ( tmp )
        return ret

    def plotProbabilityMass ( self, points : list[dict], options : dict,
                             mask_data : list[dict] = [] ) -> bool:
        """ plot the likelihood mass given in nll_points, using a KDE
        :param mask_data: list of dictionaries as we get it for the critic, e.g.
        [{"x": ..., "y": ..., "passes": True/False}, ... ]
        :returns: true if successful, else false
        """
        defaults = { "points": False, "text": False, "colors": ( "red", "darkred" ),
                     "label": "probability mass", "nlevels": 2 }
        opts = defaults
        opts.update ( options )
        points = self.normalize ( points, "max_llhd" )
        #print ( f"@@0 for {options['label']}:" )
        #for p in points[:3]:
        #    print ( f"@@1 {p}" )
        # ---- input data ----
        # example: points = [{"x": ..., "y": ..., "llhd_rel": ...}, ...]
        xs = np.array([d["x"] for d in points])
        ys = np.array([d["y"] for d in points])
        ll = np.array([d["llhd_rel"] for d in points])

        if len(points)<4:
            self.warn ( f"plotProbabilityMass, for {options['label']} we have {len(points)} points, thats too few. not plotting contours." )
        else:
            # ---- convert likelihood to probability ----
            #Z = np.maximum(Z, 0)
            # KDE
            kpoints = np.vstack([xs, ys])
            from scipy.stats import gaussian_kde
            kde = gaussian_kde(kpoints,weights=ll)

            # Grid
            nx, ny = 100, 100
            xi = np.linspace(xs.min(), xs.max(), nx)
            yi = np.linspace(ys.min(), ys.max(), ny)
            X, Y = np.meshgrid(xi, yi)

            grid_points = np.vstack([X.ravel(), Y.ravel()])
            Z = kde(grid_points).reshape(X.shape)

            from scipy.spatial import ConvexHull, Delaunay
            ## confine to the convex hull
            tri = Delaunay(np.column_stack([xs, ys]))
            hull_mask = tri.find_simplex(np.column_stack([X.ravel(), Y.ravel()])) >= 0
            hull_mask = hull_mask.reshape(X.shape)

            if len(mask_data)==0:
                Z = np.ma.array(Z, mask=~hull_mask )
            else:
                # take out all excluded by critic
                c_xs = np.array([d["x"] for d in mask_data])
                c_ys = np.array([d["y"] for d in mask_data])
                c_valid = np.array([d["passes"] for d in mask_data], dtype=bool)

                valid_grid = griddata( (c_xs, c_ys), c_valid.astype(float), (X, Y),
                    method="nearest").astype(bool)
                combined_mask = hull_mask & valid_grid
                Z = np.ma.array(Z, mask=~combined_mask )
            Z /= Z.sum()

            # ---- find contour levels for given probability mass ----
            def contour_level_for_mass(Z, mass):
                z_sorted = np.sort(Z.ravel())[::-1]
                cumsum = np.cumsum(z_sorted)
                idx = np.searchsorted(cumsum, mass)
                return z_sorted[idx]

            ## 1 and 2 sigma
            vlevels = [ float ( 1-np.exp(-.5) ),
                       float ( 1-np.exp(-2) ) ]
            level_1s = contour_level_for_mass(Z, vlevels[0] )
            level_2s = contour_level_for_mass(Z, vlevels[1] )

            colors = opts["colors"]
            #colors = ( "white", colors[1] )
            levels = [ level_2s, level_1s ]
            if opts["nlevels"]==1:
                levels = [ level_1s ]
            plt.contourf( X, Y, Z,
                levels=levels+[Z.max()],
                colors=colors, alpha=0.15 )
            # ---- plot ----
            # plt.figure(figsize=(6, 5))
            cs = plt.contour(X, Y, Z, levels=levels,
                             colors=opts["colors"], linewidths=2 )
            # Label contours
            fmt = {
                level_1s: f"{int(vlevels[0]*100):d}%",
                level_2s: f"{int(vlevels[1]*100):d}%"
            }
            plt.clabel(cs, cs.levels, inline=True, fmt=fmt, fontsize=10)
            from matplotlib.lines import Line2D
            legend_elements = [
                Line2D([0], [0], color=opts["colors"][1], lw=2, label=opts["label"] ),
    #            Line2D([0], [0], color=opts["colors"][0], lw=2, label=opts["label"] ),
            ]
            for legend_element in legend_elements:
                self.handles.append ( legend_element )

        if opts["text"] or opts["points"]:
            pointsize,fontsize=1,8
            if len(points)<10:
                pointsize,fontsize=3,10
            if opts["text"]:
                llhd_rel_min = .01
                for d in points:
                    if d["llhd_rel"] > llhd_rel_min:
                        plt.text ( d["x"], d["y"], f"{d['llhd_rel']:.1g}",
                                   fontsize=fontsize )
            else:
                pointsize=2
                plt.scatter ( xs, ys, c="black", s=pointsize )
        return len(points)>3

    def getAnaIds ( self, comb_only : bool = False,
                    dropTxname : bool = False ) -> set:
        """ get a set of all analysis ids that i can procure that have
        entries for ssm==1.0
        :param comb_only: if true, then return only (comb) ana ids
        :param dropTxname: if true, drop the txnames

        :returns: set of analysis ids
        """
        ret = set()
        for parameterpoint in self.data["parameterpoints"]:
            anlls = parameterpoint["nll"]
            for ssm, anas in anlls.items():
                for anaid,nll in anas.items():
                    if comb_only and not "(comb)" in anaid:
                        continue
                    if dropTxname:
                        p1 = anaid.find(":T")
                        anaid = anaid[:p1]
                    ret.add ( anaid )
        return ret

    def combineNLLs ( self, for_combination : dict[list[dict]] )-> list[dict]:
        """ given a dictionary of analyes and mass points, combined their NLLs
        into a joint NLL
        :param for_combination: e.g. { "CMS-SUS-20-004:(combined)": [ parameterpoints ] }
        parameterpoints is a .e.g [ { "y": ... , "x": ..., "nll": ... }, ... ]

        :returns: list of parameterpoints with combined NLLs
        """
        points = {}
        def getHash ( mx : float, my : float ):
            return int ( round(mx,4)*1e10+round(my,4)*1e5 )
        anaids = set( for_combination.keys() )
        for anaid, parameterpoints in for_combination.items():

            for parameterpoint in parameterpoints:
                h = getHash ( parameterpoint["x"], parameterpoint["y"] )
                if not h in points:
                    points[h]={}
                points[h][anaid] = parameterpoint
        ret = []
        for h,point in points.items():
            comb_point =  { "nll": 0, "llhd": 1, "anas": [] }
            for anaid,values in point.items():
                comb_point["x"]=values["x"]
                comb_point["y"]=values["y"]
                comb_point["nll"]+=values["nll"]
                comb_point["llhd"]*=values["llhd"]
                comb_point["anas"].append ( values["fullid"] )
            if len(point)<len(anaids):
                comb_point["nll"]=float("inf")
                comb_point["llhd"]=0.
            ret.append ( comb_point )
        # print ( f"@@0 ret {ret}" )
        return ret

    def getOptions( self, which : str ) -> dict:
        """ get the options for builder or critic
        :param which: one of: builder, critic, combo, all, anaids, rmCritic,
        selection
        """
        if which == "all":
            ret = {}
            for i in [ "builder", "critic", "combo", "anaids", "rmCritic",
            "selection" ]:
                tmp = self.getOptions ( i )
                if type(tmp) != dict:
                    ret[i]=tmp
                else:
                    for k,v in tmp.items():
                        ret[f"{i}:{k}"]=v
            return ret

        if which == "builder":
            ret = { "text": False, "nlevels": 1 }
            if "builder" in self.options:
                ret.update ( self.options["builder"] )
            return ret
        if which == "selection":
            ret = { "ymax": float("inf"), "xmax": float("inf") }
            if "selection" in self.options:
                ret.update ( self.options["selection"] )
            return ret

        if which == "critic":
            ret = { "c_area": "gray", "c_line": "dimgray", "hatches": "\\\\",
                    "label": "excluded", "scatter": False,
                    "text": False }
            if "critic" in self.options:
                ret.update ( self.options["critic"] )
            return ret

        if which == "combo":
            ret = { "colors": ( "0.20", "black" ), "label": "joint posterior",
                    "nlevels": 1 }
            if "combo" in self.options:
                ret.update ( self.options["combo"] )
            return ret
        if which == "anaids":
            ret = [ "CMS-EXO-20-004:(comb)", "CMS-SUS-20-004:(comb)" ]
            ret.append ( "ATLAS-EXOT-2018-06:EM10" )
            if "anaids" in self.options:
                ret = self.options["anaids"]
            return ret
        if which == "rmCritic":
            ret = False
            if "rmCritic" in self.options:
                ret = self.options["rmCritic"]
            return ret

        if which in self.options:
            return self.options[which]

        return {}

    def plot ( self ):
        """ this is the method that controls the entire plot.
        adapt according to your purpose!!
        """
        """
        critic_points = self.getCriticList( critic_type = "ul" )
        options = { "label": "excluded by ul" }
        self.plotBooleanMap ( critic_points, options )
        """
        critic_points = self.getCriticList( critic_type = "both" )
        critic_options = self.getOptions ( "critic" )
        self.plotBooleanMap ( critic_points, critic_options )

        colors = [ ( "red", "darkred" ), ( "green", "darkgreen" ),
                   ( "blue", "darkblue" ), ( "purple", "pink" ),
                   ( "brown", "orange" ) ]
        anaids = self.getAnaIds( True )
        builder_options = self.getOptions ( "builder" )
        anaids = self.getOptions("anaids")
        rmCritic = self.getOptions("rmCritic")

        for_combination = {}

        for i,anaid in enumerate ( anaids ):
            nll_points = self.getNLLList( anaid = anaid,
                   removeDisallowed = rmCritic )
            options = copy.deepcopy ( builder_options )
            options.update ( { "label": anaid, "colors": colors[i] } )
            if anaid in self.options:
                options.update ( self.options[anaid] )
            self.plotProbabilityMass ( nll_points, options, critic_points )
            for_combination[anaid] = nll_points
        comb_points = self.combineNLLs ( for_combination )
        comb_options = { "colors": ( "0.20", "black" ), "label": "joint posterior", "nlevels": 1 }
        comb_options = self.getOptions ( "combo" )

        self.plotProbabilityMass ( comb_points, comb_options, critic_points )

        # Existing scatter handles (from plt.scatter calls)
        # handles, labels = plt.gca().get_legend_handles_labels()

        # Add the area patch to the legend
        loc = "best"
        # loc = "upper left"
        plt.legend( handles=self.handles, loc=loc )
        xlabel = rf"m$\left({namer.texName(self.data['meta']['xvariable'])}\right)$ [GeV]"
        ylabel = rf"m$\left({namer.texName(self.data['meta']['yvariable'])}\right)$ [GeV]"
        if type(self.data["meta"]["xvariable"]) == tuple:
            xlabel = rf"ssm$\left({namer.texName(self.data['meta']['xvariable'])}\right)$ [GeV]"
        if type(self.data["meta"]["yvariable"]) == tuple:
            ylabel = rf"ssm$\left({namer.texName(self.data['meta']['yvariable'])}\right)$ [GeV]"
        plt.xlabel( xlabel )
        plt.ylabel( ylabel )
        plt.title ( "posteriors" )
        self.savefig()
        self.pprint ( f"plotting {self.args['inputfile']} -> {self.outputfile}" )

    def plotBooleanMap ( self, points : list[dict], options : dict  ):
        """ given a list of dicts, draw the contour
        :param points: a list of dictionaries, "x", "y", "passes"
        """
        import numpy as np
        defaults = { "c_area": "gray", "c_line": "dimgray",
            "scatter": False, "xlabel": "x", "ylabel": "y", "hatches": "////",
            "label": "excluded", "text": False }
        opts = defaults
        opts.update ( options )

        # Extract arrays
        x = np.array([d["x"] for d in points])
        y = np.array([d["y"] for d in points])
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
                color="red", edgecolor="black", s=10 )

            # passes == True → green
            plt.scatter( x[mask_true], y[mask_true],
                color="green", edgecolor="black", s=10 )

            # passes == True → green
            plt.scatter( x[mask_true], y[mask_true],
                color="green", edgecolor="black", s=10 )

        if opts["text"]:
            robs_min = 2.0
            fontsize=8
            for i,d in enumerate(points):
                #if i % 10 != 0:
                #    continue
                if d["passes"]==True:
                    continue
                if d["robs"] > robs_min:
                    plt.text ( d["x"], d["y"], f"{d['robs']:.1f}",
                               fontsize=fontsize )

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
    argparser.add_argument ( '-O', '--options',
            help="additional options, e.g. { 'critic:scatter': True }", type=str, default="{}" )
    argparser.add_argument ( '--show_options',
            help='start interactive shell', action="store_true" )
    argparser.add_argument ( '-s', '--show',
            help='show image', action="store_true" )
    args = argparser.parse_args()
    plotter = NLLPlotter ( vars(args) )
    if args.show_options:
        print ( f"\noptions" )
        print ( f"=========" )
        ops = plotter.getOptions ( "all" )
        for k,v in ops.items():
            sv = v
            if type(v) in [ str ]:
                sv = f"'{v}'"
            print ( f"'{k}': {sv}" )
        print ( )
        if args.inputfile == None:
            sys.exit()
    plotter.run()
