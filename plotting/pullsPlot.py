#!/usr/bin/env python3

""" small class to create pulls plots of analyses """

import numpy as np

from smodels.experiment.expResultObj import ExpResult
from typing import Union
from base.loggerbase import LoggerBase

class PullsPlotter ( LoggerBase ):
    def __init__ ( self ):
        super ( PullsPlotter, self ).__init__ ( "pulls" )
        self.default_options = { "sort_by_abs": False, "annotate": True,
            "pngName": "pulls@@ANALYSIS@@.png", "show_average": True,
            "sort_lexigraphically": True,
            "title": "Pulls by Signal Region, @@ANALYSIS@@" }
        self.explanations_options = { 
            "sort_by_abs": "shall we sort by absolute values", 
            "annotate": "do we annotate",
            "pngName": "filename of output file", 
            "show_average": "add a point with the average value",
            "sort_lexigraphically": "shall we sort by natsort",
            "title": "plot title" 
        }

    def createDictionary ( self, expResults : Union[list,ExpResult] ) -> dict:
        """ given a list of experimental results (or a single one), 
        create the dictionary 
        :param expResults: either a single experimental result,
        or a list
        :returns: e.g. {'SR_2016_1': 0.52, 'SR_2016_2': ... }
        if a list of expResults is given, we return e.g.:
        { 'CMS-EXO-20-004:SR_2016_1': ..., 'ATLAS-SUSY-2018.06:...' } 
        """
        from ptools.helpers import computePForDataSet, computeZFromP
        d = {}
        if type(expResults) in [ list, tuple ]:
            anaIds = []
            for eR in expResults:
                anaId = eR.globalInfo.id
                anaIds.append ( anaId )
                tmp = createDictionary ( eR )
                for k,v in tmp.items():
                    d[ f"{anaId}:{k}" ] = v
            self.anaId = ",",join ( anaIds )
            return d

        self.anaId = expResults.globalInfo.id
        for ds in expResults.datasets:
            p = computePForDataSet ( ds )
            Z = computeZFromP ( p )
            d [ ds.dataInfo.dataId ] = Z

        return d

    def getTitle ( self ):
        ret = self.options["title"]
        ret = ret.replace( "@@ANALYSIS@@", self.anaId )
        return ret

    def getPngName ( self ):
        ret = self.options["pngName"]
        ret = ret.replace( "@@ANALYSIS@@", self.anaId.replace("-","_") )
        return ret


    def plot ( self, values : dict, options : dict ):
        """ plot the pulls from the values dictionary 
        :param values: the values as obtained from self.createDictionary
        :param options: various plotting options
        """
        # Example input
        # values = {"SR1": 1.2, "SR2": -0.8, "SR3": 2.5, "SR4": 0.0, "SR5": -1.7}
        import copy
        from matplotlib import pyplot as plt
        from smodels_utils.helper.various import pngMetaInfo
        opts = copy.deepcopy ( self.default_options )
        opts.update ( options )
        self.options = opts

        if not values:
            raise ValueError("Input dictionary 'values' is empty.")

        # Optionally sort by |Z| (largest first)
        items = list(values.items())
        if opts["sort_by_abs"]:
            items = sorted( values.items(), key=lambda kv: abs(kv[1]), 
                            reverse=True)
        if opts["sort_lexigraphically"]:
            from natsort import natsorted
            items = natsorted( values.items(), key=lambda kv: kv[0] )

        labels, pulls = zip(*items)
        pulls = np.array(pulls, dtype=float)
        y = np.arange(len(labels))

        # Colors by sign
        colors = [ "tab:red" if z > 0 else ("tab:blue" if z < 0 else "gray") \
                   for z in pulls ]

        fig, ax = plt.subplots(figsize=(7, 0.5 * len(labels) + 1.5))

        # Thick vertical line at Z=0
        ax.axvline(0.0, color="k", lw=3, alpha=0.8, zorder=1)

        # Scatter points for pulls
        ax.scatter(pulls, y, s=60, c=colors, edgecolor="k", zorder=3)

        # Y-axis: region names
        ax.set_yticks(y)
        ax.set_yticklabels(labels)
        ax.invert_yaxis()  # largest (top of list) at top

        # X-axis: symmetric limits around 0 with some padding
        max_abs = max(1.0, np.max(np.abs(pulls)))
        pad = 0.2 * max_abs
        ax.set_xlim(-max_abs - pad, max_abs + pad)
     
        ax.set_xlabel("Pull Z")
        ax.set_title( self.getTitle() )



        if opts["show_average"]:
            # Plot average point (gold star)
            z_mean = float(np.mean(pulls))
            y_avg = y[-1]+1
            ax.scatter( [z_mean], [y_avg], marker="*", s=220, color="gold", 
                      edgecolor="k", zorder=4)

        # Light grid on x for readability
        ax.grid(axis="x", linestyle="--", alpha=0.3)

        # Optional annotations with numeric Z next to each point
        if opts["annotate"]:
            for x, yi in zip(pulls, y):
                # Shift text slightly away from the point
                shift = 0.03 * (max_abs + pad)
                ha = "left" if x >= 0 else "right"
                x_text = x + (shift if x >= 0 else -shift) if x != 0 else x + shift
                ax.text(x_text, yi, f"{x:.2f}", va="center", ha=ha, fontsize=9)
            if opts["show_average"]:
                x = z_mean
                ha = "left" if x >= 0 else "right"
                x_text = x + (shift if x >= 0 else -shift) if x != 0 else x + shift
                ax.text(x_text, y_avg, f"{x:.2f} (avg)", va="center", ha=ha, fontsize=9, fontweight="bold")

        fig.tight_layout()

        metadata = pngMetaInfo()
        pngName= self.getPngName()
        self.pprint ( f"saving to {pngName}" )
        plt.savefig ( pngName, metadata=metadata )

    def show ( self ):
        from smodels_utils.plotting.mpkitty import timg
        timg ( self.getPngName() )
        # return fig, ax

    def interact ( self ):
        import sys, IPython; IPython.embed( colors = "neutral" )

    def filterBySRName ( self, d : dict, srname : str,
                         invert_selection : bool = False ) -> dict:
        """ filter by the signal region names, wildcards allowed 
        :param srname: e.g. "SR_2018*"
        :param invert_selection: if true, then invert the selection

        :returns: only entries that pass
        """
        ret = {}
        import fnmatch
        for sr, value in d.items():
            if invert_selection and not fnmatch.fnmatch ( sr, srname ):
                ret[sr]=value
            if not invert_selection and fnmatch.fnmatch ( sr, srname ):
                ret[sr]=value
        return ret

    def printOptions ( self) :
        print ( f"plotting options" )
        print ( f"================" )
        for k,v in self.default_options.items():
            expl = self.explanations_options[k]
            print ( f"{k}({v}): {expl}" )

if __name__ == "__main__":
    import argparse
    plotter = PullsPlotter()
    argparser = argparse.ArgumentParser(
            description="draw pulls plot")
    argparser.add_argument ( '-a', '--analysisId',
            help='analysisId [ATLAS-EXOT-2018-06]',
            type=str, default='ATLAS-EXOT-2018-06' )
    argparser.add_argument ( '-s', '--show_options',
            action = "store_true" )
    for k,v in plotter.default_options.items():
        if type(v) == bool:
            argparser.add_argument ( f'--{k}', 
                    help=plotter.explanations_options[k],
                    action="store_true" )
        else:
            argparser.add_argument ( f'--{k}', 
                    help=plotter.explanations_options[k], default=v )
    args=argparser.parse_args()
    if args.show_options:
        plotter.printOptions()
        import sys; sys.exit()

    from base.runEnviron import RunEnviron
    environ = RunEnviron( "run.dict" )
    db = environ.database
    anaIds = [ args.analysisId ]
    eRs = db.getExpResults( analysisIDs = anaIds, dataTypes=["efficiencyMap"] )
    er = eRs[0]
    d = plotter.createDictionary( er )
    # d = plotter.filterBySRName ( d, "SR_2018*" )
    opts = {}
    plotter.plot ( d, opts )
    if True:
        plotter.show()
    if False:
        plotter.interact()
