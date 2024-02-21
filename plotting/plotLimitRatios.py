#!/usr/bin/env python3

""" simple experimental script to plot upper limit ratios.
Lets see if that helps us in identifying regions of interest.
"""

import pickle, sys, subprocess, os
from matplotlib import pyplot as plt
import numpy as np
import IPython
from builder.loggerbase import LoggerBase
from ptools.sparticleNames import SParticleNames

class LimitRatioPlotter ( LoggerBase ):
    def __init__ ( self, picklefile : os.PathLike, analysis : str,
                   outfile : str = "r_@a@.png" ):
        super ( LimitRatioPlotter, self ).__init__ ( 0 )
        self.namer = SParticleNames ( susy = False )
        self.picklefile = picklefile
        self.analysis = analysis
        f=open ( self.picklefile, "rb" )
        self.dictionary = pickle.load ( f )
        f.close()
        self.outfile = outfile.replace("@a@",analysis)

    def interact ( self ):
        IPython.embed( colors = "neutral" )

    def createDataGrid( self ):
        """ create the data(mesh?) grid that we wish to plot """
        self.grid = {}
        gridx, gridy = set(), set()
        for masspoint in self.dictionary["masspoints"]:
            mx, my = masspoint["mx"], masspoint["my"]
            oul, eul, r = float("nan"), float("nan"), float("nan")
            for name,oul_ in masspoint["oul"].items():
                if self.analysis in name:
                    oul = oul_
                    break
            for name,eul_ in masspoint["eul"].items():
                if self.analysis in name:
                    eul = eul_
                    break
            if eul>0:
                r = oul/eul
            self.debug ( f"mx={mx:.2f} my={my:.2f} oul={oul:.2f} eul={eul:.2f} r={r:.2f}" )
            self.grid[ (int(mx),int(my)) ] = r
            gridx.add ( mx )
            gridy.add ( my )
        gridx=list(gridx); gridx.sort()
        gridy=list(gridy); gridy.sort()
        if False:
            print ( "grid", self.grid )
        self.gridx, self.gridy = np.meshgrid ( gridx, gridy )
        self.gridr = []
        for y in gridy:
            tmp = []
            for x in gridx:
                r = float("nan")
                if (int(x),int(y)) in self.grid:
                    r = self.grid[(int(x),int(y))]
                tmp.append ( r )
            self.gridr.append ( tmp )

    def plot ( self ):
        fig, ax = plt.subplots()
        self.createDataGrid()
        gridr = np.array(self.gridr)[:-1, :-1]
        vmin, vmax = 0, 3.
        cmap = 'RdYlGn'
        cmap = 'PiYG'
        c = ax.pcolormesh ( self.gridx, self.gridy, gridr, cmap=cmap, 
                            vmin = vmin, vmax = vmax )
        ## draw a star for the hiscore model
        hi_x, hi_y = self.dictionary["mpid1"], self.dictionary["mpid2"]
        plt.scatter ( [ hi_x ], [ hi_y ], marker="*", s=100, color="black", label = "hiscore" )
        ax.set_title ( f"ratios of limits, {self.analysis}" )
        plt.xlabel ( f"${self.namer.texName(self.dictionary['pid1'])}$" )
        plt.ylabel ( f"${self.namer.texName(self.dictionary['pid2'])}$" )
        cb = fig.colorbar(c, ax=ax)
        cb.set_label ( "r=oul/eul" )
        self.pprint ( f"saving to {self.outfile}" )
        plt.legend()
        plt.text ( .8, -.1, "excesses are at r>1", c="gray", 
                   transform = ax.transAxes, fontsize = 10 )
        plt.tight_layout()
        plt.savefig ( self.outfile )

    def show ( self ):
        cmd = f"see {self.outfile}"
        o = subprocess.getoutput ( cmd )
        print ( o )

if __name__ == "__main__": 
    import argparse
    argparser = argparse.ArgumentParser(
            description='plot upper limit ratio plots (oUL/eUL)')
    argparser.add_argument ( '-a', '--analysis',
            help="analysis", type=str, default="ATLAS-SUSY-2019-09" )
    argparser.add_argument ( '-p', '--picklefile',
            help="the picklefile", type=str, default="llhd10000241000022.pcl" )
    argparser.add_argument ( '-I', '--interactive',
            help='interactive mode',
            action="store_true" )
    argparser.add_argument ( '-s', '--show',
            help='show plot at the end',
            action="store_true" )
    args = argparser.parse_args()
    # analysis = 'CMS-SUS-16-048'
    plotter = LimitRatioPlotter( picklefile = args.picklefile, 
            analysis = args.analysis )
    plotter.plot()
    if args.show:
        plotter.show()
    if args.interactive:
        plotter.interact()
