#!/usr/bin/env python3

""" simple experimental script to plot upper limit ratios.
Lets see if that helps us in identifying regions of interest.
"""

import pickle, sys, subprocess, os
from matplotlib import pyplot as plt
import numpy as np
import IPython
from base.loggerbase import LoggerBase
from ptools.sparticleNames import SParticleNames
from typing import Set

class LimitRatioPlotter ( LoggerBase ):
    def __init__ ( self, picklefile : os.PathLike,
                   outfile : str = "r_@a@.png" ):
        super ( LimitRatioPlotter, self ).__init__ ( 0 )
        self.namer = SParticleNames ( susy = False )
        self.picklefile = picklefile
        f=open ( self.picklefile, "rb" )
        self.dictionary = pickle.load ( f )
        f.close()
        self.outfile = outfile # .replace("@a@",analysis)

    def interact ( self ):
        IPython.embed( colors = "neutral" )

    def createDataGrid( self, analysis : str ):
        """ create the data(mesh?) grid that we wish to plot """
        self.grid = {}
        gridx, gridy = set(), set()
        for masspoint in self.dictionary["masspoints"]:
            mx, my = masspoint["mx"], masspoint["my"]
            oul, eul, r = float("nan"), float("nan"), float("nan")
            for name,oul_ in masspoint["oul"].items():
                if analysis in name:
                    oul = oul_
                    break
            for name,eul_ in masspoint["eul"].items():
                if analysis in name:
                    eul = eul_
                    break
            if type(eul)==float:
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

    def plot ( self, analysis : str ):
        """ plot for specific analysis. """
        fig, ax = plt.subplots()
        self.createDataGrid( analysis )
        gridr = np.array(self.gridr)[:-1, :-1]
        vmin, vmax = 0, 3.
        cmap = 'RdYlGn'
        cmap = 'PiYG'
        c = ax.pcolormesh ( self.gridx, self.gridy, gridr, cmap=cmap, 
                            vmin = vmin, vmax = vmax )
        ## draw a star for the hiscore model
        hi_x, hi_y = self.dictionary["mpid1"], self.dictionary["mpid2"]
        plt.scatter ( [ hi_x ], [ hi_y ], marker="*", s=100, color="black", label = "hiscore" )
        ax.set_title ( f"ratios of limits, {analysis}" )
        plt.xlabel ( f"${self.namer.texName(self.dictionary['pid1'])}$" )
        plt.ylabel ( f"${self.namer.texName(self.dictionary['pid2'])}$" )
        cb = fig.colorbar(c, ax=ax)
        cb.set_label ( "r=oul/eul" )
        outfile = self.getOutputFileName ( analysis )
        self.pprint ( f"saving to {outfile}" )
        plt.legend()
        plt.text ( .8, -.12, "excesses are at r>1", c="gray", 
                   transform = ax.transAxes, fontsize = 10 )
        plt.tight_layout()
        plt.savefig ( outfile )
        plt.clf()

    def findAllAnalyses ( self ) -> Set:
        """ find all analyses names mentioned in picklefile """
        analyses = set()
        for masspoint in self.dictionary["masspoints"]:
            llhds = list(masspoint["llhd"].values())
            for llhd in llhds:
                for name,value in llhd.items():
                    analyses.add ( name[:name.find(":")] )
            critic = list(masspoint["critic"])
            for c in critic:
                name = c.replace("(comb)","").replace("(ul)","").\
                       replace("(em)","")
                analyses.add ( name )
            ouls = list(masspoint["oul"])
            for oul in ouls:
                name = oul[:oul.find(":")]
                analyses.add ( name )
        return analyses

    def getOutputFileName ( self, analysis ):
        return self.outfile.replace("@a@",analysis )

    def show ( self, analysis ):
        cmd = f"see {self.getOutputFileName(analysis)}"
        o = subprocess.getoutput ( cmd )
        print ( o )

if __name__ == "__main__": 
    import argparse
    argparser = argparse.ArgumentParser(
            description='plot upper limit ratio plots (oUL/eUL)')
    argparser.add_argument ( '-a', '--analysis',
            help="analysis to plot. if empty, plot for all.", type=str, default=None )
    argparser.add_argument ( '-p', '--picklefile',
            help="the picklefile", type=str, default="llhdX1WX1Z.pcl" )
    argparser.add_argument ( '-I', '--interactive',
            help='interactive mode',
            action="store_true" )
    argparser.add_argument ( '-s', '--show',
            help='show plot at the end',
            action="store_true" )
    args = argparser.parse_args()
    plotter = LimitRatioPlotter( picklefile = args.picklefile )
    analyses = [ args.analysis ]
    if args.analysis is None:
        # plot for all analyses mentioned
        analyses = plotter.findAllAnalyses()
    for analysis in analyses:
        plotter = LimitRatioPlotter( picklefile = args.picklefile )
        plotter.plot( analysis )
        if args.show:
            plotter.show( analysis )
    if args.interactive:
        plotter.interact()
