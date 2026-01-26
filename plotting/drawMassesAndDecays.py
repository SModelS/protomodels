#!/usr/bin/env python3

""" code that combines the ruler and the decays plotter """

import re
import numpy as np
import matplotlib.pyplot as plt
import glob, math
#import ast
import plotly.graph_objects as go
from plotly.subplots import make_subplots
import seaborn as sns
sns.set_style('ticks',{'font.family':'Times New Roman',
                  'font.serif':'Times New Roman'})
sns.set_context('paper', font_scale=1.5)
sns.set_palette(sns.color_palette("deep"))

from smodels_utils.plotting.mpkitty import timg
from ptools import sparticleNames
from base.loggerbase import LoggerBase
namer = sparticleNames.SParticleNames ( susy=False )

def ceil_to_step(x, step=50):
    return math.ceil(x / step) * step

class MassesAndDecays ( LoggerBase ):
    def __init__ ( self, model, options ):
        super ( MassesAndDecays, self ).__init__ ( "draw" )
        self.model = model
        self.masses = model["masses"]
        self.decays = model["decays"]
        self.options = self.defaults()
        self.options.update ( options )
        self.getYRange()

    def defaults ( self ):
        ret = { "outfile": "mass_hierarchy.png" }
        ret["colors"]= { 1000022: "black", 1000023: "navy",
                         1000024: "navy", 1000025: "navy",
                         1000006: "brown" }
        ret["scale"]="symlog"
        # ret["scale"]="linear"
        return ret

    def getYRange ( self ):
        ymin, ymax = float("inf"), 0.
        for pid,mass in self.masses.items():
            if mass < ymin:
                ymin = mass
            if mass > ymax:
                ymax = mass
        ymin, ymax = .8 * ymin - 35., 1.05 * ymax + 20.
        if "ymin" in self.options:
            ymin = self.options["ymin"]
        if "ymax" in self.options:
            ymax = self.options["ymax"]
        self.yrange = ( ymin, ymax )

    def dyLabel ( self, pid ):
        """ get the dy of the label, its a bit involved """
        ret = 11
        if pid == 1000022:
            ret = -40
        if self.options["scale"] == "symlog":
            if pid == 1000022:
                ret = -20
            else:
                mass = self.masses[pid]
                ret = 36. * mass / ( self.yrange[1]+self.yrange[0] )

        return ret

    def drawMasses ( self ):
        """ for each particle draw a horizontal line at the respective
        mass value """
        ctParticles = 0
        masses = list ( self.masses.items() )
        self.xpositions = {}
        masses.sort ( key = lambda x: x[1] )
        for pid, mass in masses:
            color = self.options["colors"][pid]
            sgn = 1 if ctParticles % 2 == 0 else -1
            xpos = 45 + sgn * ctParticles * 10
            dx_line = 10 # length of the line
            dx_label = 2 # dx of the label
            if pid == 1000022:
                xpos = 45
                dx_line = 20
                dx_label = 10
            dy_label = self.dyLabel ( pid )
            self.ax.hlines ( mass, xpos, xpos+dx_line, color = color )
            self.xpositions[pid]=xpos
            pname = namer.texName ( pid, addDollars=True )
            self.ax.text( xpos+dx_label, mass+dy_label, pname, fontsize=20,
                          color = color )
            ctParticles+=1

    def drawDecays ( self ):
        for mpid,decay in self.decays.items():
            for i,(dpids,br) in enumerate(decay.items()):
                self.drawDecay ( mpid, dpids, br, i, len(decay) )

    def drawDecay ( self, mpid : int, dpids : tuple, br : float, 
                    i : int, n : int ):
        """ draw a single grey arrow connecting mother pid
        with daughter pid """
        dpid = dpids[0]
        d_products = dpids[1:]
        opposite_signs = [ (2,1), (11,11), (12,12) ]
        for os_pair in opposite_signs:
            if d_products == os_pair:
                d_products = ( d_products[0], -d_products[1] )
        if d_products == (11,-11):
            br = 3*br
        label = namer.texName ( d_products, addDollars=True, lightFlavors=False,
                                addSign = True, separator = " ")
        # print ( f"@@0 {d_products} -> {label}" )
        if br < 1.0:
            label = f"{label}:{int(round(100*br)):d}%"
        y_start = self.masses[mpid]
        y_end = self.masses[dpid]
        x_start = self.xpositions[mpid]
        x_end = self.xpositions[dpid]
        dx = 5
        if dpid == 1000022:
            dx = 10
        arrowprops = dict(color='gray', arrowstyle='<-')
        plt.annotate( '', ( x_start+5, y_start ),
                      xytext=(x_end+dx,y_end),arrowprops=arrowprops,
                      ha='center')
        dx = x_end - x_start
        dy = y_end - y_start
        x_coord =  x_start + .5 * dx
        di = i * 25
        if self.options["scale"]=="symlog":
            di = i * 13 * y_start / 100.
        di = di - ( n -1 ) * 10 * y_start / 100.
        y_coord =  y_start + .5 * dy - 7. - di
        if x_coord > 50:
            x_coord += 8
        plt.text( x_coord, y_coord, label, fontsize=15 )
        #self.pprint ( mpid, dpids, br, label )

    def init ( self ):
        fig, ax = plt.subplots()
        ax.get_xaxis().set_visible(False)
        ax.spines.top.set_visible(False)
        ax.set_ylim( self.yrange )
        self.fig, self.ax = fig, ax

    def interact ( self ):
        import sys, IPython; IPython.embed( colors = "neutral" )

    def plot ( self ):
        self.init()
        plt.ylabel('Mass [GeV]')
        self.drawMasses()
        self.drawDecays()

        if self.options["scale"]=="symlog":
            self.yrange = ( self.yrange[0], self.yrange[1]*1.2 )
            self.ax.set_yscale('symlog', linthresh=200, linscale=1. )
            ymin = ceil_to_step ( self.yrange[0], 100 )
            ymax = ceil_to_step ( self.yrange[1], 100 )
            ticks = range ( ymin, ymax+1, 100 )
            self.ax.yaxis.set_ticks(ticks)
            self.ax.yaxis.set_ticklabels([f"{t}" for t in ticks])
            self.ax.set_xlim(0,100)

        self.ax.yaxis.grid(True)

        plt.tight_layout()
        outfile = self.options["outfile"]
        plt.savefig ( outfile )
        self.pprint ( f"saving to {outfile}" )
        timg ( outfile )
        if self.options["interact"]:
            self.interact()


def getModel( modelfile : str = "truth.dict" ) -> dict:
    with open ( modelfile, "rt" ) as f:
        ret = eval(f.read())
        if type(ret) == list and type(ret[0]) == dict:
            return ret[0]
        return ret

if __name__ == "__main__":
    import argparse
    argparser = argparse.ArgumentParser(
            description="draw masses and decays plots")
    argparser.add_argument ( '-s', '--scale',
            help='scale: linear, or symlog [symlog]',
            type=str, default='symlog' )
    argparser.add_argument ( '-m', '--modelfile',
            help='path to model file [truth.dict]',
            type=str, default='truth.dict' )
    argparser.add_argument ( '-i', '--interact',
            help='enter interactive mode', action="store_true" )
    args=argparser.parse_args()
    model = getModel( args.modelfile )
    options = vars ( args )
    plotter = MassesAndDecays( model, options )
    plotter.plot()
