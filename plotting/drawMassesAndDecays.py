#!/usr/bin/env python3

""" code that combines the ruler and the decays plotter """

import re
import numpy as np
import matplotlib.pyplot as plt
import glob
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
        return ret

    def getYRange ( self ):
        ymin, ymax = float("inf"), 0.
        for pid,mass in self.masses.items():
            if mass < ymin:
                ymin = mass
            if mass > ymax:
                ymax = mass
        ymin, ymax = .8 * ymin - 30., 1.2 * ymax
        if "ymin" in self.options:
            ymin = self.options["ymin"]
        if "ymax" in self.options:
            ymax = self.options["ymax"]
        self.yrange = ( ymin, ymax )

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
            dy_label = 5 # dy of the label
            if pid == 1000022:
                xpos = 45
                dx_line = 20
                dx_label = 10
                dy_label = -40
            self.ax.hlines ( mass, xpos, xpos+dx_line, color = color )
            self.xpositions[pid]=xpos
            pname = namer.texName ( pid, addDollars=True )
            self.ax.text( xpos+dx_label, mass+dy_label, pname, fontsize=20)
            ctParticles+=1

    def drawDecays ( self ):
        for mpid,decay in self.decays.items():
            for dpids,br in decay.items():
                self.drawDecay ( mpid, dpids, br )

    def drawDecay ( self, mpid, dpids, br ):
        """ draw a single grey arrow connecting mother pid
        with daughter pid """
        dpid = dpids[0]
        d_products = dpids[1:]
        if d_products == (2,1):
            d_products = (2,-1)
        label = namer.texName ( d_products, addDollars=True, lightFlavors=False,
                                addSign = True )
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
        plt.text( x_start + .5 *dx , y_start + .5 * dy , 
                 label, fontsize=15)
        self.pprint ( mpid, dpids, br, label )

    def init ( self ):
        fig, ax = plt.subplots()
        ax.get_xaxis().set_visible(False)
        ax.spines.top.set_visible(False)
        ax.set_ylim( self.yrange )
        ax.set_xlim(0,100)
        ax.yaxis.grid()
        ax.set
        self.fig, self.ax = fig, ax

    def plot ( self ):
        self.init()
        plt.ylabel('Mass [GeV]')
        self.drawMasses()
        self.drawDecays()
        plt.tight_layout()
        plt.savefig ( self.options["outfile"] )
        timg ( self.options["outfile"] )

    def plotSN ( self ):
        model = self.model
        fig, (ax1, ax2) = plt.subplots(
            2, 1, sharex=True,
            gridspec_kw={'height_ratios': [1, 3]}
        )

        # Upper axis (true 450 GeV region)
        ax1.set_ylim(430, 570)
        ax1.hlines(model['masses'][1000025], 85, 90, color='navy')
        ax1.hlines(model['masses'][1000006], 93, 99, color='brown')

        ax1.annotate('', (87, 430), xytext=(87,model['masses'][1000025]),arrowprops=dict(color='gray', arrowstyle='-'), ha='center')
        ax1.text(86, model['masses'][1000025]+5, r'$X_{Z}^{3}$', fontsize=20)
        ax1.text(95, model['masses'][1000006]+5, r'$X_{t}^{1}$', fontsize=20)

        ax2.set_ylim(60, 180)
        ax2.hlines(model['masses'][1000022], 90,95, color='dodgerblue')
        ax2.hlines(model['masses'][1000023], 85,90, color='blue')
        ax2.hlines(model['masses'][1000024], 100,105, color='limegreen')


        ax1.spines.bottom.set_visible(False)
        ax1.get_xaxis().set_visible(False)
        ax2.get_xaxis().set_visible(False)
        ax2.spines.top.set_visible(False)

        plt.annotate('', (92, model['masses'][1000022]), xytext=(87,model['masses'][1000023]),arrowprops=dict(color='gray', arrowstyle='->'), ha='center')
        plt.annotate('', (92.5, model['masses'][1000022]), xytext=(103,model['masses'][1000024]),arrowprops=dict(color='gray', arrowstyle='->'), ha='center')
        plt.annotate('', (87, model['masses'][1000023]+15), xytext=(87,180),arrowprops=dict(color='gray', arrowstyle='->'), ha='center')
        plt.annotate('', (92, model['masses'][1000022]), xytext=(94,180),arrowprops=dict(color='gray', arrowstyle='->'), ha='center')


        plt.text(86, model['masses'][1000022]+10, r"$\nu_l,\bar{\nu}_l$", fontsize=20)
        plt.text(99, model['masses'][1000022]+30, r"$q,\bar{q}$", fontsize=20)
        plt.text(88, model['masses'][1000023]+50, r"$h$", fontsize=20)

        plt.text(86, model['masses'][1000023]+5, r'$X_{Z}^{2}$', fontsize=20)
        plt.text(101, model['masses'][1000024]+5, r'$X_{W}^{1}$', fontsize=20)
        plt.text(91, model['masses'][1000022]-13, r'$X_{Z}^{1}$', fontsize=20)

        ax1.yaxis.grid()
        ax2.yaxis.grid()
        ax2.set_xlabel('')
        plt.ylabel('Mass [GeV]')
        #ax.set_ylim(60,200)
        ax1.set_xlim(81,110)
        plt.tight_layout()
        outfile = "mass_hierarchy.png"
        plt.savefig ( outfile )
        timg ( outfile )

def getModel():
    with open ( "truth.dict", "rt" ) as f:
        return eval(f.read())

if __name__ == "__main__":
    model = getModel()
    options = {}
    plotter = MassesAndDecays( model, options )
    plotter.plot()
