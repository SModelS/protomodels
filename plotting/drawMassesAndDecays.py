#!/usr/bin/env python3

""" code that combines the ruler and the decays plotter """

import math
from typing import Union

from smodels_utils.plotting.mpkitty import timg
from ptools import sparticleNames
from base.loggerbase import LoggerBase
namer = sparticleNames.SParticleNames ( susy=False )

def ceil_to_step(x : float, step : int = 50) -> float:
    return math.ceil(x / step) * step

def floor_to_step(x : float, step : int = 50) -> float:
    return math.floor(x / step) * step

class MassesAndDecays ( LoggerBase ):
    def __init__ ( self, model, options ):
        super ( MassesAndDecays, self ).__init__ ( "draw" )
        self.model = model
        self.masses = model["masses"]
        self.decays = model["decays"]
        self.options = self.defaults()
        self.options.update ( options )
        self.importMatplot()
        self.getYRange()

    def importMatplot ( self ):
        from smodels_utils.plotting.plottingRecorder import importMatplot
        self.plt = importMatplot ( self.options["record"], False )

    def defaults ( self ):
        ret = { "outfile": "mass_hierarchy.png" }
        ret["colors"]= { 1000022: "black", 1000023: "navy",
                         1000024: "navy", 1000025: "navy",
                         1000006: "brown" }
        ret["scale"]="symlog"
        ret["record"]=False
        # ret["scale"]="linear"
        return ret

    def getYRange ( self ):
        ymin, ymax = float("inf"), 0.
        for pid,mass in self.masses.items():
            if mass < ymin:
                ymin = mass
            if mass > ymax:
                ymax = mass
        dy = ymax - ymin
        # ymin, ymax = .8 * ymin - 35., 1.05 * ymax + 20.
        # ymin, ymax = .8 * ymin - dy / 20., 1.05 * ymax + dy / 30.
        ymin, ymax = .95 * ymin - dy / 10., 1.05 * ymax + dy / 30.
        if "ymin" in self.options and self.options["ymin"] is not None:
            ymin = self.options["ymin"]
        if "ymax" in self.options and self.options["ymax"] is not None:
            ymax = self.options["ymax"]
        self.pprint ( f"setting y range to ({ymin:.1f},{ymax:.1f})" )
        self.yrange = ( ymin, ymax )
        self.delta_y = ymax - ymin # convenience

    def dyLabel ( self, pid ):
        """ get the dy of the label, its a bit involved """
        ret = 11
        if pid == 1000022:
            ret = - self.delta_y / 10.
        if self.options["scale"] == "symlog":
            if pid == 1000022:
                ret = - self.delta_y / 20.
            else:
                mass = self.masses[pid]
                ret = self.delta_y / 10. * mass / ( self.yrange[1]+self.yrange[0] )
        if self.options["scale"] == "linear":
            if pid == 1000022:
                ret = -.1 * self.delta_y
            else:
                ret = .05 * self.delta_y

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
            if mass + dy_label < self.yrange[0] - self.delta_y / 60.:
                dy_label = self.yrange[0] - 0. - mass
            self.ax.text( xpos+dx_label, mass+dy_label, pname, fontsize=20,
                          color = color )
            ctParticles+=1

    def drawDecays ( self ):
        toDraw = {}
        ## first we sum up over similar labels
        for mpid,decay in self.decays.items():
            if not mpid in toDraw:
                toDraw[mpid]={}
            for dpids,br in decay.items():
                bsm_dpid = dpids[0]
                label = self.dpidsToStr ( dpids )
                if not bsm_dpid in toDraw[mpid]:
                    toDraw[mpid][bsm_dpid]={}
                if not label in toDraw[mpid][bsm_dpid]:
                    toDraw[mpid][bsm_dpid][label]=0.
                toDraw[mpid][bsm_dpid][label]+=br
        # only then do we draw
        for mpid,decay in toDraw.items():
            for bsm_dpid,radiates in decay.items():
                for i,(label,br) in enumerate(radiates.items()):
                    self.drawDecay ( mpid, bsm_dpid, label, br, i, len(decay) )

    def dpidsToStr ( self, dpids : tuple,
            br : Union[None,float] = None ) -> str:
        """ translate the decay pids to a string
        e.g. (1000022, 11, 12) -> l+ nu
        :param br: optional branching ratio
        """
        dpid = dpids[0]
        d_products = dpids[1:]
        opposite_signs = [ (2,1), (11,11), (12,12), (13,13), (15,15) ]
        for os_pair in opposite_signs:
            if d_products == os_pair:
                d_products = ( d_products[0], -d_products[1] )
        if d_products == (11,-11) and br is not None:
            br = 3*br
        label = namer.texName ( d_products, addDollars=True, lightFlavors=False,
                                addSign = "ifboth", separator = " ")
        if br is not None and br < 1.0:
            label = rf"{label}:{int(round(100*br)):d}%"
        return label

    def drawDecay ( self, mpid : int, dpid : int, label : str, br : float,
                    i : int, n : int ):
        """ draw a single grey arrow connecting mother pid
        with daughter pid """
        if br is not None and br < 0.999:
            label = f"{label}:{int(round(100*br)):d}%"
        y_start = self.masses[mpid]
        y_end = self.masses[dpid]
        x_start = self.xpositions[mpid]
        x_end = self.xpositions[dpid]
        dx = 5
        if dpid == 1000022:
            dx = 10
        arrowprops = dict(color='gray', arrowstyle='<-')
        self.plt.annotate( '', ( x_start+5, y_start ),
                      xytext=(x_end+dx,y_end),arrowprops=arrowprops,
                      ha='center')
        dx = x_end - x_start
        dy = y_end - y_start
        x_coord =  x_start + .5 * dx
        di = i * self.delta_y / 20.
        dm = 100.
        if self.options["scale"]=="symlog":
            di = i * self.delta_y / 40. * y_start / dm
        di = di - ( n -1 ) * self.delta_y / 50. * y_start / dm
        y_coord =  y_start + .5 * dy - self.delta_y / 60. - di
        if x_coord > 50:
            x_coord += 8
        if y_coord < self.yrange[0]:
            y_coord = self.yrange[0]
        self.plt.text( x_coord, y_coord, label, fontsize=15 )

    def init ( self ):
        fig, ax = self.plt.subplots()
        if self.options["record"]:
            self.fig = self.plt.intercept ( fig, "fig" )
            self.ax = self.plt.intercept ( ax, "ax" )
        else:
            self.fig, self.ax = fig, ax
        #self.plt.subplots()
        #self.fig = self.plt.gcf()
        #self.ax = self.plt.gca()
        self.ax.get_xaxis().set_visible(False)
        self.ax.spines.top.set_visible(False)
        self.ax.set_ylim( self.yrange )

    def interact ( self ):
        import sys, IPython; IPython.embed( colors = "neutral" )

    def plot ( self ):
        self.init()
        self.plt.ylabel('Mass [GeV]')
        self.drawMasses()
        self.drawDecays()
        dm = 100
        if self.delta_y < 100.:
            dm = 50
        if self.delta_y < 50.:
            dm = 20
        if self.delta_y < 20.:
            dm = 10

        if self.options["scale"]=="symlog":
            # self.yrange = ( self.yrange[0], self.yrange[1]*1.2 )
            self.ax.set_yscale('symlog', linthresh=200, linscale=1. )
            ymin = ceil_to_step ( self.yrange[0], dm )
            ymax = ceil_to_step ( self.yrange[1], dm )
            ticks = range ( ymin, ymax+1, dm )
            self.ax.set_yticks(ticks)
            self.ax.set_yticklabels([f"{t}" for t in ticks])
            self.ax.set_xlim(0,100)
        # print ( f"@@ ymax {ymax} dm {dm} yrange {self.yrange}" )

        self.ax.get_yaxis().grid(True)

        self.plt.tight_layout()
        outfile = self.options["outfile"]
        from smodels_utils.helper.various import pngMetaInfo
        metadata = pngMetaInfo()
        self.plt.savefig ( outfile, metadata = metadata )
        self.pprint ( f"saving to {outfile}" )
        timg ( outfile )
        if self.options["interact"]:
            self.interact()


def getModel( modelfile : str = "truth.dict", 
              model_index : int = 0 ) -> dict:
    with open ( modelfile, "rt" ) as f:
        ret = eval(f.read())
        if type(ret) == list and type(ret[0]) == dict:
            if model_index >= len(ret):
                print ( f"[drawMassesAndDecays] error, list does not have {model_index+1} entries" )
                import sys; sys.exit()
            return ret[model_index]
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
    argparser.add_argument ( '--model_index',
            help='which model within the model file, in case its a list [0]',
            type=int, default=0 )
    argparser.add_argument ( '-o', '--outfile',
            help='output file [mass_hierarchy.png]',
            type=str, default='mass_hierarchy.png' )
    argparser.add_argument ( '-i', '--interact',
            help='enter interactive mode', action="store_true" )
    argparser.add_argument ( '-r', '--record',
            help='active plotting recorder', action="store_true" )
    argparser.add_argument ( '--ymin',
            help='ymin [auto]', type=float, default=None )
    argparser.add_argument ( '--ymax',
            help='ymax [auto]', type=float, default=None )
    args=argparser.parse_args()
    model = getModel( args.modelfile, args.model_index )
    options = vars ( args )
    plotter = MassesAndDecays( model, options )
    plotter.plot()
