#!/usr/bin/env python3
# coding: utf-8

import matplotlib.pyplot as plt
import sys,os,copy,glob
import numpy as np
import seaborn as sns
import pandas as pd
#sys.path.append(os.path.abspath('../smodels'))
#sys.path.append(os.path.abspath('../'))

from smodels.experiment.databaseObj import Database
from smodels.base import runtime

from smodels_utils.plotting.mpkitty import timg

from base.loggerbase import LoggerBase
from builder.protomodel import ProtoModel
from builder.manipulator import Manipulator
from tester.predictor import Predictor
from tester.combiner import Combiner
from base.runEnviron import RunEnviron
from ptools.sparticleNames import SParticleNames
from snippets import show_steps

def ceil_to_step(x : float, step : int = 50) -> float:
    return float ( np.ceil(x / step) * step )

class TenHiscores ( LoggerBase ):
    def __init__ ( self, args ):
        super ( TenHiscores, self ).__init__ ( "ten" )
        self.args = args
        self.folder = "./"

    def standardizeK ( self, K, offset : float = 0. ):
        return ( K - self.Kmin ) / self.Kstd * 2.5 + offset

    def getTicks ( self, ymin : float, ymax : float ) -> tuple:
        ymin = ceil_to_step ( ymin, 100 )
        ymin, ymax = int(ymin), int(ymax+1)
        lim1, lim2 = 600, 1400

        ticks = range ( ymin, ymax, 100 )
        if ymax > lim1:
            ticks = range ( ymin, lim1+1, 100 )
            ticks = list ( ticks ) + list ( range ( lim1+200, ymax, 200 ) )
        if ymax > 1600:
            ticks = range ( ymin, lim1+1, 100 )
            ticks = list ( ticks ) + list ( range ( lim1+200, lim2+1, 200 ) )
            ticks = list ( ticks ) + list ( range ( lim2, ymax, 400 ) )
        return tuple ( ticks )

    def offsetFor ( self, pid : int, idx : int, m : float, 
                    idxNMass : list ) -> float:
        """ compute the offset for particle at xvalue idx,
        yvalue (mass) m """
        offset = - .1
        for pid_,i_,m_ in idxNMass:
            if i_ - idx  == 0. and abs ( m_ - m ) < 8. and pid_ != pid:
                if pid_ > pid:
                    offset = .1
                elif pid_ < pid:
                    offset = -.3
        ret = idx + offset
        return ret

    def plotCombos( self, df, masses ):
        """ create the plot that plots the combinations """
        df = df[1:]
        standardizedKs = list ( self.standardizeK ( df["K"], 0. ) )

        f, axarr = plt.subplots(2,sharex=True,
                gridspec_kw = {'height_ratios':[1, 4]},figsize=(10,8))
        plt.subplots_adjust(left=0.12, bottom=0.12, right=0.97, top=None,
                            wspace=None, hspace=0)

        nsteps = 10

        pids = sorted(list(masses.keys()))

        # axarr[0].scatter(df['walkerid'],df['K'],s=80,c='gray')
        axarr[0].scatter(df['walkerid'],standardizedKs,s=50,c='gray')
        #axarr[0].scatter(df['walkerid'],df['TL'],s=80,c='black')
        axarr[0].set_ylabel(r'$K$')
        axarr[0].set_ylim( self.t_ymin, self.t_ymax )
        axarr[0].set_yticks([])
        for i,row in df.iterrows():
            index = np.where(self.walkerid_values == row['walkerid'])[0][0]
            # print(index)
            if row['K'] == max(df['K']):
                axarr[0].annotate(rf'{row["K"]:1.2f}',(index-1-.2,
                                  self.standardizeK ( row["K"], 1 )),fontsize=10)
                #axarr[0].annotate(r'$\mathbf{%1.2f}$' %row['TL'],(index-0.2,row['TL']+0.5),
                #                  fontsize=10)
            else:
                # axarr[0].annotate(r'$%1.2f$' %row['K'],(index-0.2,row['K']+0.5),fontsize=10)
                axarr[0].annotate(rf'{row["K"]:1.2f}',(index-1-.2, self.standardizeK (row['K'], 1 )),fontsize=10)
                #axarr[0].annotate(r'$%1.2f$' %row['TL'],(index-0.2,row['TL']+0.5),fontsize=10)
            bC = sorted(row['bestCombo'].split(','))
            # print(bC)
            for j, b in enumerate(bC):
                axarr[1].scatter(row['walkerid'], b, c="black", marker='s', s=150)
                index2 = self.analyses.index(b)
                #axarr[1].annotate(f'{b}',(index-0.2,index2),fontsize=10)
        axarr[1].set_xlabel('walkerid:step', fontsize=20)
        axarr[1].set_ylabel('Analyses', fontsize=20)
        axarr[1].set_ylabel('', fontsize=20)
        axarr[1].set_xticks(sorted(df['walkerid'].tolist()), labels =sorted(df['walkerid'].tolist()),  fontsize=12)
        axarr[1].set_yticks(self.analyses, labels =self.analyses,  fontsize=12)
        fname = self.outfile.replace(".png",".combs.png")
        print ( f"saving to {fname}" )
        plt.savefig ( fname )
        show_steps.show(True,self.folder)

    def plot( self ):
        max_entries_per_walk = self.args["nmax_analyses"]
        # max_entries_per_walk = 3

        # things you might want to tweak
        path = f"{self.folder}/all_hiscores/*dict"
        truthfile = f'{self.folder}/truth.dict'
        self.outfile = f"{self.folder}/hiscore.png"
        sns.set_style('ticks')
        sns.set_context('paper', font_scale=2.0)
        # sns.set_palette(sns.color_palette("Paired"))
        sns.set_palette(sns.color_palette("deep"))
        environ=RunEnviron()

        log_file = glob.glob(path)
        print ( os.path.abspath ( log_file[0] ) )

        #Set colors:
        allPids = [ 1000022, 1000006, 1000001, 1000021, 1000012, 1000023,
                    1000013, 2000006, 1000011, 1000005, 1000014, 1000004,
                    1000015, 1000016, 1000024 ]
        namer = SParticleNames ( susy = False )
        colors = sns.color_palette('deep',n_colors=len(namer.xIDs))

        #Replace default colors:
        for pid in sorted(namer.xIDs.keys()):
            if not pid in allPids:
                allPids.append(pid)
        colorDict = dict(zip(allPids,colors))
        colorDict[1000002] = colorDict # masses[pid].append(np.nan)[1000001]
        colorDict[1000003] = colorDict[1000001]
        colorDict[1000004] = colorDict[1000001]

        pTrue = None
        with open ( truthfile, "rt" ) as f:
            pTrue = eval(f.read())

        def fromDict(inputDict):
            p = ProtoModel(walkerid=0,environ=environ)
            for key,v in inputDict.items():
                setattr(p,key,copy.deepcopy(v))

            return p

        #Get highest score from each run:
        protomodelsDict = {}
        #for ff in log_file:
        pList = []
        for fname in log_file:
            with open(fname,'r') as f:
                tList = eval(f.read())
                #run = eval(os.path.basename(ff).replace('real','').replace('.dict',''))
                tmp = [fromDict(pDict) for pDict in tList[:]]
                tmp = sorted(tmp, key = lambda p: p.K, reverse=True)
                pList += tmp[:max_entries_per_walk]
        p = sorted(pList, key = lambda p: p.K, reverse=True)
        p=p[:10]
        #protomodelsDict[run] = p
        p = [ fromDict ( pTrue ) ] + p

        for proto in p:
            walkerid = proto.walkerid
            step = proto.step
            # print(walkerid,proto,'step=',step)

        # Karr = np.array([p["K"] for p in protomodelsDict.values()])
        Karr = np.array([proto.K for proto in p])
        self.Kavg = Karr.mean()
        self.Kstd = Karr.std()
        self.Kmin = min ( Karr )
        print( f'K (avg) = {self.Kavg:1.2f} +- {self.Kstd:1.2f} >= {self.Kmin:1.2f}' )

        #Get all particles which appears in all models:
        particles = []
        self.analyses = []
        #modelList = np.array(list(protomodelsDict.items()))
        #runs = modelList[:,0]
        #modelList = modelList[:,1]
        for proto in p:
            particles += proto.unFrozenParticles()
            ana = proto.description.split(',')
            self.analyses += ana
            #print(analyses)
        particles = list(set(particles))
        self.analyses = sorted(list(set(self.analyses)))

        #Build useful dataset:
        self.walkerid_values = np.array([f"{proto.walkerid}:{proto.step}" for proto in p])
        self.walkerid_values[0]="truth"
        nparticles = np.array([len(proto.unFrozenParticles()) for proto in p])
        Kvalues = np.array([proto.K if (proto.K and proto.K > 0) else 0.0 for proto in p])
        TLvalues = np.array([proto.TL if (proto.TL and proto.TL > 0) else 0.0 for proto in p])
        masses = dict([[pid,[]] for pid in particles])
        for proto in p:
            for pid in masses:
                if pid in proto.masses:
                    masses[pid].append(proto.masses[pid])
                else:
                    masses[pid].append(-100.0)
        #for pid in masses:
        #    masses[pid] = np.array(masses[pid])
        dataDict = {'walkerid': self.walkerid_values, 'K' : Kvalues,
                    'TL': TLvalues, 'nparticles' : nparticles}
        dataDict.update(masses)
        bestCombo_values = np.array([proto.description for proto in p])
        dataDict.update({"bestCombo":bestCombo_values})

        df = pd.DataFrame(dataDict)

        standardizedKs = list ( self.standardizeK ( df["K"], 0. ) )

        self.t_ymin, self.t_ymax = min(standardizedKs)-.8,max(standardizedKs)+2.2

        f, axarr = plt.subplots(2,sharex=True, gridspec_kw = {'height_ratios':[1, 4]},figsize=(10,8))
        plt.subplots_adjust(left=0.12, bottom=0.12, right=0.97, top=None, wspace=None, hspace=0)

        nsteps = 10
        pids = sorted(list(masses.keys()))

        axarr[0].scatter(df['walkerid'],standardizedKs,s=50,c='gray')
        axarr[0].set_ylabel(r'$K$', fontsize=15)
        axarr[0].set_ylim( self.t_ymin, self.t_ymax )
        axarr[0].set_yticks([])

        for i,row in df.iterrows():
            index = np.where(self.walkerid_values == row['walkerid'])[0][0]
            #print(index)
            if row['K'] == max(df['K']):
                axarr[0].annotate(rf'{row["K"]:1.2f}',(index-.2, self.standardizeK ( row["K"], 1)),fontsize=10)
            else:
                axarr[0].annotate(rf'{row["K"]:1.2f}',(index-.2, self.standardizeK (row['K'], 1 )),fontsize=10)

        amasses,idxNMass=[],[]
        for pid in pids:
            for x in df[pid]:
                if x > -99.:
                    amasses.append ( x )
            for i,m in enumerate(masses[pid]):
                    idxNMass.append ( (pid,i,m) )
        for pid in pids:
            label = namer.texName(pid,addDollars=True)
            xvalues = list ( df["walkerid"].keys() )
            sns.scatterplot(x=xvalues,y=df[pid], size=1000,
                    sizes=(1500,1500),marker='_',
                    label=label, legend=False,
                    c=[colorDict[pid]],ax=axarr[1] )
            index = [np.where(self.walkerid_values == d['walkerid'])[0][0] \
                     for j,d in df.iterrows()]
            for i,m in enumerate(masses[pid]):
                if m < 0: continue
                xcoord = self.offsetFor ( pid, index[i], m, idxNMass )
                ycoord = m+2
                axarr[1].annotate( label,
                                   (xcoord,ycoord),fontsize=10)
        ymin, ymax = .8*min(amasses), 1.2*max(amasses)
        if self.args["ymin"] is not None:
            ymin = self.args["ymin"]
        if self.args["ymax"] is not None:
            ymax = self.args["ymax"]
        axarr[1].set_ylim( ymin, ymax )
        axarr[1].set_xlabel('walkerid:step', fontsize=15)
        axarr[1].set_ylabel('mass [GeV]', fontsize=15)
        if self.args["yscale"] not in [ None, "linear" ]:
            tokens = self.args["yscale"].split(":")
            if tokens[0] == "mylog":
                linthresh = 300.
                linscale = 1.
                if len(tokens)>1:
                    linthresh = float ( tokens[1] )
                if len(tokens)>2:
                    linscale = float ( tokens[2] )
                axarr[1].set_yscale('symlog', linthresh=linthresh,
                                    linscale=linscale )
                ticks = self.getTicks ( ymin, ymax )
                axarr[1].set_yticks(ticks)
                axarr[1].set_yticklabels([f"{t}" for t in ticks])
            else:
                axarr[1].set_yscale( *tokens )
        axarr[1].set_xticks(sorted(df['walkerid'].tolist()),
                labels =sorted(df['walkerid'].tolist()),  fontsize=12)
        axarr[1].vlines(x=.5,ymin=ymin,ymax=ymax,linestyle='--',color='gray')
        axarr[0].vlines(x=.5,ymin=self.t_ymin,ymax=self.t_ymax,
                        linestyle='--',color='gray')
        self.pprint ( f"saving to {self.outfile}" )
        plt.savefig ( f"{self.outfile}" )
        timg ( f"{self.outfile}" )
        if False:
            self.plotCombos( df, masses )

if __name__ == "__main__":
    import argparse
    argparser = argparse.ArgumentParser(
            description="the script that plots the ten hiscores")
    argparser.add_argument ( '-N', '--nmax_analyses',
            help='maximum number per analysis [3]', type=int, default=3 )
    argparser.add_argument ( '--ymin',
            help='ymin [auto]', type=float, default=None )
    argparser.add_argument ( '--ymax',
            help='ymax [auto]', type=float, default=None )
    argparser.add_argument ( '--yscale',
            help='yscale argument [linear]', type=str, default=None )
    args=argparser.parse_args()
    plotter = TenHiscores( vars(args) )
    plotter.plot()
