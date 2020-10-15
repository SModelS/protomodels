#!/usr/bin/env python3
import warnings
warnings.filterwarnings("ignore")

import sys,os,copy, pickle, glob, subprocess, math
import numpy as np
import IPython
sys.path.append(os.path.abspath('../../smodels'))
sys.path.append(os.path.abspath('../'))
from builder.protomodel import ProtoModel
from tester.predictor import Predictor
from tester.combiner import Combiner
from smodels.experiment.databaseObj import Database
from smodels.tools import runtime
runtime._experimental = True
import matplotlib.pyplot as plt
import seaborn as sns
# from names import particleLabels
from ptools.sparticleNames import SParticleNames
import pandas as pd
# sns.set() #Set style
# sns.set_style('ticks')


sns.set_style('ticks',{'font.family':'Times New Roman', 'font.serif':'Times New Roman'})
sns.set_context('paper', font_scale=1.8)
# sns.set_palette(sns.color_palette("Paired"))
sns.set_palette(sns.color_palette("deep"))

cmd = 'rm -f walk*.png'
subprocess.getoutput( cmd )

#Set colors:
namer = SParticleNames ( susy = False )
#Replace default colors:
colorPalette = 'deep' #Set color palette for particles, if None use default
allpids = [ 1000001, 1000002, 1000003, 1000004, 1000005, 1000006,          
2000005, 2000006, 1000011, 1000012, 1000013, 1000014, 1000015,          
1000016, 1000021, 1000022, 1000023, 1000025, 1000024, 1000037 ]
colorDict = dict(zip( allpids,sns.color_palette(palette=colorPalette,n_colors=len(namer.names))))

import argparse
argparser = argparse.ArgumentParser( description="movie maker" )
argparser.add_argument ( '-f', '--history',
        help='history file to use [history.list]',
        type=str, default="history.list" )
argparser.add_argument ( '-m', '--maxsteps',
        help='maximum steps [1000]',
        type=int, default=1000 )
args = argparser.parse_args()

f=open(args.history,"rt")
txt=f.read()
f.close()
txt=txt.replace("nan","'nan'")
if not "]" in txt[-3:]:
    txt+="]\n"

modelList=eval(txt)

#Get all particles which appears in all steps:
particles = []
for p in modelList:
    particles += p["masses"].keys()
particles = list(set(particles))

#Build useful dataset:
steps = np.array([p["step"] for p in modelList])
nparticles = np.array([len(p["masses"]) for p in modelList])
Kvalues = np.array([p["K"] if (p["K"] and p["K"] > 0) else 0.0 for p in modelList])
Zvalues = np.array([p["Z"] if (p["Z"] and p["Z"] > 0) else 0.0 for p in modelList])
masses = dict([[pid,[]] for pid in particles])
Ks,actions,bcs=[],[],[]
for p in modelList:
    K=p["K"]
    if K == None:
        K=0.
    Ks.append(K)
    ac = []
    if "actions" in p:
        ac = p["actions"]
    actions.append(ac)
    bcs.append(p["bestCombo"])
    for pid in masses:
        if pid in p["masses"]:
            masses[pid].append(float(p["masses"][pid]))
        else:
#             masses[pid].append(np.nan)
            masses[pid].append(-100.0)
for pid in masses:
    masses[pid] = np.array(masses[pid])
dataDict = {'step' : steps, 'K' : Kvalues, 'Z' : Zvalues, 
                   'nparticles' : nparticles}
dataDict.update(masses) 
df = pd.DataFrame(dataDict)
fig = plt.figure(figsize=(10, 6))
nsteps = 1
maxstep = args.maxsteps
#maxstep = 200

maxK,stepatmax=0.,0

if maxstep > len(Ks):
    maxstep=len(Ks)

for firststep in range ( maxstep ):
    fig, (ax1, ax2) = plt.subplots( ncols=2, sharey=True, gridspec_kw={'width_ratios': [1, 10]} ) 
    if firststep % 10 == 0:
        print ( "step %d" % firststep )
    laststep=firststep+20

    nvalues = {}
    for k,v in masses.items():
        nvalues[k]=0
        for i in v:
            if float(i)>0:
                    nvalues[k]+=1
    pids = sorted(masses.keys(), key = lambda pid: nvalues[pid] )
    #pids = sorted(masses.keys(), key = lambda pid: np.sum(np.where(masses[pid][::nsteps] <= 0.)))
    ctentries=0
    K=Ks[firststep]
    if K < 0.:
        K = 0.
    if K > maxK:
        maxK = K
        stepatmax = firststep
    for pid in pids:
        if max(masses[pid][firststep:laststep:nsteps]) <= 0.0: continue
        ctentries+=1
        data = df[firststep:laststep:nsteps]
        datamax = df[stepatmax:stepatmax+1]
        tName = r'$%s$' % namer.texName(pid)
        c = colorDict[pid]
        #if ctentries>9:
        #        tName=""
        sns.scatterplot(x=data['step'],y=data[pid], size=data['K'], sizes = (80,400),
                        label= tName , color=c, legend=False, ax=ax2 )
        s= (1+1.25*maxK)*80. ## no idea why
        sns.scatterplot(x=datamax['step'],y=datamax[pid], s=s, sizes = (80,400),
                        label= tName , color=c, legend=False, ax=ax1 )
        m = np.where(masses[pid] > 0, masses[pid],np.nan) #Fix for avoid plotting to negative values
        plt.plot(df['step'][firststep:laststep:nsteps],m[firststep:laststep:nsteps],'-',linewidth=2, color = c)

    plt.ylim(0.,2500.0)
    #plt.title ( '$K_{max}=%.1f, K_{current}=%.1f$' % ( maxK, K ), loc="left" )
    ax2.set_title ( ' $\;K_{cur}=%.1f$' % ( K ), fontsize=18, pad=15., loc="left" )
    ax1.set_title ( '$K_{max}=%.1f$' % ( maxK ), fontsize=18, pad=15. )
    ax2.set_xlabel('step', fontsize=20, labelpad=-3. )
    ax1.set_xlabel('', fontsize=21 )
    ax1.set_xticks([])
    plt.ylabel('Mass [GeV]', fontsize=21 )
    ax2.set_ylabel('Mass [GeV]', fontsize=21 )
    ax1.set_ylabel('Mass [GeV]', fontsize=21, labelpad=-5. )
    dstep = 10 ## 2
    #plt.xticks(df['step'][firststep:laststep:2*nsteps])
    nextstep = math.ceil((firststep+1)/dstep)*dstep-1
    plt.xticks(df['step'][nextstep:laststep:dstep])
    # plt.xlim(-5,198)
    plt.grid(axis='x') 
    ac=actions[firststep]
    lac = len(ac)
    while len(ac)<3:
        ac.append ( "" )
    if lac>0:
        ss = math.ceil ( len(ac)/3 )
        txt="\n".join(ac[::ss])
        plt.text ( -8+firststep, -320, txt, c="gray", size=8 )
    bc=bcs[firststep]
    lbc=len(bc)
    #while len(bc)<6:
    #    bc.append( "" )
    if lbc>0:
        ss = math.ceil ( len(bc)/6 )
        txt="\n".join([ x[:x.find(":")] for x in bc[::ss] ] )
        plt.text ( 20.5+firststep, 10, txt, size=10, horizontalalignment="right", verticalalignment="bottom", c="gray" )
        
    plt.legend(loc=(.6,.7),# bbox_to_anchor=(0.6,0.5,.2,.25), 
               framealpha=1.0,ncol=3,labelspacing=0.1,
               handlelength=0.4,handletextpad=0.35,markerscale=0.8,columnspacing=1.0)
    # plt.tight_layout()
    step = firststep + 1 ## make all one-indexed, ok?
    plt.savefig('walk%.3d.png' % step, dpi=300 )
    # plt.show()
    plt.clf()

cmd = 'ffmpeg -y -i "walk%3d.png" walk.mp4'
subprocess.getoutput ( cmd )
