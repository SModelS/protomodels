#!/usr/bin/env python3
import warnings
warnings.filterwarnings("ignore")

import sys,os,copy, pickle, glob, subprocess
import numpy as np
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
sns.set_context('paper', font_scale=2.0)
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

f=open("history.list","rt")
txt=f.read()
txt=txt.replace("nan","'nan'")
f.close()

modelList=eval(txt)

#Get all particles which appears in all steps:
particles = []
for p in modelList:
    particles += p["masses"].keys()
particles = list(set(particles))
print ( "particle", particles )

#Build useful dataset:
steps = np.array([p["step"] for p in modelList])
nparticles = np.array([len(p["masses"]) for p in modelList])
Kvalues = np.array([p["K"] if (p["K"] and p["K"] > 0) else 0.0 for p in modelList])
Zvalues = np.array([p["Z"] if (p["Z"] and p["Z"] > 0) else 0.0 for p in modelList])
masses = dict([[pid,[]] for pid in particles])
Ks=[]
for p in modelList:
    K=p["K"]
    if K == None:
        K=0.
    Ks.append(K)
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
maxstep = 500
#maxstep = 10
#maxstep = 200

for firststep in range ( maxstep ):
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
    for pid in pids:
        if max(masses[pid][firststep:laststep:nsteps]) <= 0.0: continue
        ctentries+=1
        data = df[firststep:laststep:nsteps]
        tName = r'$%s$' % namer.texName(pid)
        c = colorDict[pid]
        #if ctentries>9:
        #        tName=""
        sns.scatterplot(x=data['step'],y=data[pid], size=data['K'], sizes = (80,400),
                        label= tName , color=c, legend=False)
        m = np.where(masses[pid] > 0, masses[pid],np.nan) #Fix for avoid plotting to negative values
        plt.plot(df['step'][firststep:laststep:nsteps],m[firststep:laststep:nsteps],'-',linewidth=2, color = c)

    plt.ylim(0.,2500.0)
    plt.title ( 'K=%.2f' % Ks[firststep], loc="left" )
    plt.xlabel('step', fontsize=23)
    plt.ylabel('Mass [GeV]', fontsize=23)
    plt.xticks(df['step'][firststep:laststep:1*nsteps])
    # plt.xlim(-5,198)
    plt.grid(axis='x') 
    plt.legend(loc=(.8,.75),bbox_to_anchor=(0.6,0.5,.2,.25), 
               framealpha=1.0,ncol=3,labelspacing=0.1,
               handlelength=0.4,handletextpad=0.35,markerscale=0.8,columnspacing=1.0)
    # plt.tight_layout()
    plt.savefig('walk%d.png' % firststep )
    # plt.show()
    plt.clf()

cmd = 'ffmpeg -y -i "walk%d.png" walk.mp4'
subprocess.getoutput ( cmd )
