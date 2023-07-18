#!/usr/bin/env python3
# coding: utf-8

import sys, os, copy, glob, argparse
import numpy as np
sys.path.append(os.path.abspath('../smodels'))
sys.path.append(os.path.abspath('../'))
from builder.protomodel import ProtoModel
from builder.manipulator import Manipulator
from tester.predictor import Predictor
from tester.combiner import Combiner
from walker.hiscore import Hiscore
from smodels.experiment.databaseObj import Database
from smodels.tools import runtime
from smodels.base.physicsUnits import fb
if False:
    runtime._experimental = True
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
from ptools.sparticleNames import SParticleNames
# sns.set() #Set style
# sns.set_style('ticks')
sns.set_style('ticks',{'font.family':'Times New Roman', 'font.serif':'Times New Roman'})
sns.set_context('paper', font_scale=2.0)
# sns.set_palette(sns.color_palette("Paired"))
sns.set_palette(sns.color_palette("deep"))
from smodels.tools import smodelsLogging
smodelsLogging.setLogLevel('error')

argparser = argparse.ArgumentParser(description="high score comparison plotter")
argparser.add_argument ( '-s', '--signalfiles',
        help='the run files [kinky]',
        type=str, default='kinky' )
argparser.add_argument ( '-r', '--realfiles',
        help='the real files [real]',
        type=str, default='real' )
args=argparser.parse_args()

runname = args.signalfiles
realname=args.realfiles

#Set colors:
allPids = [1000022, 1000006, 1000001, 1000021, 1000012, 1000023, 
            1000013, 2000006, 1000011, 1000005, 1000014, 1000004, 1000015, 1000016, 1000024]
namer = SParticleNames ( susy = False )
colors = sns.color_palette('deep',n_colors=len(namer.xIDs))

#Replace default colors:
for pid in sorted(namer.xIDs.keys()):
    if not pid in allPids:
        allPids.append(pid)
colorDict = dict(zip(allPids,colors))
colorDict[1000002] = colorDict[1000001]
colorDict[1000003] = colorDict[1000001]
colorDict[1000004] = colorDict[1000001]
colorDict[1000013] = colorDict[1000011]
colorDict[1000015] = colorDict[1000011]


# In[3]:


def getLikelihoods(protomodel,muvals,normalize=True):
    
    #Sort tpList
    critic = sorted(protomodel.tpList, key = lambda tp: tp[0], reverse = True)[0][2]

    #Combiner likelihood:
    combiner = Combiner(0)
    llhComb = np.array([combiner.getCombinedLikelihood(protomodel.bestCombo,mu) for mu in muvals])
    llhCombSM = combiner.getCombinedLikelihood(protomodel.bestCombo,0.)
    llhDict = {'Combined' : llhComb}
    llhDictSM = {'Combined' : llhCombSM}
    for tp in protomodel.bestCombo:
        llhDict[tp.expResult.globalInfo.id] = np.array([tp.getLikelihood(mu) for mu in muvals])
        llhDictSM[tp.expResult.globalInfo.id] = tp.getLikelihood(0.0)

    #Critic likelihood:
    if critic.getLikelihood(1.0) is not None:
        llhDictSM['Critic'] = critic.getLikelihood(0.0)
        llhDict['Critic'] = np.array([critic.getLikelihood(mu) for mu in muvals])
    else:
        llhDict['Critic'] = None
        llhDictSM['Critic'] = None
    llhDict['SM values'] = llhDictSM
    
    #Compute normalizations:
    if normalize:
        for key,llhd in llhDict.items():
            if key != 'SM values':
                norm = llhd.sum()
                llhDict[key] = llhd/norm
                llhDict['SM values'][key] *= 1/norm
    
    
    return llhDict


# In[4]:


def fromDict(inputDict):
    
    p = ProtoModel(walkerid=0)
    for key,v in inputDict.items():
        setattr(p,key,copy.deepcopy(v))
        
    return p


#Get highest score from each run:
protomodelsDict = {}
files = f'../data/{runname}*.dict'
for ff in glob.glob( files ):
    if runname == "narrow" and "fudged" in ff:
        continue
    with open(ff,'r') as f:
        pList = eval(f.read())
    run = eval(os.path.basename(ff).replace(runname,'').replace('.dict',''))
    pList = [fromDict(pDict) for pDict in pList[:]]
    p = sorted(pList, key = lambda p: p.K, reverse=True)[0]
    protomodelsDict[run] = p  
    
protomodelsDictReal = {}
# realwinner = 9
realwinner, winnerK = -1, -1
print ( "realname", realname )
for ff in glob.glob(f'../data/{realname}*.dict'):
    with open(ff,'r') as f:
        pList = eval(f.read())
    run = eval(os.path.basename(ff).replace(realname,'').replace('.dict',''))
    pList = [fromDict(pDict) for pDict in pList[:]]
    p = sorted(pList, key = lambda p: p.K, reverse=True)[0]
    protomodelsDictReal[run] = p
    if p.K > winnerK:
        winnerK = p.K
        realwinner = int(ff.replace(f"../data/{realname}","").replace(".dict",""))

if realwinner < 0:
    print ( "why did we not find a winner???" )

for run in sorted(protomodelsDict.keys()):
    print(run,protomodelsDict[run])
    
Kavg = np.array([p.K for p in protomodelsDict.values()]).mean()
Kstd = np.array([p.K for p in protomodelsDict.values()]).std()
print('K (avg) = %1.2f +- %1.2f' %(Kavg,Kstd))


# In[ ]:


#Get all particles which appears in all models:
particles = []
modelList = list(protomodelsDict.items())
modelList = sorted(modelList, key = lambda pt: pt[0])
modelList = [[-0.75,protomodelsDictReal[realwinner]]] + modelList[:] #Add winning protomodel as run 0
modelList = np.array(modelList)
runs = modelList[:11,0]
modelList = modelList[:11,1]
for p in modelList:
    particles += p.unFrozenParticles()
particles = list(set(particles))

#Build useful dataset:
nparticles = np.array([len(p.unFrozenParticles()) for p in modelList])
Kvalues = np.array([p.K if (p.K and p.K > 0) else 0.0 for p in modelList])
Zvalues = np.array([p.Z if (p.Z and p.Z > 0) else 0.0 for p in modelList])
masses = dict([[pid,[]] for pid in particles])
for p in modelList:
    for pid in masses:
        if pid in p.masses:
            masses[pid].append(p.masses[pid])
        else:
            masses[pid].append(-100.0)
for pid in masses:
    masses[pid] = np.array(masses[pid])
dataDict = {'run' : runs, 'K' : Kvalues,
                   'nparticles' : nparticles}
dataDict.update(masses) 
df = pd.DataFrame(dataDict)


# In[ ]:


f, axarr = plt.subplots(2,sharex=True, gridspec_kw = {'height_ratios':[1, 4]},figsize=(12,8))
plt.subplots_adjust(left=0.12, bottom=0.12, right=0.97, top=None, wspace=None, hspace=0)

nsteps = 10


axarr[0].scatter(df['run'],df['K'],s=80,c='gray')
axarr[0].set_ylabel(r'$K$')
axarr[0].set_ylim(1.0,15.0)
axarr[0].set_yticks([])
for i,row in df.iterrows():
    axarr[0].annotate(r'$%1.2f$' %row['K'],(row['run']-0.2,row['K']+0.8),fontsize=13)

axarr[0].vlines(x=0,ymin=2,ymax=15,linestyle='--',color='gray')

for pid in masses.keys():
    if not pid in masses: #Skip particles present in real run, but not the fakes
        continue
    data = df
    sns.scatterplot(x=data['run'],y=data[pid], size=1000,sizes=(1500,1500),marker='_',
                    label=r'$%s$' %(namer.texName(pid,addOnes=True)), legend=False,
                    color=[colorDict[pid]],ax=axarr[1])
    for i,m in enumerate(masses[pid]):
        if m < 0: continue
        if 'b' in namer.texName(pid):
            axarr[1].annotate(r'$%s$' %(namer.texName(pid,addOnes=True)),(runs[i]-0.3,m+25.),fontsize=15)
        else:
            axarr[1].annotate(r'$%s$' %(namer.texName(pid,addOnes=True)),(runs[i],m+25.),fontsize=15)
axarr[1].set_ylim(0.,1500.0)
axarr[1].set_xlabel('run (BSM)', fontsize=23)
axarr[1].set_ylabel('Mass [GeV]', fontsize=23)
axarr[1].set_xticks(np.arange(1,df['run'].max()+1))
axarr[1].vlines(x=0,ymin=0,ymax=1500,linestyle='--',color='gray')
#     plt.plot(df['run'],m,'-',linewidth=2)

axarr[1].text(-0.7,900,'highest\nscore\n(real run)',fontsize=15,horizontalalignment='center')
plt.ylim(0.,1500.0)

# plt.grid(axis='x') 
# plt.legend(loc=(0.8,0.85),framealpha=1.0,ncol=3,labelspacing=0.1,
#            handlelength=0.4,handletextpad=0.35,markerscale=0.8,columnspacing=1.0)
# plt.tight_layout()
plt.savefig( f'highScoreSignal_{args.signalfiles}.pdf')
plt.savefig( f'highScoreSignal_{args.signalfiles}.png')
print ( f"saving to highScoreSignal_{args.signalfiles}.png highScoreSignal_{args.signalfiles}.pdf" )
# plt.show()
