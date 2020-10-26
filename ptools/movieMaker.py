#!/usr/bin/env python3
import warnings
warnings.filterwarnings("ignore")

import sys, os, pickle, glob, subprocess, math, time, copy, multiprocessing, gc
import numpy as np
import IPython
sys.path.append(os.path.abspath('../../smodels'))
sys.path.append(os.path.abspath('../'))
#from smodels.tools import runtime
#runtime._experimental = True
import matplotlib.pyplot as plt
import matplotlib.patches as patches
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

import argparse
argparser = argparse.ArgumentParser( description="movie maker" )
argparser.add_argument ( '-f', '--history',
        help='history file to use [history.list]',
        type=str, default="history.list" )
argparser.add_argument ( '-o', '--outfile',
        help='output file name [walk.webm]',
        type=str, default="history.list" )
argparser.add_argument ( '-m', '--maxsteps',
        help='maximum steps [1000]',
        type=int, default=1000 )
argparser.add_argument ( '-n', '--start',
        help='step to start with [0]',
        type=int, default=0 )
argparser.add_argument ( '-D', '--dont_clean',
        help='dont clean up old files',
        action="store_true" )
argparser.add_argument ( '-t', '--timestamp',
        help='put timestamp on pic',
        action="store_true" )
args = argparser.parse_args()
prefix = args.outfile.replace(".mp4","").replace(".webm","")
    
intermediateSteps = True ## do 10 rendering steps per one random walk step

if not args.dont_clean:
    cmd = f'rm -f {prefix}*.png'
    subprocess.getoutput( cmd )

#Set colors:
namer = SParticleNames ( susy = False )
#Replace default colors:
colorPalette = 'deep' #Set color palette for particles, if None use default
allpids = [ 1000001, 1000002, 1000003, 1000004, 1000005, 1000006,
2000005, 2000006, 1000011, 1000012, 1000013, 1000014, 1000015,
1000016, 1000021, 1000022, 1000023, 1000025, 1000024, 1000037 ]
colorDict = dict(zip( allpids,sns.color_palette(palette=colorPalette,n_colors=len(namer.names))))

maxstep = args.maxsteps

f=open(args.history,"rt")
lines=f.readlines()
f.close()
txt=""
for line in lines:
    if line.startswith("#"):
        continue
    #line = line.replace("nan","'nan'")
    #line = line.replace("'nan''",'"nan"\'')
    line = line.replace("nan","-300.") ## FIXME
    line = line.replace("'nan''","-300.'") ## ugly!!
    txt+=line+"\n"
if "[" in txt[-3:]:
    txt=txt[:-3]
if not "]" in txt[-3:]:
    txt+="]\n"

modelList=eval(txt)

emptymodel = { "masses": {}, "step": 0, "bestCombo": [], "actions": [], "K": -200., "Z": -200. }
nstart=0
for i in range(19):
    nstart+=1
    em = copy.deepcopy( emptymodel )
    em["step"]=0-i
    modelList.insert(0,em)

while False: # len(modelList)<maxstep+21:
    lastModel = modelList[-1]
    lm = copy.deepcopy ( lastModel )
    lm["step"]=lm["step"]+1
    modelList.append ( lm )

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
        for a in p["actions"]:
            t=a.replace("Unfreeze","unfreeze")
            t=t.replace("Freeze","freeze")
            t=t.replace("->","$\\rightarrow$")
            ac.append ( t )
    def sortMsgs ( x ):
        if "teleport" in x:
            return 1
        if "unfreeze" in x:
            return 2
        if "freeze" in x:
            return 3
        if "mass" in x:
            return 4
        if "decay" in x:
            return 4
        if "ssm" in x:
            return 100
        return 0
    ac.sort ( key=sortMsgs )
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

maxK,stepatmax=-90.,0

currentstep = 8

if maxstep > len(Ks)-currentstep-20:
    maxstep=len(Ks)-currentstep-20

print ( "[movieMaker] setting maxstep to %d" % maxstep )

style = "Simple, tail_width=0.5, head_width=4, head_length=8"
kw = dict(arrowstyle=style, color="k")

cmd = f'ffmpeg -y -i "{prefix}%3d.png" -filter:v "setpts=6.0*PTS" {args.outfile}'
if intermediateSteps:
    cmd = f'ffmpeg -y -r 200 -i "{prefix}%5d.png" {args.outfile}'
    #cmd = f'ffmpeg -y -i "{prefix}%5d.png" {args.outfile}'
print ( "the command for movie making will be:" )
print ( cmd )

imgnr=0
maxHS = 19

def getHiscore ( firststep, currentstep, Ks ):
    """ get the hiscore at step firststep-1 """
    maxK = -90.
    stepatmax = 0
    for i in range(firststep+currentstep):
        K=Ks[i]
        if K > maxK:
            maxK = K
            stepatmax = i
            # model = modelList[i]
    return (maxK,stepatmax)

def onePic ( firststep, offs, maxK, masses, pids, lastingHS, stepatmax, imgnr, K ):
    """ make a single picture """
    fig, (ax1, ax2) = plt.subplots( ncols=2, sharey=True, 
                      gridspec_kw={'width_ratios': [1, 10]} )
    ctentries=0
    plt.text ( -3+firststep-nstart+offs, 1250, "hiscore", rotation=90., c="pink", alpha=.5, 
               size=30, horizontalalignment='center', verticalalignment='center', zorder=5 )
    if K > maxK:
        lastingHS = maxHS ## keep it for 9 frames
        stepatmax = firststep+currentstep
        maxK = K
    if lastingHS>0:
        maxK = K
        red=(1., (maxHS-lastingHS)/maxHS, (maxHS-lastingHS)/maxHS )
        plt.text ( .5+firststep-nstart+offs, 2200, "hiscore!", c=red, size=30, clip_on=False )
        lastingHS-=1 ## count down to zero
    # arrow = plt.arrow ( currentstep, -50, -currentstep-2, 0, width=20, clip_on=False, transform=ax2.transData, color="black" )
    # arrow = patches.FancyArrowPatch((.1, .5), (.1, .5), clip_on=False,
    #arrow = patches.FancyArrowPatch((firststep+currentstep-5-nstart, firststep-nstart), ( 50, 50), 
    #    clip_on=False, connectionstyle="arc3,rad=.1", **kw, transform=ax2.transData )
    # plt.gca().add_patch ( arrow )

    for pid in pids:
        if max(masses[pid][firststep:laststep:nsteps]) <= 0.0: continue
        ctentries+=1
        data = df[firststep:laststep+1:nsteps]
        datamax = df[stepatmax:stepatmax+1]
        datacur = df[firststep+currentstep:firststep+currentstep+1]
        tName = r'$%s$' % namer.texName(pid)
        c = colorDict[pid]
        #if ctentries>9:
        #        tName=""
        m = np.where(masses[pid] > 0, masses[pid],np.nan) #Fix for avoid plotting to negative values
        ## the lines
        plt.plot(df['step'][firststep:laststep+2:nsteps],m[firststep:laststep+2:nsteps],'-',linewidth=2, color = c, alpha=.5 )
        sizes = []
        def getSize ( k ):
            if k < 0.: k = 0.
            if k in [ None, float("nan"), float("inf") ]:
                k = 0.
            return (1+.6*k)*80
        for s in data["K"]:
            sizes.append( getSize ( s ) )
        sns.scatterplot(x=data['step'],y=data[pid], s=sizes, sizes = (80,400),
                        label= tName, color=c, legend=False, ax=ax2, alpha=.5 )
        # s= (1+1.25*maxK)*80. ## no idea why
        smax= getSize ( maxK ) ## no idea why
        s = getSize ( K ) ## no idea why
        sns.scatterplot(x=datacur['step'],y=datacur[pid], s=1.2*s+40, sizes = (80,400),
                        label= "", color="black", legend=False, ax=ax2, edgecolor="none",
                        linewidth=0 )
        sns.scatterplot(x=datacur['step'],y=datacur[pid], s=s, sizes = (80,400),
                        label= "", linewidth=0, edgecolor="none", color=c, 
                        legend=False, ax=ax2 )
        sns.scatterplot(x=datamax['step'],y=datamax[pid], s=smax, sizes = (80,400),
                        label= tName , color=c, legend=False, ax=ax1, zorder=10)
    plt.ylim(0.,2500.0)
    title='$K_{current}=%.1f\;\,$   $\;$ $\;\;\;\;\;\,$ $\;$ $\;$ $\;$ ' % ( K )
    if K < -50.:
        title=""
    ax2.set_title ( title, fontsize=18, pad=15., horizontalalignment="center" )
    ax1title = '$K_{max}=%.1f$' % ( maxK )
    if maxK < -50.:
        ax1title = ""
    ax1.set_title ( ax1title, fontsize=18, pad=15. )
    ax2.set_xlabel('step', fontsize=20, labelpad=-3. )
    ax1.set_xlabel('', fontsize=21 )
    ax1.set_xticks([])
    plt.ylabel('Mass [GeV]', fontsize=21 )
    ax2.set_ylabel('Mass [GeV]', fontsize=21 )
    ax1.set_ylabel('Mass [GeV]', fontsize=21, labelpad=-5. )
    dstep = 10 ## 2
    #plt.xticks(df['step'][firststep:laststep:2*nsteps])
    nextstep = math.ceil((firststep+1)/dstep)*dstep-1
    ticks=[]
    for t in df['step'][nextstep:laststep:dstep]:
        if t < 1:
            continue
            # t = 0
        ticks.append ( t )
    plt.xticks( ticks )
    # plt.xticks(df['step'][nextstep:laststep:dstep])
    # plt.xlim(-5,198)
    plt.grid(axis='x')
    ac=actions[firststep+currentstep]
    lac = len(ac)
    while len(ac)<3:
        ac.append ( "" )
    if lac>0:
        ss = math.ceil ( len(ac)/3 )
        txt="\n".join( ac[:3] )
        # txt="\n".join([ x.replace("->","$\\rightarrow$") for x in ac[::ss] ])
        plt.text ( -8+firststep-nstart+offs, -320, txt, c="gray", size=7 )
    bc=bcs[firststep+currentstep]
    lbc=len(bc)
    #while len(bc)<6:
    #    bc.append( "" )
    if lbc>0:
        ss = math.ceil ( len(bc)/6 )
        txt="\n".join([ x[:x.find(":")] for x in bc[::ss] ] )
        bbox = None
        bbox=dict(boxstyle="round", ec=(1., 0.5, 0.5), fc=(1., 0.8, 0.8),)
        plt.text ( 18+firststep-nstart, 10, txt, size=10, horizontalalignment="right", verticalalignment="bottom", c="gray", bbox=bbox, rotation=10. )

    plt.legend(loc=(.6,.7),# bbox_to_anchor=(0.6,0.5,.2,.25),
               framealpha=1.0,ncol=3,labelspacing=0.1,
               handlelength=0.4,handletextpad=0.35,markerscale=0.8,columnspacing=1.0)
    # plt.tight_layout()
    if args.timestamp:
        plt.text ( 15+firststep-nstart+offs, -280, time.asctime(), size=8, 
                   alpha=.5, c="gray" )
    step = firststep + 1 ## make all one-indexed, ok?
    off1 = firststep-nstart+.05+offs
    ax2.set_xlim ( off1, off1 + 21 )
    if intermediateSteps:
        plt.savefig('%s%.5d.png' % \
                ( prefix, imgnr ), dpi=200 )
    else:
        plt.savefig('%s%.3d.png' % (prefix, step ), dpi=200 )
    # plt.show()
    plt.clf()
    plt.close("all")
    gc.collect()
    ret = { "lastingHS": lastingHS, "maxK": maxK, "stepatmax": stepatmax }
    return ret

alloffs = [ 0. ]
if intermediateSteps:
    alloffs = np.arange(0,1.,.025)

if args.start>0:
    maxK,stepatmax = getHiscore ( args.start, currentstep, Ks )
    imgnr = len(alloffs) * args.start

for firststep in range ( args.start, maxstep ):
    if firststep % 10 == 0:
        print ( "step %d: %s" % ( firststep, time.asctime() ) )
    lastingHS = 0 ## the "hiscore!" label should last a bit
    
    laststep=firststep+20
    K=Ks[firststep+currentstep]

    nvalues = {}
    for k,v in masses.items():
        nvalues[k]=0
        for i in v:
            if float(i)>0:
                    nvalues[k]+=1
    pids = sorted(masses.keys(), key = lambda pid: nvalues[pid] )

    for offs in alloffs:
        ret = onePic ( firststep, offs, maxK, masses, pids, lastingHS, stepatmax, imgnr, K )
        lastingHS = ret["lastingHS"]
        maxK = ret["maxK" ]
        stepatmax = ret["stepatmax"]
        imgnr+=1

subprocess.getoutput ( cmd )
