#!/usr/bin/env python3

""" a bit of code to help us with parsing the dicts of the runs """

import glob, subprocess, sys, os, colorama
import numpy as np

def run( filename ):
    f = open ( filename, "rt" )
    txt = f.read()
    txt = txt.replace("nan","'nan'")
    D=eval( txt )
    f.close()
    return D

def produce():
    for tps in [ "real", "fake", "signal" ]:
        files = glob.glob ( f"{tps}*.dict" )
        files.sort()
        for fname in files:
            D = run ( fname )
            P = []
            for k,v in D[0]["masses"].items():
                if v < 5000. and k != 1000022:
                    P.append(k)
            P.sort()
            print ( fname, "K", D[0]["K"], "Z", D[0]["Z"], "particles", P )


def countSteps( Dir ):
    """ count the steps in Dir """
    import glob
    files = glob.glob(f"{Dir}/walker*log")
    steps = {}
    for f in files:
        nr = int ( f.replace("walker","").replace(".log","").replace(Dir,"").replace("/","") )
        h = open ( f, "rt" )
        lines = h.readlines()
        h.close()
        for line in lines[::-1]:
            if "Step" in line:
                laststep = line[line.find("Step")+5:]
                for c in [ "/", ":", " has", " " ]:
                    if c in laststep:
                        laststep = laststep[:laststep.find(c)]
                try:
                    laststep = int ( laststep.strip() )
                except Exception as e:
                    print ( "Error", e )
                    print ( "line", line )
                    print ( "file", f )
                    sys.exit()

                #print ( nr, laststep )
                steps[nr]=laststep
                break
    keys = list ( steps.keys() )
    keys.sort()
    tots = 0
    for k in keys:
        tots += steps[k]
    return tots

def count():
    """ count the total number of steps take in a run """
    Dir = "/scratch-cbe/users/wolfgan.waltenberger/"
    Dirs = glob.glob ( f"{Dir}/rundir.*/hiscore.hi" )
    Dirs.sort()
    for d in Dirs:
        wdir = d.replace("/hiscore.hi","")
        nsteps = countSteps ( wdir )
        if nsteps < 50000:
            print ( wdir, nsteps )

def fetch():
    """ fetch states.dict files from the individual runs """
    Dir = "/scratch-cbe/users/wolfgan.waltenberger/"
    dictfiles = "hiscores.dict"
    # dictfiles = "states.dict"
    files = glob.glob ( f"{Dir}/rundir.*/{dictfiles}" )
    hiscores={}
    for f in files:
        if "real" in f:
            print ( "suppressing copying reals for now!!!" )
        name = f.replace( Dir, "" ).replace("/"+dictfiles,"").replace("rundir.","")
        cmd = f"cp {f} {name}.dict"
        with open ( f, "rt" ) as h:
            D=eval(h.read())
            hiscores[name]=D[0]
        # print ( cmd )
        subprocess.getoutput ( cmd )

    """
    statesfiles = "states.dict"
    files = glob.glob ( f"{Dir}/rundir.*/{statesfiles}" )
    states={}
    for f in files:
        name = f.replace( Dir, "" ).replace("/"+statesfiles,"").replace("rundir.","")
        with open ( f, "rt" ) as h:
            D=eval(h.read())
            states[name]=D[0]
    print ( "when fetching:" )
    for name in hiscores.keys():
        if not "real" in name:
            continue
        if name in states:
            print ( "fetch:", name,hiscores[name]["K"],states[name]["K"] )
        else:
            print ( "fetch:", name,hiscores[name]["K"],"absent!!" )
    """

def belowMassWall ( D, name ):
    """ check if anything is below the mass wall """
    walledpids = [ 1000001, 1000002, 1000003, 1000004, 1000021 ]
    wallmass = 310.
    masses = D["masses"]
    ctr=0
    for pid in walledpids:
        if pid in masses and masses[pid] < wallmass:
            print ( f"{name} has {pid} below mass wall: m={masses[pid]}" )
            ctr+=1
    return ctr
            

def inCorridor ( D, name ):
    """ check if the model contained in dictionary D has a stop in the corridor """
    masses = D["masses"]
    if not 1000006 in masses and not 2000006 in masses:
        return 0
    ctr=0
    for mpid in [ 1000006, 2000006 ]:
        if mpid in masses:
            dm = masses[mpid]-masses[1000022]
            if 150 < dm < 200:
                print ( f"{name} has {mpid} in corridor: dm={dm}" )
                ctr+=1
    return ctr

def analyzeStats ( globber ):
    print ( f"analyzing {globber}.dict" )
    files = glob.glob ( f"{globber}.dict" )
    pids = {}
    Ks = {}
    inc, bmw = 0, 0
    for f in files:
        with open ( f, "rt" ) as h:
            lines = h.read()
        lines = lines.replace("nan","'nan'")
        D = eval(lines)
        Ks[ D[0]["K"] ] = f
        for ps in D[0]["masses"].keys():
            if not ps in pids:
                pids[ps]=0
            pids[ps]+=1
        inc += inCorridor ( D[0], f )
        bmw += belowMassWall ( D[0], f )
    maxK = 0
    col, res = "", ""
    if inc + bmw > 0:
        col = colorama.Fore.RED
        res = colorama.Fore.RESET
    print ( f"of {len(files)} files: {col}{inc} in corridor, {bmw} below mass wall{res}" )
    if len(Ks)>0:
        maxK = max(Ks.keys())
        print ( "winner is", Ks[ maxK ],"with K=", maxK )
        print ( "Ks were between", min(Ks.keys()), np.mean(list(Ks.keys())), maxK )
    print ( "particles that are in", pids )

def addAsterisk ( pattern ):
    wcs =  [ "*", "?" ]
    wcs += list(map(str,range(0,10)))
    addAsterisk=True
    for char in wcs:
        if char in pattern:
            addAsterisk=False
    if addAsterisk:
        pattern = pattern + "*"
    return pattern

def getBest( pattern ):
    """ find out best scoring models """
    files = glob.glob ( f"{pattern}.dict" )
    Ks = {}
    maxK, maxf = 0, None
    for f in files:
        with open ( f, "rt" ) as h:
            txt = h.read()
            txt = txt.replace("nan","'nan'")
            L = eval( txt )
            Ks[f]=L[0]["K"]
            if Ks[f] > maxK:
                maxK=Ks[f]
                maxf = f
    print ( f"getBest {pattern}: {Ks}" )
    print ( f"the winner is {maxf} with K={maxK}" )
    return Ks

def defineModel1 ( ):
    """ define the winning model, copy it to $RUNDIR/pmodel1.py """
    Ks = getBest ( "real*" )
    bestf,bestK = None, 0.
    for k,v in Ks.items():
        if v > bestK:
            bestK = v
            bestf = k
    rundir = os.environ["RUNDIR"]
    cmd = "cp %s %s/%s" % ( bestf, rundir, "pmodel1.py" )
    with open ( bestf, "rt" ) as r:
        lines = r.readlines()
    with open ( f"{rundir}/pmodel1.py", "wt" ) as f:
        f.write ( lines[0][1:-2]+ "\n" )
        f.close()
    #print ( cmd )
    #subprocess.getoutput ( cmd )
    print ( "now run in ~/git/protomodels/ptools:" )
    print ( "S=signal; for i in `seq 1 19`; do ./expResModifier.py -R %s.${S}${i} -d %s/original.pcl -s ${S}${i} -P %s/pmodel1.py; ./expResModifier.py -R %s.${S}${i} -d %s.${S}${i}/${S}${i}.pcl --remove_orig --nofastlim --onlyvalidated --nosuperseded --dontsample -s ${S}${i} --symlink -o %s.${S}${i}/filtered${i}.pcl; done" % ( rundir, rundir, rundir, rundir, rundir, rundir ) )
    #print ( "for i in `seq 1 19`; do ./expResModifier.py -R %s.signal${i} -d %s/original.pcl -s signal${i} -P %s/pmodel1.py; ./expResModifier.py -R %s.signal${i} -d %s.signal${i}/signal${i}.pcl --remove_orig --nofastlim --onlyvalidated --nosuperseded --dontsample -s signal${i} --symlink -o %s.signal${i}/filtered${i}.pcl; done" % ( rundir, rundir, rundir, rundir, rundir, rundir ) )
    print ( )

if __name__ == "__main__":
    import argparse
    argparser = argparse.ArgumentParser(description="summarize all runs" )
    argparser.add_argument ( '-c', '--count',
               help='count the number of steps already taken', action='store_true' )
    argparser.add_argument ( '-d', '--define',
               help='define signal model', action='store_true' )
    argparser.add_argument ( '-p', '--produce',
               help='produce the stats', action='store_true' )
    argparser.add_argument ( '-a', '--analyze',
               help='analyze the stats', action='store_true' )
    argparser.add_argument ( '-g', '--globber',
               help='globber to use for stats analysis', type=str, default="real*" )
    args=argparser.parse_args()

    if args.define:
        defineModel1()
    if args.count:
        count()
    if args.produce:
        fetch()
        produce()
        getBest( "real*" )
        getBest( "signal*" )
    if args.analyze:
        analyzeStats ( addAsterisk ( args.globber ) )
    
