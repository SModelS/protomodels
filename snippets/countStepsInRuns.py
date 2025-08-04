#!/usr/bin/env python3

from ptools.updateHiscores import countSteps
import glob, os, colorama, time, argparse
import numpy as np

def main():
    argparser = argparse.ArgumentParser(
        description='count the number of finished jobs in protomodels production' )
    argparser.add_argument ( '-p', '--pattern', type=str, default="", 
        help="show only the ones that contain <pattern>" )
    os.chdir ( os.environ['CODEDIR'] )
    args = argparser.parse_args()
    Dirs = glob.glob ( f"rundir.{args.pattern}*/" )
    oldDict = {}
    t0 = time.time()
    oldt = t0
    oldfile = f"{args.pattern}.old"
    if os.path.exists ( oldfile ):
        with open ( oldfile, "rt" ) as f:
            oldDict = eval ( f.read() )
            if "t" in oldDict:
                oldt = oldDict["t"]
    else:
        oldt = t0 - 1000000
    dt = t0 - oldt
    Dirs.sort()
    ntot = 0
    assumed = 0
    njobs = 0
    nfinishedjobs = 0
    Dicts={}
    K=0
    dfound="???"
    Ks=[]
    for d in Dirs:
        if "default" in d or "skeleton" in d:
            continue
        njobs += 50
        os.chdir ( d )
        dictfile = f"hiscores.dict"
        if os.path.exists ( dictfile ):
            try:
                with open ( dictfile, "rt" ) as h:
                    D=eval(h.read())
                    tmp = D[0]["K"]
                    if tmp > K:
                        K = tmp
                        dfound = d
                    Ks.append ( tmp )
            except:
                pass
        n,steps = countSteps( printout=False )
        Dicts[d]=steps
        for k,v in steps.items():
            if v == 1000:
                nfinishedjobs += 1
        ntot += n
        assumed += 50000
        nold = None
        if d in oldDict :
            nold = 0
            for i,v in oldDict[d].items():
                nold += v
        line = f"{d:>16}: {int(n):5}"
        if n == 50000:
            line = f"{colorama.Fore.GREEN}{line}{colorama.Fore.RESET}"
        if n < 50000 and nold == n and dt > 1200:
            line = f"{colorama.Fore.RED}{line}{colorama.Fore.RESET}"
        if n < 50000 and nold == n and dt <= 1200:
            line = f"{colorama.Fore.YELLOW}{line}{colorama.Fore.RESET}"
        if nold != None:
            line += f" (was {int(nold):5})"
        print ( line )
        os.chdir ( os.environ['HOME'] )
    if len(Ks)==0:
        Ks=[0.]
    print ( f"The current {len(Ks)} Ks are at [{min(Ks):.2f},{np.mean(Ks):.2f}+/-{np.std(Ks):.2f},{K:.2f}], best found in {dfound.replace('rundir.', '').replace('/', '')}" )
    print ( "total %d/%d jobs finished: %d%s" % ( nfinishedjobs, njobs, (nfinishedjobs*100/njobs), "%" ) )
    print ( "total %dk/%dk: %d%s" % ( ntot/1000, assumed/1000, ntot/assumed*100, "%" ) )
    Dicts["t"]=time.time()
    if dt > 1800: ## older than 30 minutes?
        with open ( oldfile, "wt" ) as f:
            f.write ( f"{Dicts}\n" )
            f.close()
main()
