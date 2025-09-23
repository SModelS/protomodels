#!/usr/bin/env python3

""" summarize the data in hiscores.dict to have an overview """

import os, time
from os import PathLike
from ptools.sparticleNames import SParticleNames 
from smodels_utils.helper.terminalcolors import *
from typing import Union
import numpy as np

def summarizeJobsForThisDir():
    from clip.cliphelpers import getJobStatus, readJobIds
    jobids = readJobIds()
    statuses = getJobStatus ( jobids )
    reverse = {}
    for k,v in statuses.items():
        if not v in reverse:
            reverse[v]=set()
        reverse[v].add(k)
    nrunning = 0
    if "running" in reverse:
        nrunning = len(reverse["running"])
    print ( f"In this directory: {GREEN}{nrunning}{RESET} running jobs, {YELLOW}{len(statuses)}{RESET} total" )
    print ( )
    return 1

def summarizeHiscores ( dictfile : PathLike = "hiscores.dict",
    extended : bool = False, nmax : Union[None,int] = None ) -> int:
    """ summarize the content of the dict file 

    :param dictfile: path to dictionary file
    :param extended: extended output, add description timestamp
    """
    printTruth()
    nlines = 0
    if not os.path.exists ( dictfile ):
        print ( f"[printSimpleHiscoreList] {dictfile} does not exist" )
        nlines += 1
        return nlines
    with open( dictfile, "rt" ) as f:
        txt=f.read().replace('"inf"',"float('inf')").replace('"nan"',"float('nan')")
        txt=txt.replace("'inf'",'float("inf")').replace("'nan'",'float("nan")')
        f.close()
        try:
            D=eval(txt)
        except SyntaxError as e:
            print ( f"[printSimpleHiscoreList.summarizeHiscores] could not read {dictfile}: {e}" )
            print ( f"  message {e.msg}" )
            print ( f"  line number {e.lineno}" )
            print ( f"  offset {e.offset}" )
            print ( f"  text {e.text}" )
            nlines += 1
            return nlines
    if nmax == None:
        nmax = 10
        if extended:
            nmax = 3
    for i,entry in enumerate ( D ):
        if extended and i >= nmax:
            break
        wid = 0
        K, TL = entry['K'], entry['TL']
        if "walkerid" in entry:
            wid = entry['walkerid']
        particles = entry["masses"].keys()
        sparticles = ""
        for ip, p in enumerate ( particles ):
            if ip != 0:
                sparticles += ", "
            name = SParticleNames( False).asciiName(p)
            mass = entry["masses"][p]
            sparticles += f"{CYAN}{name}{RESET}={mass:.1f}"
        timestamp = ""
        if "timestamp" in entry:
            timestamp = entry["timestamp"]
            r1 = timestamp.find(" ")
            r2 = timestamp.rfind(" ")
            timestamp = timestamp[r1:r2]
        if extended:
            step = entry["step"]
            print ( f"#{i}({wid:3d}): K={GREEN}{K:.3f}{RESET} TL={TL:.3f}; {sparticles}" )
            print ( f"       `---: {entry['description']}" )
            print ( f"       `---:{timestamp}" )
            print ( f"       `---: step {step}" )
            print ( )
            nlines += 5
        else:
            sK = "None" if K == None else f"{K:.3f}"
            print ( f"#{i}({wid:3d}): K={GREEN}{sK}{RESET}; TL={TL:.3f}; {sparticles} {timestamp}" )
            nlines += 1
    return nlines

def runSlurmWalk() -> int:
    rundir = os.getcwd()
    cmd = f"slurm_walk.py -R {rundir} -q"
    import subprocess
    o = subprocess.getoutput ( cmd )
    print ( f"{RED}Running Status:{RESET} {time.asctime()}" )
    print ( "=========================================" )
    print ( o )
   #  print ( )
    return 4

def printTruth():
    """ if a truth.dict file is found, pretty print it """
    dictfile = "truth.dict"
    if not os.path.exists ( dictfile ):
        return
    with open( dictfile, "rt" ) as f:
        txt=f.read().replace('"inf"',"float('inf')").replace('"nan"',"float('nan')")
        txt=txt.replace("'inf'",'float("inf")').replace("'nan'",'float("nan")')
        f.close()
    d = eval(txt)
    K, TL = d["K"], d["TL"]
    sK = "None" if K == None else f"{K:.3f}"
    particles = d["masses"].keys()
    sparticles = ""
    for ip, p in enumerate ( particles ):
        if ip != 0:
            sparticles += ", "
        name = SParticleNames( False).asciiName(p)
        mass = d["masses"][p]
        sparticles += f"{CYAN}{name}{RESET}={RED}{mass:.1f}{RESET}"
    print ( f"{RED}Truth:{RESET}   K={GREEN}{sK}{RESET}; TL={TL:.3f}; {sparticles}" )

if __name__ == "__main__":
    import argparse
    argparser = argparse.ArgumentParser(
        description='summarize the data in hiscores.dict to have an overview' )
    argparser.add_argument ( '-H', '--hiscores', type=str,
        help="path to hiscores.dict file [./hiscores_global.dict]", 
        default="./hiscores_global.dict" )
    argparser.add_argument ( '-x', '--extended', action="store_true",
        help="extended info" )
    argparser.add_argument ( '-l', '--loop', action="store_true",
        help="loop" )
    argparser.add_argument ( '-n', '--nmax', type=int, default=None,
        help="print maximally this number of entries [None]" )
    args = argparser.parse_args()
    import colorama
    colorama.init()
    while True:
        nlines = runSlurmWalk()
        nlines += summarizeJobsForThisDir()
        nlines += summarizeHiscores ( args.hiscores, args.extended, args.nmax )
        if not args.loop:
            break
        time.sleep(10.)
        print( colorama.Cursor.UP()*(nlines+2) )
