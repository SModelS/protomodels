#!/usr/bin/env python3

""" summarize the data in hiscores.dict to have an overview """

import os, time
from os import PathLike
from ptools.sparticleNames import SParticleNames 
from colorama import Fore as ansi
from typing import Union
import numpy as np

def summarizeHiscores ( dictfile : PathLike = "hiscores.dict",
    extended : bool = False, nmax : Union[None,int] = None ) -> int:
    """ summarize the content of the dict file 

    :param dictfile: path to dictionary file
    :param extended: extended output, add description timestamp
    """
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
            sparticles += f"{ansi.CYAN}{name}{ansi.RESET}={mass:.1f}"
        timestamp = ""
        if "timestamp" in entry:
            timestamp = entry["timestamp"]
            r1 = timestamp.find(" ")
            r2 = timestamp.rfind(" ")
            timestamp = timestamp[r1:r2]
        if extended:
            step = entry["step"]
            print ( f"#{i}({wid:3d}): K={ansi.GREEN}{K:.3f}{ansi.RESET} TL={TL:.3f}; {sparticles}" )
            print ( f"       `---: {entry['description']}" )
            print ( f"       `---:{timestamp}" )
            print ( f"       `---: step {step}" )
            print ( )
            nlines += 5
        else:
            sK = "None" if K == None else f"{K:.3f}"
            print ( f"#{i}({wid:3d}): K={ansi.GREEN}{sK}{ansi.RESET}; TL={TL:.3f}; {sparticles} {timestamp}" )
            nlines += 1
    return nlines

def runSlurmWalk() -> int:
    cmd = "slurm_walk.py -q"
    import subprocess
    o = subprocess.getoutput ( cmd )
    print ( f"{ansi.RED}Running Status:{ansi.RESET} {time.asctime()}" )
    print ( "=========================================" )
    print ( o )
    print ( )
    return 4

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
        nlines += summarizeHiscores ( args.hiscores, args.extended, args.nmax )
        if not args.loop:
            break
        time.sleep(10.)
        print( colorama.Cursor.UP()*(nlines+2) )
