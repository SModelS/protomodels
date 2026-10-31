#!/usr/bin/env python3

""" summarize the data in hiscores.dict to have an overview """

import os, time
import sys
from os import PathLike
from typing import Union
import numpy as np
pmgit = f"{os.environ['HOME']}/git/protomodels/"
if os.path.exists ( pmgit ):
    sys.path.insert(0,pmgit )
from ptools.sparticleNames import SParticleNames 
from ptools.helpers import formatObject
from smodels_utils.helper.terminalcolors import *

def summarizeJobsForThisDir():
    from clip.cliphelpers import getJobStatus, readJobIds
    jobids = readJobIds()
    statuses = getJobStatus ( jobids )
    reverse = {}
    for k,v in statuses.items():
        if not v in reverse:
            reverse[v]=set()
        reverse[v].add(k)
    nrunning, npending = 0, 0
    if "running" in reverse:
        nrunning = len(reverse["running"])
    if "pending" in reverse:
        npending = len(reverse["pending"])
    print ( f"In this directory: {GREEN}{nrunning}{RESET} running jobs, {YELLOW}{npending}{RESET} pending jobs, {RED}{len(statuses)}{RESET} total" )
    print ( )
    return 1

def getHiscores ( dictfile : PathLike = "hiscores_global.dict" ) -> list:
    """ get the hiscores
    :param dictfile: if filename, get hiscores from file, if foldername,
    get the best model from each file

    :returns: list of hiscores
    """
    nlines = 0
    if not os.path.exists ( dictfile ):
        print ( f"[printSimpleHiscoreList] {dictfile} does not exist" )
        nlines += 1
        return []
    if os.path.isdir ( dictfile ):
        import glob
        files = glob.glob ( dictfile + "/*dict" )
        files.sort()
        D = []
        for f in files:
            d1 = getHiscores ( f )
            for d in d1:
                D.append ( d )
            #if len(d1)>0:
            #    D.append ( d1[0] )
        return D
    with open( dictfile, "rt" ) as f:
        txt=f.read().replace('"inf"',"float('inf')").replace('"nan"',"float('nan')")
        txt=txt.replace("'inf'",'float("inf")').replace("'nan'",'float("nan")')
        f.close()
        try:
            D=eval(txt)
            return D
        except SyntaxError as e:
            print ( f"[printSimpleHiscoreList.summarizeHiscores] could not read {dictfile}: {e}" )
            print ( f"  message {e.msg}" )
            print ( f"  line number {e.lineno}" )
            print ( f"  offset {e.offset}" )
            print ( f"  text {e.text}" )
            nlines += 1
            return []
    return []

def sortEntries ( D : list, nmax : Union[None,int], 
                  nmax_analysis : Union[None,int] ) -> list:
    """ sort entries in D according to K """
    dc = {}
    for d in D:
        k = d["K"]
        while k in dc:
            k += 1e-6
        dc[k]=d
    keys = list ( dc.keys() )
    keys.sort( reverse=True )
    ret=[]
    walkerids = {} ## count how oftern walkerid is in list
    for i,k in enumerate(keys):
        doAdd = True
        if nmax_analysis not in  [ None, 0 ]:
            walkerid = dc[k]["walkerid"]
            if not walkerid in walkerids:
                walkerids[walkerid]=0
            if walkerids[walkerid]>=nmax_analysis:
                doAdd = False
            walkerids[walkerid]+=1
        if doAdd:
            ret.append ( dc[k] )
        if nmax not in [ None, 0 ] and len(ret) >= nmax:
            break
    return ret

def formatTimestamp ( timestamp : str ) -> str:
    """ format the timestamp of a hiscore entry, and grey it out if the
    entry is recent.

    If a uphi.log file exists, an entry counts as recent if it is more
    recent than the timestamp stored in uphi.log. Otherwise an entry
    counts as recent if it is less than two hours old.

    :param timestamp: the timestamp string of the hiscore entry
    :returns: the (possibly greyed out) timestamp, without the year
    """
    from datetime import datetime
    from zoneinfo import ZoneInfo
    when = datetime.strptime ( timestamp, "%a %b %d %H:%M:%S %Y" )
    r1 = timestamp.find(" ")
    r2 = timestamp.rfind(" ")
    ret = timestamp[r1:r2]
    ret = ret.strip()
    tz = ZoneInfo("Europe/Vienna")
    now = datetime.now(tz)
    when = when.replace(year=now.year, tzinfo=tz)
    ref = None
    if os.path.exists ( "uphi.log" ):
        import ast
        try:
            with open ( "uphi.log", "rt" ) as f:
                ref = float ( ast.literal_eval ( f.read() )["timestamp"] )
        except Exception as e:
            print ( f"[printSimpleHiscoreList.formatTimestamp] could not read uphi.log: {e}" )
    if ref != None:
        recent = when.timestamp() > ref
    else:
        hours_ago = (now - when).total_seconds() / 3600
        recent = hours_ago < 2
    if recent:
        ret = f"{GREY}{ret}{RESET}"
    return ret

def summarizeHiscores ( dictfile : PathLike = "hiscores_global.dict",
    extended : bool = False, nmax : Union[None,int] = None,
    nmax_analysis : Union[None,int] = None ) -> int:
    """ summarize the content of the dict file 

    :param dictfile: path to dictionary file
    :param extended: extended output, add description timestamp

    :returns: number of lines printed
    """
    printTruth()
    nlines = 0
    D = getHiscores ( dictfile )
    if len(D) == 0:
        return 1
    if nmax == None:
        nmax = 10
        if extended:
            nmax = 3
    if True:
        D = sortEntries ( D, nmax, nmax_analysis )
    for i,entry in enumerate ( D ):
        if i >= nmax: #  and extended:
            break
        wid = 0
        K, TL = entry['K'], entry['TL']
        if "walkerid" in entry:
            wid = entry['walkerid']
        particles = list ( entry["masses"].keys() )
        particles.sort ( key = lambda x: entry["masses"][x] )
        sparticles = ""
        for ip, p in enumerate ( particles ):
            if ip != 0:
                sparticles += ", "
            name = SParticleNames( False).asciiName(p)
            mass = entry["masses"][p]
            color = CYAN
            if "widths" in entry and p in entry["widths"]:
                w = entry["widths"][p]
                if 10**-19<=w<=10**-10:
                    color = ORANGE
            sparticles += f"{color}{name}{RESET}={mass:.1f}"
        timestamp = ""
        if "timestamp" in entry:
            timestamp = formatTimestamp ( entry["timestamp"] )
        if extended:
            step = entry["step"]
            print ( f"#{i+1:2d}({wid:3d}): K={GREEN}{K:6.3f}{RESET} TL={TL:6.3f}; {sparticles}" )
            print ( f"       `---: {entry['description']}" )
            print ( f"       `---:{timestamp}" )
            print ( f"       `---: step {step}" )
            print ( )
            nlines += 5
        else:
            sK = formatObject ( K, "6.3f" )
            print ( f"#{i+1:2d}({wid:3d}): K={GREEN}{sK}{RESET}; TL={formatObject(TL,'6.3f')}; {sparticles} {timestamp}" )
            nlines += 1
    return nlines

def runSlurmWalk() -> int:
    rundir = os.getcwd()
    path_to_clip_tools = f"{os.environ['HOME']}/git/smodels-utils/clip/"
    cmd = f"{path_to_clip_tools}/slurm_walk.py -R {rundir} -q"
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
    sK = formatObject ( K, '6.3f' )
    particles = list ( d["masses"].keys() )
    particles.sort ( key = lambda x: d["masses"][x] )
    sparticles = ""
    for ip, p in enumerate ( particles ):
        if ip != 0:
            sparticles += ", "
        name = SParticleNames( False).asciiName(p)
        mass = d["masses"][p]
        sparticles += f"{CYAN}{name}{RESET}={RED}{mass:.1f}{RESET}"
    print ( f"{RED}Truth:{RESET}    K={GREEN}{sK}{RESET}; TL={formatObject(TL,'6.3f')}; {sparticles}" )

def reportExceptionsLog():
    fname = "exceptions.log"
    # fname = "old_exceptions.log"
    if os.path.exists ( fname ):
        print ( f"{RED}ALERT: There is an {fname}!!!{RESET}" )

if __name__ == "__main__":
    import argparse
    argparser = argparse.ArgumentParser(
        description='summarize the data in e.g. hiscores_global.dict to have an overview' )
    argparser.add_argument ( '-H', '--hiscores', type=str,
        help="path to hiscores.dict file, if directory (e.g. 'all_hiscores/') then print best of each file in directory [./hiscores_global.dict]", 
        default="./hiscores_global.dict" )
    argparser.add_argument ( '-x', '--extended', action="store_true",
        help="extended info" )
    argparser.add_argument ( '-l', '--loop', action="store_true",
        help="loop" )
    argparser.add_argument ( '-N', '--nmax_analysis', type=int, default=None,
        help="print maximally this number of entries per analysis [None]" )
    argparser.add_argument ( '-n', '--nmax', type=int, default=None,
        help="print maximally this number of entries [None]" )
    args = argparser.parse_args()
    import colorama
    colorama.init()
    reportExceptionsLog()
    while True:
        nlines = runSlurmWalk()
        nlines += summarizeJobsForThisDir()
        nlines += summarizeHiscores ( args.hiscores, args.extended, args.nmax,
               args.nmax_analysis )
        if not args.loop:
            break
        time.sleep(10.)
        print( colorama.Cursor.UP()*(nlines+2) )
    reportExceptionsLog()
