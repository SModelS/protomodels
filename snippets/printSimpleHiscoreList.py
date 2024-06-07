#!/usr/bin/env python3

""" summarize the data in hiscores.dict to have an overview """

import os
from os import PathLike
from ptools.sparticleNames import SParticleNames 
from colorama import Fore as ansi

def summarizeHiscores ( dictfile : PathLike = "hiscores.dict",
    extended : bool = False  ):
    """ summarize the content of the dict file 

    :param dictfile: path to dictionary file
    :param extended: extended output, add description timestamp
    """
    if not os.path.exists ( dictfile ):
        print ( f"[printSimpleHiscoreList] {dictfile} does not exist" )
        return
    f=open( dictfile, "rt" )
    D=eval(f.read() )
    f.close()
    for i,entry in enumerate ( D ):
        if extended and i > 2:
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
            sparticles += f"{name}={mass:.1f}"
        timestamp = ""
        if "timestamp" in entry:
            timestamp = entry["timestamp"]
            r1 = timestamp.find(" ")
            r2 = timestamp.rfind(" ")
            timestamp = timestamp[r1:r2]
        if extended:
            print ( f"#{i}({wid:3d}): K={ansi.GREEN}{K:.3f}{ansi.RESET} TL={TL:.3f}; {sparticles}" )
            print ( f"       `---: {entry['description']}" )
            print ( f"       `---: {timestamp}" )
            print ( )
        else:
            print ( f"#{i}({wid:3d}): K={ansi.GREEN}{K:.3f}{ansi.RESET}; TL={TL:.3f}; {sparticles} {timestamp}" )


if __name__ == "__main__":
    import argparse 
    argparser = argparse.ArgumentParser(
        description='summarize the data in hiscores.dict to have an overview' )
    argparser.add_argument ( '-H', '--hiscores', type=str,
        help="path to hiscores.dict file [./hiscores.dict]", 
        default="./hiscores.dict" )
    argparser.add_argument ( '-x', '--extended', action="store_true",
        help="extended info" )
    args = argparser.parse_args()
    summarizeHiscores ( args.hiscores, args.extended )
