#!/usr/bin/env python3

""" summarize the data in hiscores.dict to have an overview """

from os import PathLike
from ptools.sparticleNames import SParticleNames 
from colorama import Fore

def summarizeHiscores ( dictfile : PathLike = "hiscores.dict",
    extended : bool = False  ):
    """ summarize the content of the dict file 

    :param dictfile: path to dictionary file
    :param extended: extended output, add description timestamp
    """
    f=open( dictfile, "rt" )
    D=eval(f.read() )
    f.close()
    for i,entry in enumerate ( D ):
        if extended and i > 5:
            break
        K, Z, wid = entry['K'], entry['Z'], entry['walkerid']
        particles = entry["masses"].keys()
        sparticles = ""
        for ip, p in enumerate ( particles ):
            if ip != 0:
                sparticles += ", "
            name = SParticleNames( False).asciiName(p)
            mass = entry["masses"][p]
            sparticles += f"{name}={mass:.1f}"
        timestamp = entry["timestamp"]
        if extended:
            print ( f"#{i}({wid:3d}): K={Fore.GREEN}{K:.3f}{Fore.RESET}; {sparticles}" )
            print ( f"       `---: {entry['description']} {timestamp}" )
            print ( )
        else:
            print ( f"#{i}({wid:3d}): K={Fore.GREEN}{K:.3f}{Fore.RESET}; {sparticles} {timestamp}" )


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
