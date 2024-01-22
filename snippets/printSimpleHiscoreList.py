#!/usr/bin/env python3

""" summarize the data in hiscores.dict to have an overview """

from os import PathLike
from ptools.sparticleNames import SParticleNames 

def summarizeHiscores ( dictfile : PathLike = "hiscores.dict" ):
    """ summarize the content of the dict file """
    f=open( dictfile, "rt" )
    D=eval(f.read() )
    f.close()
    for i,entry in enumerate ( D ):
        K, Z = entry['K'], entry['Z']
        particles = entry["masses"].keys()
        sparticles = ""
        for ip, p in enumerate ( particles ):
            if ip != 0:
                sparticles += ", "
            name = SParticleNames( False).asciiName(p)
            mass = entry["masses"][p]
            sparticles += f"{name}={mass:.1f}"
        timestamp = entry["timestamp"]
        print ( f"#{i}: K={K:.2f}; {sparticles}; {timestamp}" )


if __name__ == "__main__":
    import argparse 
    argparser = argparse.ArgumentParser(
        description='summarize the data in hiscores.dict to have an overview' )
    argparser.add_argument ( '-H', '--hiscores', type=str,
        help="path to hiscores.dict file [./hiscores.dict]", 
        default="./hiscores.dict" )
    args = argparser.parse_args()
    summarizeHiscores ( args.hiscores )
