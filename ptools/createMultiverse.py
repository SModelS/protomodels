#!/usr/bin/env python3

import subprocess, os

def isFudged ( f ):
    return abs(f-1.0)>1e-5

def create ( n = 100, f = 1.0, overwrite = False ):
    """
    :param n: number of universes
    :param f: fudge factor
    :param overwrite: overwrite older versions of the files
    """
    directory = "dicts"
    if isFudged(f):
        directory = f"dictsf{int(f*100)}"
    if not os.path.exists ( directory ):
        os.mkdir ( directory )
    for i in range(1,n+1):
        out = f"{i:03d}"
        if os.path.exists ( f"{directory}/{out}.dict" ) and not overwrite:
            continue
        if isFudged(f):
            out += f"f{int(f*100)}"
        cmd = f"./ptools/expResModifier.py -C -d official -o none -s {out}"
        if isFudged(f):
            cmd += f" -f {f}"
        print ( f"{i}: {cmd}" )
        o = subprocess.getoutput ( cmd )
        if len(o)==0:
            o = "done"
        print ( f"   `-: {o}" )

        cmd = f"mv {out}.dict {directory}"
        print ( f"{i}: {cmd}" )
        o = subprocess.getoutput ( cmd )
        if len(o)==0:
            o = "done"
        print ( f"   `-: {o}" )

if __name__ == "__main__":
    import argparse
    argparser = argparse.ArgumentParser( description="tool to create a consistent set of fake universes at once" )
    argparser.add_argument ( '-n', help='number of fake universes [100]',
            type=int, default=100 )
    argparser.add_argument ( '-f', '--fudge', help='fudge factor [1.0]',
            type=float, default=1.0 )
    argparser.add_argument ( '-o', '--overwrite', help='overwrite old files',
            action='store_true' )
    args = argparser.parse_args()
    create( args.n, args.fudge, args.overwrite )
