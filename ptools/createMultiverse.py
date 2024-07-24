#!/usr/bin/env python3

import subprocess, os

def isFudged ( f ):
    return abs(f-1.0)>1e-5

def create ( nmin = 1, nmax = 100, f = 1.0, overwrite = False ):
    """
    :param n: number of universes
    :param f: fudge factor
    :param overwrite: overwrite older versions of the files
    """
    directory = "dicts"
    if isFudged(f):
        directory = f"dictsf{int(f*100)}"
    if not os.path.exists ( directory ):
        try:
            os.mkdir ( directory )
        except FileExistsError as e:
            pass
    for i in range(nmin,nmax+1):
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
    argparser.add_argument ( '-N', help='highest id of fake universe to be created [100]',
            type=int, default=100 )
    argparser.add_argument ( '-n', help='lowest id of fake universe to be created [1]',
            type=int, default=1 )
    argparser.add_argument ( '-p', help='number of processes [5]',
            type=int, default=5 )
    argparser.add_argument ( '-f', '--fudge', help='fudge factor [1.0]',
            type=float, default=1.0 )
    argparser.add_argument ( '-o', '--overwrite', help='overwrite old files',
            action='store_true' )
    args = argparser.parse_args()
    nprocesses = 5
    dn = int ( ( args.N+1-args.n) / nprocesses )
    for p in range(nprocesses):
        pid = os.fork()
        if pid != 0:
            nmin = args.n + p * dn
            nmax = args.n + (p+1)*dn - 1
            print ( nmin, nmax )
            create( nmin, nmax, args.fudge, args.overwrite )
