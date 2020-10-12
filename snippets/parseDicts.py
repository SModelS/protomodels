#!/usr/bin/env python3

""" a bit of code to help us with parsing the dicts of the runs """

import glob, subprocess, sys

def run( filename ):
    f = open ( filename, "rt" )
    D=eval(f.read())
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
    files = glob.glob ( f"{Dir}/rundir.*/states.dict" )
    for f in files:
        name = f.replace( Dir, "" ).replace("/states.dict","").replace("rundir.","")
        cmd = f"cp {f} {name}.dict"
        print ( cmd )
        subprocess.getoutput ( cmd )

def getBest( pattern="real" ):
    """ find out best scoring models """
    files = glob.glob ( f"{pattern}*.dict" )
    Ks = {}
    for f in files:
        with open ( f, "rt" ) as h:
            L = eval(h.read())
            K = L[0]["K"]
            Ks[K]=f
    keys = list ( Ks.keys() )
    keys.sort( reverse = True )
    print ()
    print ( f"Scores for {pattern}" )
    for k in keys[:20]:
        print ( k, Ks[k] )

if __name__ == "__main__":
    import argparse
    argparser = argparse.ArgumentParser(description="summarize all runs" )
    argparser.add_argument ( '-c', '--count',
               help='count the number of steps already taken', action='store_true' )
    argparser.add_argument ( '-p', '--produce',
               help='produce the stats', action='store_true' )
    args=argparser.parse_args()

    if args.count:
        count()
    if args.produce:
        fetch()
        produce()
        getBest()
        getBest("ewk")
        getBest("fake")
