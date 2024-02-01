#!/usr/bin/env python3

""" script to create table of runs, with test statistics and particle contents """

import glob
from smodels_utils.helper.sparticleNames import SParticleNames

def getDicts():
    basedir = "/scratch-cbe/users/wolfgan.waltenberger/"
    Tp = "rundir.real"
    files = glob.glob ( basedir + Tp + "*/states.dict" )
    Dicts = {}
    for File in files:
        nr = File.replace( basedir, "" ).replace ( "/states.dict", "" )
        nr = nr.replace ( Tp, "" )
        nr = int(nr)
        with open ( File, "rt" ) as f:
            lines=f.readlines()
            f.close()
            # print ( "l", lines[0][1:-2] )
            Dicts[nr]=eval(lines[0][1:-2])
            for x in [ "ssmultipliers", "decays", "timestamp", "dbver", "Z", "codever", \
                       "step" ]:
                Dicts[nr].pop( x )
    return Dicts

def writeTex ( D ):
    with open ( "./realruns.tex", "wt" ) as f:
        f.write ( "\\begin{tabular}{r|r|l}\n" )
        f.write ( "run & \\K & Particles \\\\\n" )
        f.write ( "\\hline\n" )
        keys = list ( D.keys() )
        keys.sort()
        namer = SParticleNames ( susy=False )
        for key in keys:
            v = D[key]
            pids = list ( v["masses"].keys() )
            def sorter ( x ):  ## we sort such that stop is first and the lbp is last
                ## which is not the real order but a close bet
                if x == 1000006:
                    x-=100000
                if x == 1000022:
                    x+= 1000000
                return x
            pids.sort( key= sorter )
            # particles = ", ".join ( map ( namer.texName, pids ) )
            particles = namer.texName ( pids, addDollars=True )
            # print ( "p", particles )
            f.write ( "\\#%s & %.2f & %s \\\\\n" % \
                       ( key, v["K"], particles ) )
            ### \#1 & 7.64 & $X_{t}, X_{c}, X^1_Z$ \\
        f.write ( "\\end{tabular}\n" )
        f.close()

def main():
    D = getDicts()
    writeTex ( D )

main()
