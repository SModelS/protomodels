#!/usr/bin/env python3

""" check for first missing pic in movie making """

import glob, sys

files = glob.glob ( "smooth*.png" )
nrs = set()
for f in files:
    nr = f.replace("smooth","").replace(".png","")
    nrs.add ( int(nr) )

# print ( nrs )
for nr in range(50000):
    if not nr in nrs:
        print ( "missing", nr, "step", nr/40 )
        sys.exit()
