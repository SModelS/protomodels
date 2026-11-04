#!/usr/bin/env python3
from pbase.constants import LSP, X2Z, X1W, X3Z, X2W

fname = "{os.environ['HOME']}/rundir.fudged1/history34.list" 
with open ( fname, "rt" ) as h:
    txt = h.read()
    txt = txt.replace("nan",'"nan"')
    if txt[-2] == ",":
        txt+="]\n"
    L = eval ( txt )

print ( len(L), "steps." )

for model in L:
    if not 1000006 in model["masses"]:
        continue
    step = model["step"]
    dm = model["masses"][1000006] - model["masses"][LSP]
    inCorridor = 150 < dm < 200
    if inCorridor:
        print ( "step", step , "dm", dm )

