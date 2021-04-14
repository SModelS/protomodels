#!/usr/bin/env python3

""" tiny snippet to flash the most prominent p-values """

def run():
    fname = "../ptools/dbnormal1.dict"
    f = open ( fname, "rt" )
    tmp = f.readlines()
    f.close()
    lines = []
    meta = eval ( tmp[0] )
    for line in tmp[1:]:
        if not line.startswith ( "#" ):
            lines.append ( line )
    D = eval ( "\n".join ( lines ) )
    pvalues = {}
    for anaid,values in D.items():
        if "orig_p" in values:
            pvalues[ values["orig_p"] ] = anaid
    pkeys = list ( pvalues.keys() )
    pkeys.sort()
    for k in pkeys[:]:
        if k>.05:
            break
        anaid = pvalues[k]
        txns = D[anaid]["txns"]
        if anaid == "ATLAS-CONF-2013-062:SL2m":
            txns = "txnames are all zero"
        print ( "p=%.3f, id=%s, txnames=%s" % ( k, anaid, txns ) )

run()
