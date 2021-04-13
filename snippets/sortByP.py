#!/usr/bin/env python3

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
    for k in pkeys[:20]:
        print ( "p=%.3f, id=%s, txnames=%s" % ( k,pvalues[k], D[pvalues[k]]["txns"] ) )

run()
