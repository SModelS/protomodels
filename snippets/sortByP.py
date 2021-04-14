#!/usr/bin/env python3

""" tiny snippet to flash the most prominent p-values """

import subprocess

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
    f=open("pvalues.tex","wt")
    f.write( "\\documentclass{article}\n" )
    f.write( "\\usepackage{lscape}\n" )
    f.write( "\\usepackage{pdflscape}\n" )
    f.write( "\\begin{document}\n" )
    f.write( "\\thispagestyle{empty}\n" )
    f.write( "\\begin{landscape}\n" )
    f.write( "\\begin{table}\n" )
    f.write( "\\begin{tabular}{c|c|c|c}\n" )
    f.write( "$p$ & analysis & signal region & topologies \\\\\n" )
    f.write( "\\hline\n" )
    for k in pkeys[:]:
        if k>.05:
            break
        anaid = pvalues[k]
        txns = D[anaid]["txns"]
        skip = False
        if anaid == "ATLAS-CONF-2013-062:SL2m":
            txns = "txnames are all zero"
            skip = True
        print ( "p=%.3f, id=%s, txnames=%s" % ( k, anaid, txns ) )
        if skip:
            continue
        p1 = anaid.find(":")
        ana, ds = anaid[:p1], anaid[p1+1:]
        ptxns = txns
        dsmax = 40
        if len ( ds ) > dsmax:
            ds = ds[:dsmax-4]+"..."
        if len(ptxns)>dsmax:
            ptxns = txns[:dsmax-4]+"..."
        ds = ds.replace("_","\\_" )
        line = "%.3f & %s & %s & %s \\\\\n" % ( k, ana, ds, ptxns )
        f.write ( line )
    f.write ( "\\end{tabular}\n" )
    f.write ( "\\end{table}\n" )
    f.write( "\\end{landscape}\n" )
    f.write ( "\\end{document}\n" )
    f.close()
    cmd = "pdflatex pvalues.tex"
    print ( cmd )
    subprocess.getoutput ( cmd )
    cmd = "convert -alpha deactivate -trim -background white -antialias -density 600  pvalues.pdf pvalues.png"
    subprocess.getoutput ( cmd )

run()
