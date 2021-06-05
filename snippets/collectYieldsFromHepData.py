#!/usr/bin/env python3

""" from csv tables in hep data, create dictionary of nobs, nbg, bgerr """

import glob, sys, numpy
sys.path.insert(0,"../")
from ptools.helpers import computeP

def parseFile ( fname ):
    f=open( fname, "rt" )
    lines = f.readlines()
    f.close()
    isIn = "before"
    # data = { "obs": [], "bg": [] }
    data={}
    for line in lines:
        if "Total BG" in line:
            isIn = "bg"
            continue
        if "Observed [Events]" in line:
            isIn = "obs"
            continue
        if line.startswith ( "#" ):
            continue
        if isIn == "before":
            continue
        line = line.strip()
        tokens = line.split(",")
        if len(tokens)<2:
            continue
        tokens = list ( map(float, tokens) )
        tokens[0]=int(tokens[0])
        binnr = tokens[0]
        tokens = tuple ( tokens )

        if not binnr in data:
            data[binnr]={}
        if isIn == "obs":
            data[binnr]["origN"]=tokens[1]
        if isIn == "bg":
            data[binnr]["expectedBG"]=tokens[1]
            data[binnr]["bgError"]=max( abs(tokens[2]), abs(tokens[3]) )

    for binnr, values in data.items():
        if binnr % 10 == 0:
            print ( ".", flush=True, end="" )
        p = computeP ( values["origN"], values["expectedBG"], values["bgError"] )
        values["orig_p"] = p
    return data

def read():
    anaid = "CMS-SUS-19-006"
    f=open(anaid+".py","rt" )
    D=eval(f.read())
    ps= []
    for anaid,values in D.items():
        ps.append( values["orig_p"] )
    print ( numpy.mean(ps ) )
    return D

def main():
    dirname = "/home/walten/Downloads/HEPData-ins1749379-v1-csv/"
    # fname = f"{dirname}Yieldsandpost-fitBG,Njets2-3.csv"
    files = glob.glob ( f"{dirname}/Yieldsandpost-fit*.csv" )
    alld = {}
    for fname in files:
        data = parseFile ( fname )
        alld.update ( data )
    entries = {}
    anaid = "CMS-SUS-19-006"
    for binnr,values in alld.items():
        label = anaid+":SR"+str(binnr)
        values["lumi"]=137.
        values["fudge"]=1.
        entries[label]=values
    f=open ( anaid+".py", "wt" )
    f.write ( "{" )
    for label,values in entries.items():
        f.write ( "'%s': %s,\n" % ( label, str(values) ) )
    f.write ( "}\n" )
    f.close()
    print ( f"{len(alld)} entries" )

read()
#main()
