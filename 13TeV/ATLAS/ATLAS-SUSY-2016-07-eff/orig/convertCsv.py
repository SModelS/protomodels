#!/usr/bin/env python3

import glob

def convertUL ( oldn, newn ):
    with open(oldn,"rt") as f:
        lines=f.readlines()
    f=open(newn,"wt")
    f.write ( "# produced by converting %s\n" % oldn )
    for line in lines:
        line = line.strip()
        if len(line)==0:
            continue
        if line.startswith("#" ):
            continue
            # f.write ( line )
        #print ( ">%s<" % line )
        tokens = list(map(float,line.split(",")))
        # print ( tokens )
        # tokens = list(map(float,line.split(",")))
        l = "%s,%s,%s\n" % ( tokens[0],60.+tokens[1]*(tokens[0]-60.),tokens[2] ) 
        f.write ( l )
        
    f.close()

def convertExcl ( oldn, newn ):
    with open(oldn,"rt") as f:
        lines=f.readlines()
    f=open(newn,"wt")
    f.write ( "# produced by converting %s\n" % oldn )
    for line in lines:
        line = line.strip()
        if len(line)==0:
            continue
        if line.startswith("#" ):
            continue
            # f.write ( line )
        if "GEV" in line:
            continue
        #print ( ">%s<" % line )
        tokens = list(map(float,line.split(",")))
        # print ( tokens )
        # tokens = list(map(float,line.split(",")))
        l = "%s,%s\n" % ( tokens[0],60.+tokens[1]*(tokens[0]-60.) ) 
        #l = "%s,%s,%s\n" % ( tokens[0],60.+tokens[1]*(tokens[0]-60.),tokens[1] ) 
        f.write ( l )
        
    f.close()

if __name__ == "__main__":
    for idx in [ 4, 6, 8, 9 ]:
        convertUL ( "X-sectionU.L.&bestSR%d:Meff.csv" % idx, "UL_SR%d.csv" % idx )
        convertExcl ( "Exclusioncontour(exp.)%d:Meff.csv" % idx, "expSR%d.csv" % idx )
        convertExcl ( "Exclusioncontour(obs.)%d:Meff.csv" % idx, "obsSR%d.csv" % idx )
