#!/usr/bin/env python3

""" take all models in Pmodels, agglomerate them into a simple list """

def extractModels( after : int = 10 ):
    import glob, sys
    from datetime import datetime
    files = glob.glob ( "logs/walker_*.log" )
    models = {}
    startTime, endTime = None, None
    for fl in files:
        with open ( fl, "rt" ) as f:
            lines = f.readlines()
            stepnr = 0
            weDone=False
            for line in lines:
                if "It is" in line:
                    p1 = line.find("-")
                    p2 = line.find("]")
                    fmt = "%H:%M:%S"
                    startTime = datetime.strptime(line[p1+1:p2], fmt)
                if weDone == True:
                    break
                p2 = line.find(" begins")
                p1 = line.find("Step " )
                if p2 > 0 and p1 > 0:
                    tmp = line[p1+5:p2]
                    tmp.strip()
                    stepnr = int ( tmp )
                p3 = line.find("Protomodel:")
                if p3 > 0: #  and stepnr == after:
                    tmp = line[p3+11:]
                    p1 = line.find("-")
                    p2 = line.find("]")
                    fmt = "%H:%M:%S"
                    endTime = datetime.strptime(line[p1+1:p2], fmt)
                    m = eval ( tmp )
                    m["dt"] = (endTime - startTime).seconds
                    models[fl]=m
                    
                if stepnr == after:
                    weDone = True
                    break
    from ptools.helpers import py_dump
    py_dump ( list(models.values()), f"pmodels_{after}.list" )

if __name__ == "__main__":
    extractModels( 3 )
    extractModels( 5 )
    extractModels( 10 )
    extractModels( 20 )
    extractModels( 50 )
