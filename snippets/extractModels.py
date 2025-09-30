#!/usr/bin/env python3

""" take all models in Pmodels, agglomerate them into a simple list """


def extractModels( after : int = 10 ):
    import glob
    files = glob.glob ( "logs/walker_*.log" )
    models = []
    for fl in files:
        with open ( fl, "rt" ) as f:
            lines = f.readlines()
            stepnr = 0
            for line in lines:
                p2 = line.find(" begins")
                p1 = line.find("Step " )
                if p2 > 0 and p1 > 0:
                    tmp = line[p1+5:p2]
                    tmp.strip()
                    stepnr = int ( tmp )
                p3 = line.find("Protomodel:")
                if p3 > 0 and stepnr == after:
                    tmp = line[p3+11:]
                    m = eval ( tmp )
                    models.append ( m )
    from ptools.helpers import py_dump
    py_dump ( models, f"pmodels_{after}.list" )

if __name__ == "__main__":
    extractModels( 5 )
    extractModels( 10 )
    extractModels( 20 )
    extractModels( 50 )
