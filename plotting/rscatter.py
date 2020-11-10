#!/usr/bin/env python3

import sys
sys.path.append(".")

from roughviz.charts.scatter import Scatter
import glob, subprocess, os

def getData ():
    """ collect the data to scatter plot """
    x,y=[],[]
    Dir="../protomodels/data/"
    Dir = os.getcwd()
    Dir = Dir.replace("/protomodels/plotting","")
    Dir = Dir.replace("/py-roughviz","")
    Dir = Dir + "/protomodels/data"
    files = glob.glob ( f"{Dir}/real*dict" )
    for f in files:
        h = open ( f, "rt" )
        D = eval ( h.read() )
        h.close()
        for i in D:
            D0m = i["masses"]
            if not 1000006 in D0m:
                continue
            x.append ( D0m[1000006] )
            y.append ( D0m[1000022] )
    print ( f"{len(files)} files, {len(x)} points." )
    return x,y

x,y=  getData()
scatter = Scatter(data={"x": x, "y": y}, radius=30, margin= {"top": 60, "right": 20, "bottom": 70, "left": 170}, interactive=False )
scatter.set_title("Distribution of masses of runs", fontsize=3)
scatter.set_xlabel("Xt [GeV]", fontsize=2)
scatter.set_ylabel("XZ [GeV]", fontsize=2 )

scatter.show()

with open ( "bla.html", "wt" ) as f:
    f.write ( scatter.output+"\n" )
with open ( "blu.html", "wt" ) as f:
    more = scatter.output
    print ( "tp", type(more) )
    f.write ( more )

# print ( scatter.output )
