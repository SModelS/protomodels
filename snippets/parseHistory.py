#!/usr/bin/env python3

f=open("../ptools/history.list","rt")
txt=f.read()
f.close()

g=open("end.list","wt")
g.write ( txt[:-1] )
g.write ( "]\n" )
g.close()

f=open("end.list","rt")
L=eval(f.read())
Kmax,maxstep=-90.,None
for l in L:
    K,step=l["K"],l["step"]
    if K != None and K > Kmax:
        Kmax,maxstep=K,step
print ( Kmax, maxstep )
