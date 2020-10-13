#!/usr/bin/env python3

def parse( inputfile = "../ptools/history.list" ):
    f=open( inputfile,"rt")
    txt=f.read()
    f.close()

    g=open("end.list","wt")
    ntxt=txt
    if not "]" in txt[-3:]:
        ntxt = txt[:-1] +"]"
    g.write ( ntxt+ "\n" )
    g.close()

    f=open("end.list","rt")
    L=eval(f.read())
    Kmax,maxstep=-90.,None
    for l in L:
        K,step=l["K"],l["step"]
        if K != None and K > Kmax:
            Kmax,maxstep=K,step
    print ( Kmax, maxstep )

if __name__ == "__main__":
    import argparse
    argparser = argparse.ArgumentParser( description="parse the history files" )
    argparser.add_argument ( '-i', '--inputfile', help="input history file",
                            type=str, default="../ptools/history.list" )
    args = argparser.parse_args()

    parse( args.inputfile )
