#!/usr/bin/env python3

""" script that draws the dileptonic stuff """

def load():
    with open ( "dileptons.dict", "rt" ) as f:
        d = eval ( f.read() )
        f.close()
        return d

def convert( d : dict ):
    ret = { "ll_ratio": [] }
    for k in d[0.]:
        ret[k]=[]
    for i,(ll_ratio, values) in enumerate(d.items()):
        if ll_ratio == "meta":
            continue
        ret["ll_ratio"].append ( ll_ratio )
        for k,v in values.items():
            if not k in ret:
                ret[k]=[ float("nan" ) ]*i
            ret[k].append ( v )
    return ret

def draw( d : dict ):
    from matplotlib import pyplot as plt
    v = convert ( d )
    print  ( v["ll_ratio"] )
    plt.plot ( v["ll_ratio"], v["K"], label="K" )
    plt.plot ( v["ll_ratio"], v["TL"], label="TL" )
    plt.plot ( v["ll_ratio"], v["robsmax"], label="robsmax" )
    plt.legend ()
    plt.savefig ( "dileptons.png" )
    from smodels_utils.plotting.mpkitty import timg
    timg ( "dileptons.png" )


def run():
    d = load()
    draw ( d )

if __name__ == "__main__":
    run()
