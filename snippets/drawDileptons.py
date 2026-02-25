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

def getExcludedRegions ( mask : dict ) -> list:
    import itertools

    false_ranges = []
    for k, g in itertools.groupby(mask.items(), key=lambda x: x[1]):
        if not k:  # False group
            xs = [item[0] for item in g]
            false_ranges.append((min(xs), max(xs)))
    return false_ranges

def draw( d : dict ):
    from matplotlib import pyplot as plt
    v = convert ( d )
    # print  ( v["ll_ratio"] )
    fig, ax1 = plt.subplots()
    lK, = plt.plot ( v["ll_ratio"], v["K"], label="K" )
    lTL, = plt.plot ( v["ll_ratio"], v["TL"], label="TL" )
    ax1.set_ylabel ( "K,TL")
    ax2 = ax1.twinx()
    lrmax, = ax2.plot ( v["ll_ratio"], v["robsmax"], label="max($r_{obs}$)", c="green" )
    ax2.set_ylabel ( "r" )
    ax2.set_ylim(0, 13)
    ax1.set_xlabel ( r"br($X_Z^2 \rightarrow l l X_Z^1$)"  )
    mask_ul, mask_llhd = {}, {}
    for ll_ratio, values in d.items():
        if ll_ratio == "meta": 
            continue
        mask_ul[ll_ratio] = values["allowed"]
        mask_llhd[ll_ratio] = values["llhd_allowed"]

    excluded_ul = getExcludedRegions ( mask_ul ) 
    excluded_llhd = getExcludedRegions ( mask_llhd ) 
    for x_start, x_end in excluded_ul:
        lexcl_ul=ax1.axvspan(x_start, x_end, color='gray', alpha=0.2,label="excluded by fast critic")
    for x_start, x_end in excluded_llhd:
        lexcl_llhd=ax1.axvspan(x_start, x_end, color='gray', alpha=0.5,label="excluded by fast+slow critic")
    lines = [ lK, lTL, lrmax,lexcl_ul, lexcl_llhd ]
    labels = [l.get_label() for l in lines]
    ax1.legend(lines, labels, loc="best")
    plt.savefig ( "dileptons.png" )
    from smodels_utils.plotting.mpkitty import timg
    timg ( "dileptons.png" )


def run():
    d = load()
    draw ( d )

if __name__ == "__main__":
    run()
