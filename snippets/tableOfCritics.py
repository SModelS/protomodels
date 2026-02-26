#!/usr/bin/env python3

""" script that creates the latex table of the critics, for br(ll) where
the critic starts to veto """

def load():
    with open ( "dileptons.dict", "rt" ) as f:
        d = eval ( f.read() )
        f.close()
        return d

def discuss ( ll_ratio, values ):
    print ( r"\begin{tabular}{l|r|r}" )
    print ( r"\bf{Analysis} & \bf{ $r_\mathrm{obs}$ } & \bf{ $r_\mathrm{exp}$ }\\" )
    print ( r"\hline" )
    items = list ( values.items() )
    items = []
    for k,v in values.items():
        if type(v)==dict: #  and v["rexp"]>0.7:
            items.append ( (k,v) )
    items.sort ( key = lambda x: x[1]["rexp"], reverse=True )
    for k,v in items:
        dId = f"{v['anaid']}:{v['dataid']}"
        if v['dataid'] == None:
            dId = f"{v['anaid']}:UL"
        dId = dId.replace("_",r"\_")
        print ( rf"{dId:41s} & {v['robs']:.2f} & {v['rexp']:.2f} \\" )
    print ( r"\end{tabular}" )

def create ( d ):
    for ll_ratio, values in d.items():
        if values["allowed"]==True:
            continue
        discuss ( ll_ratio, values )
        return

def run():
    d = load()
    create ( d )

if __name__ == "__main__":
    run()
