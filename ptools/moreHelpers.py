#!/usr/bin/env python3

""" various helper functions that do not fit in any of the more
    specific modules """

from typing import Union, Text, List, Tuple

def namesForSetsOfPids ( names: List[str|int] ) -> Tuple[int]:
    """ short names for various sets of pids. used as abbreviations
    for forbiddenparticles

    :param names: e.g. [ "stops","sbottoms",2000037"
    :returns: list of pids as ints e.g. "1000006,1000005"
    """
    dic = { "stops": [ 1000006, 2000006 ], 
            "sbottoms": [ 1000005, 2000005 ],
            "squarks": [ 1000001, 2000001, 1000002, 2000002, 1000003, 2000003, \
                         1000004, 2000004 ],
            "gluinos": [ 1000021 ],
            "charginos": [ 1000024, 1000037 ],
            "neutralinos": [ 1000023, 1000025, 1000036 ],
            "sneutrinos": [ 1000012, 1000014, 1000016 ],
            "sleptons": [ 1000011, 1000013, 1000015 ],
            "electroweakinos": [ 1000024, 1000037, 1000023, 1000025, 1000036 ]
    }
    ret = []
    for n in names:
        if n in dic.keys():
            for pid in dic[n]:
                ret.append ( pid )
        else:
            try:
                n=int(n)
            except ValueError as e:
                print ( f"[moreHelpers.namesForSetsOfPids] do not know {n}, will ignore" )
                continue
            ret.append ( n )
    return ret

def namesForSetsOfTopologies ( name : Union[Text,List,Tuple,None] ) \
        -> Tuple[Text,Union[Text,None]]:
    """ some abbreviations for sets of topologies,
    e.g. electroweakino -> TChiWZ, TChiWH, .... 
    :param name: abbreviation
    :returns: string with comma separated list of topos, and description
              if not an abbreviation, returns originalname, None
    """
    if name is None:
        return "all", None
    if "," in name:
        name = name.split(",")
    if type(name) in [ list, tuple ]:
        topos, descriptions = [], []
        for n in name:
            topo,descr = namesForSetsOfTopologies ( n )
            topos.append ( topo )
            if descr == None:
                descr = topo
            descriptions.append ( descr )
        return ( ",".join(topos), "+".join(descriptions) )

    shorts, description = { }, {}
    shorts["gauginos"]="TChiWZ,TChiWH,TChiZZ,TChiHH,TChiWW,TChiZH,TChiZ,TChiH"
    shorts["gauginos_offshell"]=shorts["gauginos"]+",TChiWZoff,TChiWWoff"
    shorts["electroweakinos"]='TChiH,TChiChipmSlepStau,TChipChimSlepSnu,TChiHH,TChiWW,TChiWZ,TChipChimSlepSlep,TChipChimStauSnu,TChiZH,TChiWH,TChipChimgg,TChiZZ,TChiChipmStauStau,TChiChipmSlepSlep,TChiChipmSlepL,TChiChipmStauL'
    shorts["electroweakinos_offshell"]=shorts["electroweakinos"]+",TChiWZoff,TChiWWoff"
    shorts["stops"]="T2tt,T2ttoff,T2bbffff,T2bbWW,T2bbWWoff,T6bbWW,T6bbWWoff"
    shorts["sbottoms"]="T2bb,T6ttWW,T6ttWWoff"
    shorts["colored"]="T1,T2,TGQ,T3GQ,T5GQ,TGQqtt,TGQbtq,TGQbbq,T1bbbb,T1tttt,T1bbbboff,T1ttttoff,T1btbt,T6WW"
    description["gauginos"]="ewkinos + onshell gauge bosons"
    description["electroweakinos"]="TChi* (only on-shell)"
    description["electroweakinos_offshell"]="TChi*"
    description["stops"]="stops, on- and off-shell"
    description["sbottoms"]="sbottoms"
    description["colored"]="light squarks and gluinos"
    if name in shorts:
        d = None
        if name in description:
            d= description[name]
        return shorts[name], d
    try:
        import Levenshtein
        rmin, kmin = 1., None
        for k,v in shorts.items():
            r = Levenshtein.distance(name,k)/(len(k)+len(name))
            if r < rmin:
                rmin, kmin = r, k
        if rmin < 0.2:
            print ( f"I assume you meant {kmin}, not {name}?" )
            d = None
            if kmin in description:
                d= description[kmin]
            return shorts[kmin], d
        #if rmin < 0.3:
        #    print ( f"could not find {name}, did you mean {kmin}?" )
        #    return name, None
    except ImportError as e:
        pass
    return name, None

def findLargestExcess ( db ):
    """ find the largest excess in any efficiency map type result
        in the given database
    :param db: a SModelS database object
    :returns: the dataset object
    """
    results = db.getExpResults ( dataTypes = [ "efficiencyMap" ] )
    excesses = {}
    for expRes in results:
        datasets = expRes.datasets
        for dataset in datasets:
            nobs = dataset.dataInfo.observedN
            nbg = dataset.dataInfo.expectedBG
            bgErr = dataset.dataInfo.bgError
            S = 0.
            toterr = math.sqrt ( bgErr**2 + nbg )
            if toterr > 0.:
                S = ( nobs - nbg ) / toterr
            if S < 1.:
                continue
            if not S in excesses:
                excesses[S]=[]
            excesses[S].append ( dataset )

    def pprint ( excesses ):
        keys = list ( excesses.keys() )
        keys.sort()
        for k in keys[-5:]:
            ds = excesses[k]
            if len(ds)!=1:
                print ( "error cannot handle" )
                continue
            ds = ds[0]
            obsN = ds.dataInfo.observedN
            eBG = ds.dataInfo.expectedBG
            print ( "Z=%.2f: %15s, %s: %d/%.2f" % \
                    ( k, ds.globalInfo.id, str(ds.dataInfo.dataId), obsN, eBG ) )

    pprint ( excesses )
    print ( "[helpers.findLargestExcess] found %d eff maps" % len(results) )
    return excesses

if __name__ == "__main__":
    """"
    In = "electroweakinos"
    Out = namesForSetsOfTopologies ( In )
    """
    In = [ "stops","sbottoms",2000037 ]
    Out = namesForSetsOfPids ( In )
    print ( In, "->", Out )
