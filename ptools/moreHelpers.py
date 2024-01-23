#!/usr/bin/env python3

""" various helper functions that do not fit in any of the more
    specific modules """

from typing import Union, Text, List, Tuple

def namesForSetsOfTopologies ( name : Union[Text,List,Tuple] ) \
        -> Tuple[Text,Union[Text,None]]:
    """ some abbreviations for sets of topologies,
    e.g. electroweakino -> TChiWZ, TChiWH, .... 
    :param name: abbreviation
    :returns: string with comma separated list of topos, and description
              if not an abbreviation, returns originalname, None
    """
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
    shorts["electroweakinos_offshell"]="TChiWZ,TChiWH,TChiWZoff,TChiZZ,TChiHH,TChiWW,TChiZ,TChiH,TChiWWoff,TChiZH"
    shorts["electroweakinos"]="TChiWZ,TChiWH,TChiZZ,TChiHH,TChiWW,TChiZH,TChiZ,TChiH"
    shorts["stops"]="T2tt,T2ttoff,T2bbffff,T2bbWW,T2bbWWoff,T6bbWW,T6bbWWoff"
    shorts["sbottoms"]="T2bb,T6ttWW,T6ttWWoff"
    description["electroweakinos"]="ewkinos + onshell gauge bosons"
    description["stops"]="stops, on- and off-shell"
    description["sbottoms"]="sbottoms"
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
    In = "T1"
    In = "electroweakinos"
    Out = namesForSetsOfTopologies ( In )
    print ( In, "->", Out )
