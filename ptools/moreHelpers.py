#!/usr/bin/env python3

""" various helper functions that do not fit in any of the more
    specific modules """

def namesForSetsOfTopologies ( name ):
    """ some abbreviations for sets of topologies,
    e.g. electroweakino -> TChiWZ, TChiWH, .... 
    :param name: abbreviation
    :returns: comma separate list of topos, or original name if nothing found
    """
    shorts = { }
    shorts["electroweakino"]="TChiWZ,TChiWH,TChiWZoff,TChiZZ,TChiHH,TChiWW,TChiZ,TChiH,TChiWWoff,TChiZH"
    if name in shorts:
        return shorts[name]
    return name

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
