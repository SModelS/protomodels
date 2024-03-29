#!/usr/bin/env python3

""" 
.. module:: bamFromPrettyNames
   :synopsis: construct a binary acceptance matrix from our pretty names

.. moduleauthor:: Wolfgang Waltenberger <wolfgang.waltenberger@gmail.com>

"""

from typing import Dict, List, Text
from smodels.experiment.databaseObj import Database
from smodels_utils.helper.databaseManipulations import filterSqrtsFromList, \
         filterSupersededFromList, filterCollaborationFromList
from smodels_utils.helper.various import getCollaboration

def getPrettyNames():
    dbpath = "../../smodels-database/"
    db = Database ( dbpath, progressbar=True )
    ers = db.getExpResults()
    ers = filterSqrtsFromList ( ers, 8 )
    ers = filterCollaborationFromList ( ers, "CMS" )
    ers = filterSupersededFromList ( ers )
    # ers = db.expResultList
    anaIds = { er.globalInfo.id for er in ers }
    prettyNames = { x.globalInfo.id : x.globalInfo.prettyName for x in ers if hasattr ( x.globalInfo, "prettyName" ) }
    return prettyNames

def add ( D : Dict, anaid1 : str, anaid2 : str ):
    D[anaid1].append ( anaid2 )
    D[anaid2].append ( anaid1 )

def exclusive ( string1 : str, string2 : str, pname1 : str, pname2 : str ) -> bool:
    """ check if string_i appears in one and only one of pname_i """
    if string1 in pname1 and string2 in pname2:
        return True
    if string1 in pname2 and string2 in pname1:
        return True
    return False

def inBoth ( string1 : str,  pname1 : str, pname2 : str ) -> bool:
    """ check if string1 appears in both pnames """
    if string1 in pname1 and string1 in pname2:
        return True
    return False

def ANotB ( A : Text, B : List[Text], pname1 : Text, pname2 : Text ):
    """ A is in one, none of the Bs are in the other """
    if A in pname1:
        for b in B:
            if b in pname2:
                return False
        return True
    if A in pname2:
        for b in B:
            if b in pname1:
                return False
        return True
    return False

def pprint ( ret : Dict ):
    """ pretty print the final matrix """
    for k,v in ret.items():
        print ( f"{k}: {v}" )

def save ( ret ):
    f = open ( "matrix.py", "wt" )
    f.write ( "    ww = {}\n" )
    for k,v in ret.items():
        vals = ", ".join ( [ f'"{x}"' for x in v ] )
        f.write ( f'    ww["{k}"]= [ {vals} ]\n' )
    f.write ( f'    return ww\n' )
    # f.write ( str(ret)+"\n" )
    f.close()

def synonyms ( pname : str ) -> str:
    """ replace synonymous terms """
    syns = {}
    pname = pname.replace(", EWK","")
    pname = pname.replace(" EWK","")
    syns["2 h(b b)"] = "0 l"
    syns["2 h(bb)"] = "0 l"
    syns["W h(b b)"] = "0 l"
    syns["W h(bb)"] = "0 l"
    syns["hadr."] = "0 l"
    syns["multi-jets"] = "0 l"
    syns["1-2 b-jets, MCT"] = "0 l"
    syns[">= 5 (1 b-)jets"] = "0 l"
    syns["2 b- or 2 c-jets"] = "0 l"
    syns["jets + boosted h(b b)"] = "0 l"
    syns["jets + boosted h(bb)"] = "0 l"
    syns["jets + top- and W-tag"] = "0 l"
    syns["jets + boosted Z"] = "0 l"
    syns["2 SFOS l"] = "2 l"
    syns["2 SS l"] = "2 l"
    syns["2 OS l"] = "2 l"
    syns["2 hadronic taus"] = "0 l"
    syns["2 b-jets"] = "0 l"
    syns[">= 2 c-jets"] = "0 l"
    syns["3 l"] = "multi-l"
    syns["2 SS or 3 l + jets"] = "multi-l"
    syns["jets + 2 SS or >= 3 l"] = "multi-l"
    syns["SFOS l"] = "2 l"
    syns["soft OS l"] = "soft l"
    for k,v in syns.items():
        if pname == k: # these need to match 100%
            pname = v
    matches = {}
    matches["2 OS l"] = "2 l"
    for k,v in matches.items():
        pname = pname.replace(k,v) # these will be replaced also within the pname
    return pname

def constructBAM():
    prettyNames = getPrettyNames()
    ret = { x: [] for x in prettyNames.keys() }
    ctUnknowns=0
    for anaid1, pname1 in prettyNames.items():
        pname1 = synonyms ( pname1 )
        coll1 = getCollaboration ( anaid1 )
        for anaid2, pname2 in prettyNames.items():
            if anaid2 == anaid1:
                continue
            if anaid2 < anaid1:
                continue
            coll2 = getCollaboration ( anaid2 )
            if coll1 != coll2:
                continue
            pname2 = synonyms ( pname2 )
            if exclusive ( "0 l", "1 l", pname1, pname2 ):
                add ( ret, anaid1, anaid2 )
                continue
            if exclusive ( "0 l", "2 l", pname1, pname2 ):
                add ( ret, anaid1, anaid2 )
                continue
            if exclusive ( "0 l", "2-3 l", pname1, pname2 ):
                add ( ret, anaid1, anaid2 )
                continue
            if exclusive ( "0 l", "multi-l", pname1, pname2 ):
                add ( ret, anaid1, anaid2 )
                continue
            if exclusive ( "1 l", "2 l", pname1, pname2 ):
                add ( ret, anaid1, anaid2 )
                continue
            if exclusive ( "1 l", "2-3 l", pname1, pname2 ):
                add ( ret, anaid1, anaid2 )
                continue
            if exclusive ( "1 l", "multi-l", pname1, pname2 ):
                add ( ret, anaid1, anaid2 )
                continue
            if exclusive ( "2 l", "multi-l", pname1, pname2 ):
                add ( ret, anaid1, anaid2 )
                continue
            if ANotB ( "HSCP", [ "HSCP", "LLP" ], pname1, pname2 ):
                add ( ret, anaid1, anaid2 )
                continue
            if ANotB ( "disappearing tracks", [ "disappearing tracks" ], pname1, pname2 ):
                add ( ret, anaid1, anaid2 )
                continue
            if ANotB ( "non-prompt jets", [ "non-prompt jets" ], pname1, pname2 ):
                add ( ret, anaid1, anaid2 )
                continue
            if ANotB ( "displaced vertices", [ "displaced vertices" ], pname1, pname2 ):
                add ( ret, anaid1, anaid2 )
                continue
            if ANotB ( "LLP", [ "HSCP", "LLP" ], pname1, pname2 ):
                add ( ret, anaid1, anaid2 )
                continue
            if ANotB ( "gamma", [ "gamma" ], pname1, pname2 ):
                add ( ret, anaid1, anaid2 )
                continue
            if ANotB ( "gamma gamma", [ "gamma gamma" ], pname1, pname2 ):
                add ( ret, anaid1, anaid2 )
                continue
            ## now the ones we know are correlated
            if exclusive ( "2 l", "2-3 l", pname1, pname2 ):
                continue
            if exclusive ( "multi-l", "2-3 l", pname1, pname2 ):
                continue
            checkForBoth = [ "0 l", "1 l", "2 l", "multi-l", "gamma gamma", 
                "disappearing tracks", "displaced vertices", "2-3 l",
                "HSCP", "gamma", "tau" ]
            iB = False
            for o in checkForBoth:
                if inBoth ( o, pname1, pname2 ):
                    iB = True
            if iB:
                continue
            if exclusive ( "HSCP", "LLP", pname1, pname2 ):
                continue
            if ANotB ( "soft l", [ "HSCP", "LLP" ], pname1, pname2 ):
                ## combination is correlated with everything but HSCP and LLP
                continue
            if ANotB ( "comb", [ "HSCP", "LLP" ], pname1, pname2 ):
                ## combination is correlated with everything but HSCP and LLP
                continue
            if ANotB ( "tau", [ "tau" ], pname1, pname2 ):
                ## combination is correlated with everything but HSCP and LLP
                continue
            ctUnknowns+=1
            print ( f"{ctUnknowns}: what about this: {pname1} <-> {pname2}" )
            # print ( f"    : {anaid1}, {anaid2}" )
    # pprint ( ret )
    save ( ret )

if __name__ == "__main__":
    constructBAM()
