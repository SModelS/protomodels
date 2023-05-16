#!/usr/bin/env python3

""" Code that decides which analyses can be combined and which cannot """

from smodels.tools import runtime
runtime._cap_likelihoods = True
from smodels.theory.theoryPrediction import TheoryPrediction
import fnmatch

moreComments = { ## collect a few more comments on analyses
    "CMS-SUS-18-002": "lepton veto",
    "ATLAS-SUSY-2018-31": "lepton veto",
    "CMS-SUS-17-006": "lepton veto, no photon veto",
    "CMS-SUS-19-006": "veto on photon pt > 100",
    "CMS-SUS-16-043": "one lepton",
    "CMS-SUS-16-046": "photon pt > 25",
    "CMS-SUS-18-002": "photon pt > 100",
#    "CMS-SUS-17-005": "one lepton",
    "CMS-SUS-17-004": "combination, many channels",
    "CMS-SUS-17-010": "2 OS leptons",
    "CMS-SUS-17-009": "2 OS leptons",
    "CMS-SUS-16-039": "3+ leptons",
    "CMS-SUS-16-033": "photon veto",
    "CMS-SUS-16-009": "all hadronic, no photon veto",
    "CMS-SUS-16-041": "3+ leptons",
    "ATLAS-SUSY-2017-01": "all channels",
    "ATLAS-SUSY-2015-01": "lepton veto",
    "ATLAS-SUSY-2016-26": "lepton veto",
    "ATLAS-SUSY-2016-27": "no lepton vetos",
    "ATLAS-SUSY-2016-28": "0 or 1 leptons",
    "ATLAS-SUSY-2016-19": "different channels",
    "ATLAS-SUSY-2018-16": "2 low pt OS leptons",
    "ATLAS-SUSY-2017-02": "lepton veto, no photon veto",
    "ATLAS-SUSY-2016-07": "lepton veto, no photon veto",
    "ATLAS-SUSY-2016-15": "lepton veto, no photon veto",
    "ATLAS-SUSY-2018-04": "b-jet veto, no lepton veto",
    "ATLAS-SUSY-2017-03": "2 or 3 leptons",
}

def getExperimentName ( globI ):
    """ returns name of experiment of exp result """
    if "CMS" in globI.id:
        return "CMS"
    if "ATLAS" in globI.id:
        return "ATLAS"
    return "???"

def getInfoFromAnaId ( anaId, results ):
    """ given the analysis id as a string, get the globalInfo object """
    ret = None
    for i in results:
        if i.globalInfo.id == anaId:
            ret = i.globalInfo
            break
    if type(ret) == type(None):
        print ( "[analysisCombiner] supplied an analysis id string (%s), but couldnt find such an analysis in the database!" % anaId )
    return ret

def canCombine ( predA, predB, results = None ):
    """ can we combine predA and predB? predA and predB can be
        individual predictions, or lists of predictions.
    :param predA: lists of predictions, or single predictions (given either
                  as a TheoryPrediction or a GlobalInfo object)
    :param results: optionally supply a list of results. In this case, predA
                    and predB can be the analysis ids as strings
    """
    if type(predA) == str:
        predA = getInfoFromAnaId ( predA, results )
    if type(predB) == str:
        predB = getInfoFromAnaId ( predB, results )
    if type(predA) == type(None) or type(predB) == type(None):
        return False
    if type(predA) == list:
        for pA in predA:
            ret = canCombine ( pA, predB )
            if ret == False:
                return False
        return True
    if type(predB) == list:
        for pB in predB:
            ret = canCombine ( predA, pB )
            if ret == False:
                return False
        return True
    elA, elB = None, None
    if type(predA)==TheoryPrediction:
        elA = predA.elements
        predA = predA.expResult.globalInfo
    if type(predB)==TheoryPrediction:
        elB = predB.elements
        predB = predB.expResult.globalInfo
    return canCombineUsingMatrix ( predA, predB, elA, elB )

def hasOverlap ( elA, elB, globA = None, globB = None ):
    """ is there an overlap in the elements in A and in B? """
    if elA is None or elB is None:
        return False
    for eA in elA:
        for eB in elB:
            # print ( "el %s is el %s? %s! %s!" % ( eA, eB, eA.__cmp__ ( eB ), eA == eB ) )
            if eA == eB: ## an element of elA is in elB
                return True
    return False

def canCombineUsingMatrix ( globA, globB, elA, elB ):
    """ method that defines what we allow to combine, using combinationsmatrix.py
         """
    if globA.sqrts != globB.sqrts:
        return True
    if getExperimentName(globA) != getExperimentName(globB):
        return True
    if hasOverlap ( elA, elB, globA, globB ):
        ## overlap in the constraints? then for sure a no!
        return False
    anaidA = globA.id
    anaidB = globB.id
    from combinationsmatrix import getMatrix
    allowCombination = getMatrix()
    if anaidA in allowCombination.keys():
        for i in allowCombination[anaidA]:
            if len ( fnmatch.filter ( [anaidB ], i ) ) == 1:
                return True
    if anaidB in allowCombination.keys():
        for i in allowCombination[anaidB]:
            if len ( fnmatch.filter ( [anaidA ], i ) ) == 1:
                return True
        if anaidA in allowCombination[anaidB]:
            return True
    return False

def getSummary():
    from smodels.experiment.databaseObj import Database
    # dbpath = "official"
    dbpath = "<rundir>/database.pcl"
    print ( "[analysisCombiner] checking %s" % dbpath )
    db = Database ( dbpath )
    results = db.getExpResults()
    ana1 = "CMS-SUS-16-042"
    ana2 = "CMS-SUS-16-033"
    canC = canCombine ( ana1, ana2, results )
    print ( "[analysisCombiner] can combine %s with %s: %s" % ( ana1, ana2, str(canC) ) )
    ctr,combinable=0,0
    for x,e in enumerate(results):
        for y,f in enumerate(results):
            if y <= x:
                continue
            ctr += 1
            isUn = canCombine ( e.globalInfo, f.globalInfo )
            combinable+=isUn
    print ( "[analysisCombiner] can combine %d/%d pairs of results" % ( combinable, ctr ) )

def checkOneAnalysis():
    import argparse
    green, reset = "", ""
    try:
        import colorama
        green, reset = colorama.Fore.GREEN, colorama.Fore.RESET
    except Exception as e:
        pass
    #import IPython
    argparser = argparse.ArgumentParser(
            description='print the correlations of one specific analysis')
    argparser.add_argument ( '-d', '--dbpath',
            help='specify path to database [<rundir>/database.pcl]',
            type=str, default="<rundir>/database.pcl" )
    argparser.add_argument ( '-a', '--analysis',
            help='print for <analysis>',
            type=str, default="CMS-SUS-19-006" )
    args = argparser.parse_args()
    from smodels.experiment.databaseObj import Database
    print ( "[analysisCombiner] checking %s" % args.dbpath )
    db = Database ( args.dbpath )
    results = db.getExpResults()
    info = getInfoFromAnaId ( args.analysis, results )
    sqrts = info.sqrts
    collaboration = getExperimentName ( info )
    prettyName = info.prettyName
    if args.analysis in moreComments:
        prettyName += " (%s)" % moreComments[args.analysis]
    # IPython.embed()
    print ( "correlations for %s: %s" % (  args.analysis, prettyName ) )
    combs, nocombs = set(), set()
    pnames = {}
    for er in results:
        if er.globalInfo.sqrts != sqrts:
            continue
        if getExperimentName (er.globalInfo ) != collaboration:
            continue
        Id = er.globalInfo.id
        Id = Id.replace("-eff","").replace("-agg","")
        if Id == "CMS-SUS-19-006-2":
            Id = "CMS-SUS-19-006"
        if Id == args.analysis:
            continue
        pname = er.globalInfo.prettyName
        if Id in moreComments:
            pname += " (%s)" % moreComments[Id]
        pnames[Id]=pname
        cc = canCombine ( info, er.globalInfo )
        if cc:
            combs.add ( Id )
        else:
            nocombs.add ( Id )
    print ( f"{green}can combine with:{reset} " )
    for Id in combs:
        pname = pnames[Id]
        print ( " `- %s: %s" % ( Id, pname ) )
    print ( f"{green}cannot combine with:{reset} " )
    for Id in nocombs:
        pname = pnames[Id]
        print ( " `- %s: %s" % ( Id, pname ) )

if __name__ == "__main__":
    checkOneAnalysis()
    # getSummary()
