#!/usr/bin/env python3

""" Code that decides which analyses can be combined and which cannot """

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

def canCombine ( predA, predB, strategy="aggressive", results = None ):
    """ can we combine predA and predB? predA and predB can be
        individual predictions, or lists of predictions.
    :param predA: lists of predictions, or single predictions (given either
                  as a TheoryPrediction or a GlobalInfo object)
    :param strategy: combination strategy, can be conservative, moderate, aggressive
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
            ret = canCombine ( pA, predB, strategy )
            if ret == False:
                return False
        return True
    if type(predB) == list:
        for pB in predB:
            ret = canCombine ( predA, pB, strategy )
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
    if strategy == "conservative":
        return canCombineConservative ( predA, predB, elA, elB )
    if strategy == "moderate":
        return canCombineModerate ( predA, predB, elA, elB )
    if strategy != "aggressive":
        print ( "Error: strategy ``%s'' unknown" % strategy )
        return None
    return canCombineAggressive ( predA, predB, elA, elB )

def canCombineModerate ( globA, globB, elA, elB ):
    """ method that defines what we allow to combine, moderate version.
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
    allowCombination = { "ATLAS-SUSY-2013-02": [ "ATLAS-SUSY-2013-11" ],
                         "CMS-SUS-13-012": [ "CMS-SUS-13-007" ] }
    if anaidA in allowCombination.keys():
        if anaidB in allowCombination[anaidA]:
            return True
    if anaidB in allowCombination.keys():
        if anaidA in allowCombination[anaidB]:
            return True
    return False

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

def canCombineAggressive ( globA, globB, elA, elB ):
    """ method that defines what we allow to combine, aggressive version.
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
    allowCombinationATLAS8TeV = { "ATLAS-SUSY-2013-02": [ "ATLAS-SUSY-2013-04", "ATLAS-SUSY-2013-11" ],
              "ATLAS-CONF-2013-024": [ "ATLAS-CONF-2013-037", "ATLAS-CONF-2013-048", "ATLAS-CONF-2013-062", "ATLAS-CONF-2013-093", "ATLAS-SUSY-2013-09", "ATLAS-SUSY-2013-11", "ATLAS-SUSY-2013-12", "ATLAS-SUSY-2013-15", "ATLAS-SUSY-2013-19", "ATLAS-SUSY-2013-23" ],
              "ATLAS-CONF-2013-037": [ "ATLAS-CONF-2013-024", "ATLAS-CONF-2013-047", "ATLAS-CONF-2013-048", "ATLAS-CONF-2013-053", "ATLAS-CONF-2013-054", "ATLAS-SUSY-2013-02", "ATLAS-SUSY-2013-04", "ATLAS-SUSY-2013-05", "ATLAS-SUSY-2013-09", "ATLAS-SUSY-2013-11", "ATLAS-SUSY-2013-12", "ATLAS-SUSY-2013-16", "ATLAS-SUSY-2013-19" ],
              "ATLAS-CONF-2013-047": [ "ATLAS-CONF-2013-037", "ATLAS-CONF-2013-048", "ATLAS-CONF-2013-062", "ATLAS-CONF-2013-093", "ATLAS-SUSY-2013-11", "ATLAS-SUSY-2013-15", "ATLAS-SUSY-2013-09", "ATLAS-SUSY-2013-11", "ATLAS-SUSY-2013-12", "ATLAS-SUSY-2013-15", "ATLAS-SUSY-2013-19", "ATLAS-SUSY-2013-23" ],
              "ATLAS-CONF-2013-048": [ "ATLAS-CONF-2013-024", "ATLAS-CONF-2013-037", "ATLAS-CONF-2013-047", "ATLAS-CONF-2013-053", "ATLAS-CONF-2013-054", "ATLAS-CONF-2013-062", "ATLAS-CONF-2013-093", "ATLAS-SUSY-2013-02", "ATLAS-SUSY-2013-04", "ATLAS-SUSY-2013-05", "ATLAS-SUSY-2013-12", "ATLAS-SUSY-2013-15", "ATLAS-SUSY-2013-16", "ATLAS-SUSY-2013-18", "ATLAS-SUSY-2013-23", "ATLAS-SUSY-2013-02", "ATLAS-SUSY-2013-04", "ATLAS-SUSY-2013-05", "ATLAS-SUSY-2013-04", "ATLAS-SUSY-2013-02", "ATLAS-SUSY-2013-05", "ATLAS-SUSY-2013-15", "ATLAS-SUSY-2013-16", "ATLAS-SUSY-2013-18", "ATLAS-SUSY-2013-21", "ATLAS-SUSY-2014-03" ],
              "ATLAS-CONF-2013-053": [ "ATLAS-CONF-2013-062", "ATLAS-CONF-2013-093", "ATLAS-CONF-2013-007", "ATLAS-CONF-2013-089", "ATLAS-SUSY-2013-09", "ATLAS-SUSY-2013-11", "ATLAS-SUSY-2013-12", "ATLAS-SUSY-2013-15", "ATLAS-SUSY-2013-19", "ATLAS-SUSY-2013-23", "ATLAS-CONF-2013-037", "ATLAS-CONF-2013-048", "ATLAS-CONF-2013-062", "ATLAS-CONF-2013-093", "ATLAS-SUSY-2013-11", "ATLAS-SUSY-2013-15" ],
              "ATLAS-CONF-2013-054": [ "ATLAS-CONF-2013-062", "ATLAS-CONF-2013-093", "ATLAS-CONF-2013-007", "ATLAS-CONF-2013-089", "ATLAS-SUSY-2013-02", "ATLAS-SUSY-2013-05", "ATLAS-SUSY-2013-09", "ATLAS-SUSY-2013-11", "ATLAS-SUSY-2013-12", "ATLAS-SUSY-2013-15", "ATLAS-SUSY-2013-16", "ATLAS-SUSY-2013-19", "ATLAS-SUSY-2013-23", "ATLAS-CONF-2013-024", "ATLAS-CONF-2013-037", "ATLAS-CONF-2013-047", "ATLAS-CONF-2013-048", "ATLAS-CONF-2013-062", "ATLAS-CONF-2013-093", "ATLAS-SUSY-2013-11", "ATLAS-SUSY-2013-15", "ATLAS-SUSY-2013-21" ],
              "ATLAS-SUSY-2013-02": [ "ATLAS-CONF-2013-007", "ATLAS-CONF-2013-089", "ATLAS-CONF-2013-037", "ATLAS-CONF-2013-048", "ATLAS-CONF-2013-062", "ATLAS-CONF-2013-093", "ATLAS-SUSY-2013-09", "ATLAS-SUSY-2013-11", "ATLAS-SUSY-2013-12", "ATLAS-SUSY-2013-15", "ATLAS-SUSY-2013-19", "ATLAS-SUSY-2013-23" ],
              "ATLAS-SUSY-2013-04": [ "ATLAS-CONF-2013-007", "ATLAS-CONF-2013-089", "ATLAS-CONF-2013-037", "ATLAS-CONF-2013-048", "ATLAS-CONF-2013-062", "ATLAS-CONF-2013-093", "ATLAS-SUSY-2013-09", "ATLAS-SUSY-2013-11", "ATLAS-SUSY-2013-12", "ATLAS-SUSY-2013-15", "ATLAS-SUSY-2013-19", "ATLAS-SUSY-2013-23" ],
              "ATLAS-SUSY-2013-09": [ "ATLAS-SUSY-2013-02", "ATLAS-SUSY-2013-04", "ATLAS-SUSY-2013-05", "ATLAS-SUSY-2013-12", "ATLAS-SUSY-2013-15", "ATLAS-SUSY-2013-16", "ATLAS-SUSY-2013-18", "ATLAS-SUSY-2013-19", "ATLAS-SUSY-2013-23" ],
              "ATLAS-SUSY-2013-11": [ "ATLAS-SUSY-2013-02", "ATLAS-SUSY-2013-04", "ATLAS-SUSY-2013-05", "ATLAS-SUSY-2013-12", "ATLAS-SUSY-2013-15", "ATLAS-SUSY-2013-16", "ATLAS-SUSY-2013-18", "ATLAS-SUSY-2013-23" ],
              "ATLAS-SUSY-2013-12": [ "ATLAS-SUSY-2013-02", "ATLAS-SUSY-2013-04", "ATLAS-SUSY-2013-05", "ATLAS-SUSY-2013-08", "ATLAS-SUSY-2013-09", "ATLAS-SUSY-2013-11", "ATLAS-SUSY-2013-15", "ATLAS-SUSY-2013-16", "ATLAS-SUSY-2013-18", "ATLAS-SUSY-2013-19", "ATLAS-SUSY-2013-23", "ATLAS-CONF-2013-007", "ATLAS-CONF-2013-089" ],
              "ATLAS-SUSY-2013-15": [ "ATLAS-SUSY-2013-02", "ATLAS-SUSY-2013-04", "ATLAS-SUSY-2013-05", "ATLAS-SUSY-2013-09", "ATLAS-SUSY-2013-11", "ATLAS-SUSY-2013-12", "ATLAS-SUSY-2013-16", "ATLAS-SUSY-2013-19", "ATLAS-CONF-2013-007", "ATLAS-CONF-2013-089" ],
              "ATLAS-SUSY-2013-16": [ "ATLAS-SUSY-2013-02", "ATLAS-SUSY-2013-09", "ATLAS-SUSY-2013-11", "ATLAS-SUSY-2013-12", "ATLAS-SUSY-2013-15", "ATLAS-SUSY-2013-19", "ATLAS-SUSY-2013-23" ],
              "ATLAS-SUSY-2013-18": [ "ATLAS-SUSY-2013-09", "ATLAS-SUSY-2013-11", "ATLAS-SUSY-2013-12", "ATLAS-SUSY-2013-19", "ATLAS-SUSY-2013-23" ],
              "ATLAS-SUSY-2014-03": [ "ATLAS-CONF-2013-037", "ATLAS-CONF-2013-048", "ATLAS-CONF-2013-062", "ATLAS-CONF-2013-093", "ATLAS-SUSY-2013-11", "ATLAS-SUSY-2013-15", "ATLAS-SUSY-2013-21" ],
              "ATLAS-SUSY-2013-21": [ "ATLAS-CONF-2013-037", "ATLAS-CONF-2013-048", "ATLAS-CONF-2013-062", "ATLAS-CONF-2013-093", "ATLAS-SUSY-2013-11", "ATLAS-SUSY-2013-15" ],
              "ATLAS-SUSY-2013-18": [ "ATLAS-CONF-2013-048", "ATLAS-SUSY-2013-11" ],
              "ATLAS-SUSY-2013-16": [ "ATLAS-CONF-2013-037", "ATLAS-CONF-2013-048", "ATLAS-CONF-2013-062", "ATLAS-CONF-2013-093", "ATLAS-SUSY-2013-11", "ATLAS-SUSY-2013-15" ],
              "ATLAS-SUSY-2013-15": [ "ATLAS-CONF-2013-024", "ATLAS-CONF-2013-047", "ATLAS-CONF-2013-048", "ATLAS-CONF-2013-053", "ATLAS-CONF-2013-054", "ATLAS-SUSY-2013-02", "ATLAS-SUSY-2013-04", "ATLAS-SUSY-2013-05", "ATLAS-SUSY-2013-11", "ATLAS-SUSY-2013-16", "ATLAS-SUSY-2013-21", "ATLAS-SUSY-2014-03" ],
              "ATLAS-SUSY-2013-11": [ "ATLAS-CONF-2013-024", "ATLAS-CONF-2013-037", "ATLAS-CONF-2013-047", "ATLAS-CONF-2013-053", "ATLAS-CONF-2013-054", "ATLAS-CONF-2013-062", "ATLAS-CONF-2013-093", "ATLAS-SUSY-2013-02", "ATLAS-SUSY-2013-04", "ATLAS-SUSY-2013-05", "ATLAS-SUSY-2013-15", "ATLAS-SUSY-2013-16", "ATLAS-SUSY-2013-18", "ATLAS-SUSY-2013-21", "ATLAS-SUSY-2014-03" ],
              "ATLAS-SUSY-2013-05": [ "ATLAS-CONF-2013-007", "ATLAS-CONF-2013-089", "ATLAS-SUSY-2013-09", "ATLAS-SUSY-2013-11", "ATLAS-SUSY-2013-12", "ATLAS-SUSY-2013-15", "ATLAS-SUSY-2013-19", "ATLAS-SUSY-2013-23", "ATLAS-CONF-2013-037", "ATLAS-CONF-2013-048", "ATLAS-CONF-2013-062", "ATLAS-CONF-2013-093", "ATLAS-SUSY-2013-11", "ATLAS-SUSY-2013-15", "ATLAS-SUSY-2013-21" ],
              "ATLAS-CONF-2013-061": [ "ATLAS-SUSY-2013-21", ] ,

    }
    allowCombinationCMS8TeV = {
                         "CMS-SUS-13-012": [ "CMS-SUS-13-007", "CMS-SUS-13-013" ],
                         "CMS-SUS-12-024": [ "CMS-SUS-13-007", "CMS-SUS-13-013" ],
                         "CMS-SUS-13-007": [ "CMS-SUS-12-024", "CMS-SUS-13-012", "CMS-SUS-13-013", "CMS-EXO-12-026", "CMS-PAS-SUS-13-016", "CMS-PAS-SUS-13-018", "CMS-PAS-SUS-13-023" ],
                         "CMS-SUS-13-002": [ "CMS-EXO-12-026", "CMS-PAS-SUS-13-016", "CMS-PAS-SUS-13-018", "CMS-PAS-SUS-13-023", "CMS-SUS-12-024", "CMS-SUS-12-028", "CMS-SUS-13-007", "CMS-SUS-13-011", "CMS-SUS-13-012", "CMS-SUS-13-013", "CMS-SUS-13-019", "CMS-SUS-14-021" ],
                         "CMS-SUS-14-021": [ "CMS-EXO-12-026", "CMS-PAS-SUS-13-016", "CMS-SUS-13-002", "CMS-SUS-13-007", "CMS-SUS-13-011", "CMS-SUS-13-013" ],
                         "CMS-SUS-14-010": [ "CMS-EXO-12-026" ],
                         "CMS-SUS-13-013": [ "CMS-PAS-SUS-13-016", "CMS-EXO-12-026", "CMS-PAS-SUS-13-018", "CMS-PAS-SUS-13-023", "CMS-SUS-12-024", "CMS-SUS-12-028", "CMS-SUS-13-002", "CMS-SUS-13-007", "CMS-SUS-13-011", "CMS-SUS-13-012", "CMS-SUS-13-019", "CMS-SUS-14-021" ],
                         "CMS-SUS-13-012": [ "CMS-EXO-12-026", "CMS-PAS-SUS-13-016", "CMS-SUS-13-002", "CMS-SUS-13-007", "CMS-SUS-13-011", "CMS-SUS-13-013" ],
                         "CMS-SUS-12-028": [ "CMS-EXO-12-026", "CMS-PAS-SUS-13-016", "CMS-SUS-13-002", "CMS-SUS-13-007", "CMS-SUS-13-011", "CMS-SUS-13-013" ],
                         "CMS-SUS-12-024": [ "CMS-EXO-12-026", "CMS-PAS-SUS-13-016", "CMS-SUS-13-002", "CMS-SUS-13-007", "CMS-SUS-13-011", "CMS-SUS-13-013" ],
                         "CMS-PAS-SUS-13-016": [ "CMS-SUS-13-013", "CMS-EXO-12-026", "CMS-PAS-SUS-13-018", "CMS-PAS-SUS-13-023", "CMS-SUS-12-024", "CMS-SUS-12-028", "CMS-SUS-13-002", "CMS-SUS-13-007", "CMS-SUS-13-011", "CMS-SUS-13-012", "CMS-SUS-13-019", "CMS-SUS-14-021" ],
                         "CMS-SUS-13-015": [ "CMS-EXO-12-026", "CMS-PAS-SUS-13-016", "CMS-SUS-13-002", "CMS-SUS-13-007", "CMS-SUS-13-011", "CMS-SUS-13-013" ],
                         "CMS-EXO-13-006": [ "CMS-SUS-*", "CMS-PAS-SUS-*" ],
                         "CMS-EXO-12-026": [ "CMS-SUS-*", "CMS-PAS-SUS-*" ],
    }

    allowCombinationCMS13TeV = {
        "CMS-SUS-16-009": [ "CMS-SUS-17-009", "CMS-PAS-EXO-16-036", "CMS-SUS-16-051", "CMS-SUS-17-001", "CMS-SUS-16-042", "CMS-SUS-16-041", "CMS-SUS-16-039", "CMS-SUS-16-037", "CMS-SUS-16-035", "CMS-SUS-16-034" ], ## all hadronic
        "CMS-SUS-17-004": [ "CMS-PAS-EXO-16-036" ], # Chargino-neutralino production with WZ topology, combination paper
        "CMS-SUS-17-005": [ "CMS-SUS-19-006", "CMS-SUS-19-006-2", "CMS-SUS-16-009", "CMS-SUS-17-006", "CMS-SUS-16-033", "CMS-SUS-17-009", "CMS-PAS-EXO-16-036", "CMS-SUS-16-051", "CMS-SUS-17-001", "CMS-SUS-16-042", "CMS-SUS-16-041", "CMS-SUS-16-039", "CMS-SUS-16-037", "CMS-SUS-16-035", "CMS-SUS-16-034"  ], # multijets + Etmiss, top tagging
        "CMS-SUS-17-006": [ "CMS-SUS-17-009", "CMS-PAS-EXO-16-036", "CMS-SUS-16-051", "CMS-SUS-17-001", "CMS-SUS-16-042", "CMS-SUS-16-041", "CMS-SUS-16-039", "CMS-SUS-16-037", "CMS-SUS-16-035", "CMS-SUS-16-034" ], # High momentum Higgs Boson+ Etmiss, no leptons in SR
        "CMS-SUS-17-009": [ "CMS-SUS-16-009", "CMS-SUS-17-005", "CMS-SUS-17-006", "CMS-SUS-18-002", "CMS-SUS-16-033", "CMS-SUS-16-036", "CMS-SUS-16-032", "CMS-SUS-16-045", "CMS-SUS-16-046", "CMS-SUS-16-047", "CMS-SUS-16-049", "CMS-SUS-16-050" ], # leptons + Etmiss
        "CMS-SUS-17-010": [ "CMS-SUS-16-033", "CMS-SUS-19-006", "CMS-SUS-19-006-2", "CMS-SUS-16-049", "CMS-PAS-EXO-16-036", "CMS-SUS-16-051", "CMS-SUS-16-042", "CMS-SUS-16-041", "CMS-SUS-16-039", "CMS-SUS-16-037", "CMS-SUS-16-035", "CMS-SUS-16-036", "CMS-SUS-17-005", "CMS-SUS-16-032", "CMS-SUS-18-002", "CMS-SUS-16-046", "CMS-SUS-19-006", "CMS-PAS-SUS-17-004", "CMS-SUS-16-009", "CMS-SUS-16-043", "CMS-SUS-16-047", "CMS-PAS-SUS-16-022", "CMS-SUS-17-006" ], # hadronic stop, 2 OS leptons
        "CMS-SUS-18-002": [ "CMS-SUS-16-033", "CMS-SUS-19-006", "CMS-SUS-19-006-2", "CMS-SUS-17-005", "CMS-SUS-17-009", "CMS-PAS-EXO-16-036", "CMS-SUS-16-051", "CMS-SUS-17-001", "CMS-SUS-16-042", "CMS-SUS-16-041", "CMS-SUS-16-039", "CMS-SUS-16-037", "CMS-SUS-16-035", "CMS-SUS-16-034" ], # photon, jets, b-jets+ Etmiss, top tagging, lepton veto
        "CMS-SUS-19-006": [ "CMS-SUS-17-009", "CMS-PAS-EXO-16-036", "CMS-SUS-16-051", "CMS-SUS-17-001", "CMS-SUS-16-042", "CMS-SUS-16-041", "CMS-SUS-16-039", "CMS-SUS-16-037", "CMS-SUS-16-035", "CMS-SUS-16-034", "CMS-PAS-SUS-16-022", "CMS-PAS-SUS-17-004", "CMS-SUS-16-043", "CMS-SUS-18-002", "CMS-SUS-17-005" ], # 0L + multijets with MHT
        "CMS-PAS-EXO-16-036": [ "CMS-PAS-SUS-*", "CMS-SUS-*" ], # hscp search
        "CMS-PAS-SUS-16-022": [ "CMS-PAS-EXO-16-036", "CMS-PAS-SUS-16-052", "CMS-SUS-16-032", "CMS-SUS-16-033", "CMS-SUS-16-034", "CMS-SUS-16-035", "CMS-SUS-16-036", "CMS-SUS-16-037", "CMS-SUS-16-042", "CMS-SUS-16-045", "CMS-SUS-16-046", "CMS-SUS-16-047", "CMS-SUS-16-049", "CMS-SUS-16-050", "CMS-SUS-16-051", "CMS-SUS-17-001" ], # >= 3 leptons + Etmiss
        "CMS-SUS-16-051": [ "CMS-PAS-EXO-16-036", "CMS-PAS-SUS-16-022", "CMS-PAS-SUS-16-052", "CMS-PAS-SUS-17-004", "CMS-SUS-16-032", "CMS-SUS-16-033", "CMS-SUS-16-034", "CMS-SUS-16-035", "CMS-SUS-16-036", "CMS-SUS-16-039", "CMS-SUS-16-041", "CMS-SUS-16-045", "CMS-SUS-16-046", "CMS-SUS-16-047", "CMS-SUS-16-049", "CMS-SUS-16-050", "CMS-SUS-17-001", "CMS-PAS-SUS-16-052-agg" ], # 1L stop
        "CMS-SUS-16-050": [ "CMS-PAS-EXO-16-036", "CMS-PAS-SUS-16-022", "CMS-PAS-SUS-17-004", "CMS-SUS-16-034", "CMS-SUS-16-035", "CMS-SUS-16-037", "CMS-SUS-16-039", "CMS-SUS-16-041", "CMS-SUS-16-042", "CMS-SUS-16-051" ], # 0L + top tag
        "CMS-SUS-16-049": [ "CMS-PAS-EXO-16-036", "CMS-PAS-SUS-16-022", "CMS-PAS-SUS-17-004", "CMS-SUS-16-034", "CMS-SUS-16-035", "CMS-SUS-16-037", "CMS-SUS-16-039", "CMS-SUS-16-041", "CMS-SUS-16-042", "CMS-SUS-16-051" ], # All hadronic stop
        "CMS-SUS-16-036": [ "CMS-PAS-EXO-16-036", "CMS-PAS-SUS-16-022", "CMS-PAS-SUS-17-004", "CMS-SUS-16-034", "CMS-SUS-16-035", "CMS-SUS-16-037", "CMS-SUS-16-039", "CMS-SUS-16-041", "CMS-SUS-16-042", "CMS-SUS-16-051" ], # 0L + jets + Etmiss (using MT2)
    }
    allowCombinationATLAS13TeV = {
        "ATLAS-SUSY-2016-15": [ "ATLAS-SUSY-2016-16", "ATLAS-SUSY-2016-24", "ATLAS-SUSY-2018-06", "ATLAS-SUSY-2018-32", "ATLAS-SUSY-2015-02", "ATLAS-SUSY-2015-09" ], # 0L stop
        "ATLAS-SUSY-2016-16": [ "ATLAS-SUSY-2016-15", "ATLAS-SUSY-2016-24", "ATLAS-SUSY-2017-02", "ATLAS-SUSY-2016-17", "ATLAS-SUSY-2016-14", "ATLAS-SUSY-2018-06", "ATLAS-SUSY-2018-32", "ATLAS-SUSY-2015-01", "ATLAS-SUSY-2015-09", "ATLAS-SUSY-2016-07", "ATLAS-SUSY-2015-06" ], # 1L stop
        "ATLAS-SUSY-2016-24": [  "ATLAS-SUSY-2016-15", "ATLAS-SUSY-2016-16", "ATLAS-SUSY-2016-28", "ATLAS-SUSY-2017-02", "ATLAS-SUSY-2018-04", "ATLAS-SUSY-2015-01", "ATLAS-SUSY-2015-02", "ATLAS-SUSY-2016-26", "ATLAS-SUSY-2015-02", "ATLAS-SUSY-2015-06", "ATLAS-SUSY-2016-07" ], # 2+ leptons (e,mu) + jets + Etmiss
        "ATLAS-SUSY-2016-27": [], #  "ATLAS-SUSY-2016-16", "ATLAS-SUSY-2016-24", "ATLAS-SUSY-2018-06", "ATLAS-SUSY-2018-32", "ATLAS-SUSY-2015-02", "ATLAS-SUSY-2015-09"  ], # jets + photon + Etmiss
        "ATLAS-SUSY-2017-03": [ "ATLAS-SUSY-2016-16", "ATLAS-SUSY-2016-28", "ATLAS-SUSY-2019-08", "ATLAS-SUSY-2017-02", "ATLAS-SUSY-2018-31", "ATLAS-SUSY-2016-15" ], # Multilepton EWK searches (2 or 3 leptons)
        "ATLAS-SUSY-2016-26": [ "ATLAS-SUSY-2019-08", "ATLAS-SUSY-2016-17", "ATLAS-SUSY-2016-33", "ATLAS-SUSY-2017-03", "ATLAS-SUSY-2016-16", "ATLAS-SUSY-2016-14" ], # >=2 c jets + MET, lepton veto
        "ATLAS-SUSY-2016-28": [ "ATLAS-SUSY-2016-24", "ATLAS-SUSY-2018-06", "ATLAS-SUSY-2018-32", "ATLAS-SUSY-2015-02", "ATLAS-SUSY-2015-09" ], # 2b
        "ATLAS-SUSY-2017-01": [], # EWK WH(bb) + Etmiss, all lepton multiplicities, with and without photon
        "ATLAS-SUSY-2017-02": [ "ATLAS-SUSY-2016-16", "ATLAS-SUSY-2016-24", "ATLAS-SUSY-2018-06", "ATLAS-SUSY-2018-32", "ATLAS-SUSY-2015-02", "ATLAS-SUSY-2015-09" ], # 0L + jets + Etmiss
        "ATLAS-SUSY-2018-04": [ "ATLAS-SUSY-2016-28", "ATLAS-SUSY-2018-31", "ATLAS-SUSY-2016-16", "ATLAS-SUSY-2015-02", "ATLAS-SUSY-2015-09" ], # 2 hadronic taus
        "ATLAS-SUSY-2018-06": [ "ATLAS-SUSY-2016-15", "ATLAS-SUSY-2016-16", "ATLAS-SUSY-2016-28", "ATLAS-SUSY-2017-02", "ATLAS-SUSY-2018-04", "ATLAS-SUSY-2015-01", "ATLAS-SUSY-2015-02", "ATLAS-SUSY-2016-26", "ATLAS-SUSY-2015-02", "ATLAS-SUSY-2015-06", "ATLAS-SUSY-2016-07" ], # 3 leptons EW-ino
        "ATLAS-SUSY-2018-16": [], # EW production in models with compressed mass spectra, low momentum leptons
        "ATLAS-SUSY-2016-33": [ "ATLAS-SUSY-2016-16", "ATLAS-SUSY-2016-15", "ATLAS-SUSY-2016-28", "ATLAS-SUSY-2018-06", "ATLAS-SUSY-2018-31", "ATLAS-SUSY-2017-02" ], # 2 OSSF leptons + Etmiss
        "ATLAS-SUSY-2018-31": [ "ATLAS-SUSY-2016-14", "ATLAS-SUSY-2016-16", "ATLAS-SUSY-2016-24", "ATLAS-SUSY-2018-06", "ATLAS-SUSY-2018-32", "ATLAS-SUSY-2015-02", "ATLAS-SUSY-2015-09" ], # higgs + b-jets + Etmiss, lepton veto
        "ATLAS-SUSY-2018-32": [ "ATLAS-SUSY-2016-15", "ATLAS-SUSY-2016-16", "ATLAS-SUSY-2016-28", "ATLAS-SUSY-2017-02", "ATLAS-SUSY-2018-06", "ATLAS-SUSY-2015-01", "ATLAS-SUSY-2015-02", "ATLAS-SUSY-2016-26", "ATLAS-SUSY-2015-02", "ATLAS-SUSY-2015-06", "ATLAS-SUSY-2016-07" ], # 2 leptons + Etmiss
        "ATLAS-SUSY-2019-08": [ "ATLAS-SUSY-2016-14", "ATLAS-SUSY-2016-17", "ATLAS-SUSY-2016-33", "ATLAS-SUSY-2016-15", "ATLAS-SUSY-2016-24", "ATLAS-SUSY-2018-31", "ATLAS-SUSY-2017-02", "ATLAS-SUSY-2018-06", "ATLAS-SUSY-2018-32", "ATLAS-SUSY-2015-01", "ATLAS-SUSY-2015-09", "ATLAS-SUSY-2016-07", "ATLAS-SUSY-2015-06" ], # l + higgs + Etmiss
        "ATLAS-SUSY-2015-01": [ "ATLAS-SUSY-2018-04", "ATLAS-SUSY-2015-02", "ATLAS-SUSY-2015-09", "ATLAS-SUSY-2016-14", "ATLAS-SUSY-2016-17", "ATLAS-SUSY-2016-33", "ATLAS-SUSY-2017-03", "ATLAS-SUSY-2015-02" ], # 2 b-jets + Etmiss
        "ATLAS-SUSY-2015-02": [ "ATLAS-SUSY-2015-01", "ATLAS-SUSY-2015-09", "ATLAS-SUSY-2016-14", "ATLAS-SUSY-2016-17", "ATLAS-SUSY-2016-26", "ATLAS-SUSY-2016-33", "ATLAS-SUSY-2017-03", "ATLAS-SUSY-2015-06", "ATLAS-SUSY-2016-07" ], # single lepton stop
        "ATLAS-SUSY-2015-09": [ "ATLAS-SUSY-2015-01", "ATLAS-SUSY-2015-02", "ATLAS-SUSY-2016-26", "ATLAS-SUSY-2015-02", "ATLAS-SUSY-2015-06", "ATLAS-SUSY-2016-07" ], # jets + 2 SS leptons or >= 3 leptons
        "ATLAS-SUSY-2016-08": [ ], # ?
        "ATLAS-SUSY-2016-07": [ "ATLAS-SUSY-2015-02", "ATLAS-SUSY-2015-09", "ATLAS-SUSY-2016-14", "ATLAS-SUSY-2016-17", "ATLAS-SUSY-2016-33", "ATLAS-SUSY-2017-03", "ATLAS-SUSY-2015-02" ], # 0L + jets + Etmiss
        "ATLAS-SUSY-2015-06": [ "ATLAS-SUSY-2015-02", "ATLAS-SUSY-2015-09", "ATLAS-SUSY-2016-14", "ATLAS-SUSY-2016-17", "ATLAS-SUSY-2016-33", "ATLAS-SUSY-2017-03", "ATLAS-SUSY-2015-02" ], # 0 leptons + 2-6 jets + Etmiss
    }
    allowCombination = {}
    allowCombination.update ( allowCombinationATLAS8TeV )
    allowCombination.update ( allowCombinationCMS8TeV )
    allowCombination.update ( allowCombinationATLAS13TeV )
    allowCombination.update ( allowCombinationCMS13TeV )
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

def canCombineConservative ( globA, globB, elA, elB ):
    """ method that defines what we allow to combine, conservative version.
         """
    if globA.sqrts != globB.sqrts:
        return True
    if getExperimentName(globA) != getExperimentName(globB):
        return True
    return False

def getSummary():
    from smodels.experiment.databaseObj import Database
    # dbpath = "official"
    dbpath = "<rundir>/database.pcl"
    print ( "[analysisCombiner] checking %s" % dbpath )
    db = Database ( dbpath )
    results = db.getExpResults()
    strategy="aggressive"
    ana1 = "CMS-SUS-16-042"
    ana2 = "CMS-SUS-16-033"
    canC = canCombine ( ana1, ana2, strategy, results )
    print ( "[analysisCombiner] can combine %s with %s: %s" % ( ana1, ana2, str(canC) ) )
    ctr,combinable=0,0
    for x,e in enumerate(results):
        for y,f in enumerate(results):
            if y <= x:
                continue
            ctr += 1
            isUn = canCombine ( e.globalInfo, f.globalInfo, strategy )
            combinable+=isUn
    print ( "[analysisCombiner] can combine %d/%d pairs of results" % ( combinable, ctr ) )

def checkOneAnalysis():
    import argparse
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
        cc = canCombine ( info, er.globalInfo, "aggressive" )
        # cc = canCombine ( pred, er.globalInfo )
        if cc:
            combs.add ( Id )
        else:
            nocombs.add ( Id )
    print ( "can combine with: " )
    for Id in combs:
        pname = pnames[Id]
        print ( " `- %s: %s" % ( Id, pname ) )
    print ( "cannot combine with: " )
    for Id in nocombs:
        pname = pnames[Id]
        print ( " `- %s: %s" % ( Id, pname ) )

if __name__ == "__main__":
    checkOneAnalysis()
    # getSummary()
