#!/usr/bin/env python3

"""
.. module:: checkSignal
   :synopsis: simple snippet to check if signal was injected correctly

.. moduleauthor:: Wolfgang Waltenberger <wolfgang.waltenberger@gmail.com>

"""

import sys, os
sys.path.insert(0,"../")
from smodels.base.physicsUnits import GeV

def analyseUL ( orig_db, new_db, analysis : str, txname : str, masses : list ):
    """ analyse this upper limit result """
    #dT = ["upperLimit"] 
    dT = ["all"] 
    orig_ers = orig_db.getExpResults ( analysisIDs = [analysis], txnames=[txname],
                                      dataTypes=dT)
    if len(orig_ers) == 0:
        print ( f"no results for {analysis}:{txname}:{dT}" )
        return
    orig_er = orig_ers[0]
    new_er = new_db.getExpResults ( analysisIDs = [analysis], txnames=[txname],
                                    dataTypes=dT)[0]
    orig_ul = orig_er.getUpperLimitFor(mass=masses,txname=txname)
    new_ul = new_er.getUpperLimitFor(mass=masses,txname=txname)
    print ( analysis, masses, txname )
    print ( f"orig_ul {orig_ul}, new_ul {new_ul}" )
    if False:
        import sys, IPython; IPython.embed( colors = "neutral" ); sys.exit()

def analyseCombined ( orig_db, new_db, analysis : str, slhafile : str ):
    """ analyse this upper limit result """
    from smodels.decomposition import decomposer
    from smodels.base.physicsUnits import fb, GeV, TeV
    from smodels.matching.theoryPrediction import theoryPredictionsFor
    from smodels.tools.particlesLoader import load
    from smodels.share.models.SMparticles import SMList
    from smodels.base.model import Model
    BSMList = load()
    model = Model(BSMparticles=BSMList, SMparticles=SMList)
    model.updateParticles(inputFile=slhafile,
                          ignorePromptQNumbers = ['eCharge','colordim','spin'])
    # Set main options for decomposition
    sigmacut = 0.001 * fb
    mingap = 5.*GeV

    # Decompose model
    topDict = decomposer.decompose(model, sigmacut,
                                   massCompress=True, invisibleCompress=True,
                                   minmassgap=mingap)
    #dT = ["upperLimit"] 
    dT = ["all"] 
    orig_ers = orig_db.getExpResults ( analysisIDs = [analysis], txnames=["all"],
                                      dataTypes=dT)
    orig_er = orig_ers[0]
    
    allPredictions = theoryPredictionsFor(orig_db, topDict, combinedResults=True)
    print ( f"original {orig_db.databaseVersion}" )
    for p in allPredictions:
        print ( p.getUpperLimit() )
    new_ers = new_db.getExpResults ( analysisIDs = [analysis], txnames=["all"],
                                    dataTypes=dT)
    new_er = new_ers[0]
    newPredictions = theoryPredictionsFor(new_db, topDict, combinedResults=True)
    print ( f"new {orig_db.databaseVersion}" )
    for p in newPredictions:
        print ( p.getUpperLimit() )
    import sys, IPython; IPython.embed( colors = "neutral" ); sys.exit()


def checkSignal():
    from smodels.experiment.databaseObj import Database
    # dbpath = "official"
    # dbpath = f"{os.environ['HOME']}/git/smodels-database"  )
    dbpath = "./original.pcl"
    orig_db = Database ( dbpath )
    print ( orig_db )
    # ./expResModifier.py -R ./ -d original.pcl -s stop1 -P model_stop.dict
    new_db = Database ( "./signal.pcl" )
    analysis = "ATLAS-SUSY-2018-05-ewk"
    #slhafile = "TChiWZ_300_50.slha"
    #analyseCombined ( orig_db, new_db, analysis, slhafile )
    analysis = "ATLAS-SUSY-2018-05"
    txname = "TChiWZ"
    masses = [[300*GeV, 100*GeV],[300*GeV, 100*GeV]]
    analyseUL ( orig_db, new_db, analysis, txname, masses )
    masses = [[550*GeV, 120*GeV],[550*GeV, 120*GeV]]
    analyseUL ( orig_db, new_db, analysis, txname, masses )
    masses = [[1000*GeV, 400*GeV],[1000*GeV, 400*GeV]]
    analyseUL ( orig_db, new_db, analysis, txname, masses )



if __name__ == "__main__":
    checkSignal()
