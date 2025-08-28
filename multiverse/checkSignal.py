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
    dT = ["upperLimit"] 
    orig_er = orig_db.getExpResults ( analysisIDs = [analysis], txnames=[txname],
                                      dataTypes=dT)[0]
    new_er = new_db.getExpResults ( analysisIDs = [analysis], txnames=[txname],
                                    dataTypes=dT)[0]
    orig_ul = orig_er.getUpperLimitFor(mass=masses,txname=txname)
    new_ul = new_er.getUpperLimitFor(mass=masses,txname=txname)
    print ( analysis, masses, txname )
    print ( f"orig_ul {orig_ul}, new_ul {new_ul}" )
    # import sys, IPython; IPython.embed( colors = "neutral" )

def checkSignal():
    from smodels.experiment.databaseObj import Database
    # dbpath = "official"
    # dbpath = f"{os.environ['HOME']}/git/smodels-database"  )
    dbpath = "./original.pcl"
    orig_db = Database ( dbpath )
    print ( orig_db )
    # ./expResModifier.py -R ./ -d original.pcl -s stop1 -P model_stop.dict
    new_db = Database ( "./stop1.pcl" )
    analysis = "ATLAS-SUSY-2018-08"
    txname = "T2tt"
    masses = [[500*GeV, 100*GeV],[500*GeV, 100*GeV]]
    analyseUL ( orig_db, new_db, analysis, txname, masses )
    masses = [[550*GeV, 120*GeV],[550*GeV, 120*GeV]]
    analyseUL ( orig_db, new_db, analysis, txname, masses )
    masses = [[1000*GeV, 400*GeV],[1000*GeV, 400*GeV]]
    analyseUL ( orig_db, new_db, analysis, txname, masses )



if __name__ == "__main__":
    checkSignal()
