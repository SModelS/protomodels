#!/usr/bin/env python3

"""
.. module:: convert
   :synopsis: used to create info.txt and the <txname>.txt files.

"""
import sys
import os
import argparse

argparser = argparse.ArgumentParser(description =  
'create info.txt, txname.txt, twiki.txt and sms.py')
argparser.add_argument ('-utilsPath', '--utilsPath', 
help = 'path to the package smodels_utils',\
type = str )
argparser.add_argument ('-smodelsPath', '--smodelsPath', 
help = 'path to the package smodels_utils',\
type = str )
args = argparser.parse_args()

if args.utilsPath:
    utilsPath = args.utilsPath
else:
    databaseRoot = '../../../'
    sys.path.append(os.path.abspath(databaseRoot))
    from utilsPath import utilsPath
    utilsPath = databaseRoot + utilsPath
if args.smodelsPath:
    sys.path.append(os.path.abspath(args.smodelsPath))

sys.path.append(os.path.abspath(utilsPath))
from smodels_utils.dataPreparation.inputObjects import MetaInfoInput,DataSetInput
from smodels_utils.dataPreparation.databaseCreation import databaseCreator
from smodels_utils.dataPreparation.massPlaneObjects import x, y, z



#+++++++ global info block ++++++++++++++
info = MetaInfoInput('ATLAS-SUSY-2016-33')
info.url = "https://atlas.web.cern.ch/Atlas/GROUPS/PHYSICS/PAPERS/SUSY-2016-33/"
info.sqrts = 13
info.lumi = 36.1
info.prettyName = "2 OSSF leptons + Etmiss"
info.private = False
info.arxiv =  'https://arxiv.org/abs/1805.11381'
info.contact = 'ATLAS collaboration'
info.publication ='in process of being published'
info.comment = "For T6ZZ we assume that the validation was performed with the squarks of only one 'chirality'"

#+++++++ dataset block ++++++++++++++
dataset = DataSetInput('data')
dataset.setInfo(dataType = 'upperLimit', dataId = None)

#+++++++ next txName block ++++++++++++++
T5ZZ = dataset.addTxName('T5ZZ')
T5ZZ.checked = ''
T5ZZ.constraint ="[[[jet,jet],[Z]],[[jet,jet],[Z]]]"
T5ZZ.conditionDescription = None
T5ZZ.condition = None
T5ZZ.source = "ATLAS"
#+++++++ next mass plane block ++++++++++++++
T5ZZ_1 = T5ZZ.addMassPlane(2*[[x, 0.5*(x+y), y]])
T5ZZ_1.figure = 'Fig.13a'
T5ZZ_1.figureUrl ='https://atlas.web.cern.ch/Atlas/GROUPS/PHYSICS/PAPERS/SUSY-2016-33/fig_13a.png'
T5ZZ_1.dataUrl ='https://www.hepdata.net/record/ins1675352?version=1&table=Cross%20section%20UL%206'
T5ZZ_1.setSources(dataLabels= ['expExclusion', 'obsExclusion', 'upperLimits'],
                 dataFiles= ['orig/Exp_Excl_3.csv', 'orig/Obs_Excl_3.csv', 'orig/Obs_UL_6.csv'],
                 units = [ None, None, 'fb' ],
                 dataFormats= ['csv', 'csv', 'csv'])
#+++++++ next mass plane block ++++++++++++++
T5ZZ_2 = T5ZZ.addMassPlane(2*[[x, y, 1.]])
T5ZZ_2.figure = 'Fig.14a'
T5ZZ_2.figureUrl ='https://atlas.web.cern.ch/Atlas/GROUPS/PHYSICS/PAPERS/SUSY-2016-33/fig_14a.png'
T5ZZ_2.dataUrl = 'https://www.hepdata.net/record/ins1675352?version=1&table=Cross%20section%20UL%201'
T5ZZ_2.setSources(dataLabels= ['expExclusion', 'obsExclusion', 'upperLimits'],
                 dataFiles= ['orig/Exp_Excl_5.csv', 'orig/Obs_Excl_5.csv', 'orig/Obs_UL_1.csv'],
                 units = [ None, None, 'fb' ],
                 dataFormats= ['csv', 'csv', 'csv'])
#+++++++ next mass plane block ++++++++++++++
T5ZZ_3 = T5ZZ.addMassPlane(2*[[x, y+100., y]])
T5ZZ_3.figure = 'Fig.15'
T5ZZ_3.figureUrl ='https://atlas.web.cern.ch/Atlas/GROUPS/PHYSICS/PAPERS/SUSY-2016-33/fig_15.png'
T5ZZ_3.dataUrl = 'https://www.hepdata.net/record/ins1675352?version=1&table=Cross%20section%20UL%203'
T5ZZ_3.setSources(dataLabels= ['expExclusion', 'obsExclusion', 'upperLimits'],
                 dataFiles= ['orig/Exp_Excl_7.csv', 'orig/Obs_Excl_7.csv', 'orig/Obs_UL_3.csv'],
                 units = [ None, None, 'fb' ],
                 dataFormats= ['csv', 'csv', 'csv'])

# +++++++ next txName block ++++++++++++++
T6ZZ = dataset.addTxName('T6ZZ')
T6ZZ.checked = ''
T6ZZ.constraint ="[[[jet],[Z]],[[jet],[Z]]]"
T6ZZ.conditionDescription = None
T6ZZ.condition = None
T6ZZ.source = "ATLAS"
#+++++++ next mass plane block ++++++++++++++
T6ZZ_1 = T6ZZ.addMassPlane(2*[[x, y, 1.]])
T6ZZ_1.figure = 'Fig.14b'
T6ZZ_1.figureUrl = "https://atlas.web.cern.ch/Atlas/GROUPS/PHYSICS/PAPERS/SUSY-2016-33/fig_14b.png"
T6ZZ_1.dataUrl = 'https://www.hepdata.net/record/ins1675352?version=1&table=Cross%20section%20UL%202'
T6ZZ_1.setSources(dataLabels= ['expExclusion', 'obsExclusion', 'upperLimits'],
                 dataFiles= ['orig/Exp_Excl_2.csv', 'orig/Obs_Excl_2.csv', 
                             'orig/Obs_UL_2.csv'],
                 units = [ None, None, 'fb' ],
                 dataFormats= ['csv', 'csv', 'csv'])

databaseCreator.create()
