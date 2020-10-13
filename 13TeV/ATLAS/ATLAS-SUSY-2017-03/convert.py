#!/usr/bin/env python

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
info = MetaInfoInput('ATLAS-SUSY-2017-03')
info.url = "https://atlas.web.cern.ch/Atlas/GROUPS/PHYSICS/PAPERS/SUSY-2017-03/"
info.sqrts = 13
info.lumi = 36.1
info.prettyName = "Multilepton EWK searches"
info.private = False
info.arxiv =  'https://arxiv.org/abs/1806.02293'
info.contact = 'ATLAS collaboration'
info.publication ='Submitted to PRD'

#+++++++ dataset block ++++++++++++++
dataset = DataSetInput('data')
dataset.setInfo(dataType = 'upperLimit', dataId = None)

#+++++++ next txName block ++++++++++++++
TChiWZ = dataset.addTxName('TChiWZ')
TChiWZ.checked = ''
TChiWZ.constraint ="[[[W]],[[Z]]]"
TChiWZ.conditionDescription = None
TChiWZ.condition = None
TChiWZ.source = "ATLAS"

#+++++++ next mass plane block ++++++++++++++
TChiWZ_1 = TChiWZ.addMassPlane(2*[[x, y]])
TChiWZ_1.figure = 'Fig.13c'
TChiWZ_1.figureUrl ='https://atlas.web.cern.ch/Atlas/GROUPS/PHYSICS/PAPERS/SUSY-2017-03/fig_13c.png'
TChiWZ_1.dataUrl ='https://www.hepdata.net/record/ins1676551?version=1&table=Cross-section%20UL%20combined'
TChiWZ_1.setSources(dataLabels= ['expExclusion', 'obsExclusion', 'upperLimits'],
                 coordinates= [ {x: 0, y: 1, 'value': None}, {x: 0, y: 1, 'value': None},  {x : 1, y: 0, 'value' :2} ],
		 units = [None, None, 'fb'],
		 dataFiles= ['orig/Exp_Excl_Comb.csv', 'orig/Obs_Excl_Comb.csv', 'orig/Obs_UL_Comb.csv'],
                 dataFormats= ['csv', 'csv', 'csv'])


databaseCreator.create()

