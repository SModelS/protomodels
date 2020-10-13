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
info = MetaInfoInput('ATLAS-SUSY-2016-15')
info.url 		 = 'https://atlas.web.cern.ch/Atlas/GROUPS/PHYSICS/PAPERS/SUSY-2016-15/'
info.sqrts 		 = 13
info.lumi 		 = 36.1
info.prettyName  = '0L stop'
info.private 	 = False
info.arxiv 		 = 'https://arxiv.org/abs/1709.04183'
info.contact 	 = 'atlas-phys-susy-conveners@cern.ch'
info.publication = 'https://link.springer.com/article/10.1007/JHEP12(2017)085'

#+++++++ dataset block ++++++++++++++
dataset = DataSetInput('data')
dataset.setInfo(dataType = 'upperLimit', dataId = None)

T2tt 							= dataset.addTxName('T2tt')
T2tt.checked 					= 'False'
T2tt.constraint 				= '[[[t]],[[t]]]'
T2tt.conditionDescription 		= None
T2tt.condition 					= None
T2tt.massConstraint				= [['dm >= 169.0'], ['dm >= 169.0']]
T2tt.source 					= 'ATLAS'
#+++++++ next txName block ++++++++++++++
T2ttoff 						= dataset.addTxName('T2ttoff')
T2ttoff.checked 				= 'False'
T2ttoff.constraint 				= '[[[b, W]],[[b, W]]]'
T2ttoff.conditionDescription 	= None
T2ttoff.condition 				= None
T2ttoff.massConstraint 			= [['80 <= dm < 169.0'], ['80 <= dm < 169.0']]
T2ttoff.source 					= 'ATLAS'
#+++++++ next mass plane block ++++++++++++++
T2bbffff						= dataset.addTxName('T2bbffff')
T2bbffff.checked				= 'False'
T2bbffff.constraint				= "9./4.*[[['b', 'jet','jet']],[['b', 'jet','jet']]]"
T2bbffff.conditionDescription 	= None
T2bbffff.condition				= None
T2bbffff.source					= 'ATLAS'
T2bbffff.massConstraint			= [['dm < 80'], ['dm < 80']]
#+++++++ next mass plane block ++++++++++++++
T2tt_1 							= T2tt.addMassPlane(2*[[x, y]])
T2tt_1.figure 					= 'figaux_12'
T2tt_1.figureUrl 				= 'https://atlas.web.cern.ch/Atlas/GROUPS/PHYSICS/PAPERS/SUSY-2016-15/figaux_12.png'
T2tt_1.dataUrl					= 'https://www.hepdata.net/record/ins1623207?version=7&table=X-section U.L. direcTT'
T2tt_1.setSources(dataLabels	= ['expExclusion', 'obsExclusion', 'upperLimits'],
                 	dataFiles	= ['orig/ExpectedexclusioncontourdirectTT.csv', 'orig/ObservedexclusioncontourdirectTT.csv', 'orig/X-sectionU.L.direcTT.csv'],
					coordinates = [ {x: 0, y: 1, 'value': None}, {x: 0, y: 1, 'value': None},  {x : 1, y: 0, 'value' :2} ],
                 		units 	= [ None, None, 'fb' ],
                 	dataFormats	= ['csv', 'csv', 'csv'])

T2ttoff.addMassPlane(T2tt_1)
T2bbffff.addMassPlane(T2tt_1)




databaseCreator.create()


