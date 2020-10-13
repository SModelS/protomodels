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
info 			 = MetaInfoInput('ATLAS-SUSY-2016-16')
info.url 		 = 'https://atlas.web.cern.ch/Atlas/GROUPS/PHYSICS/PAPERS/SUSY-2016-16/'
info.sqrts 		 = 13
info.lumi 		 = 36.1
info.prettyName  = '1L stop'
info.private 	 = False
info.arxiv 		 = 'https://arxiv.org/abs/1711.11520'
info.contact 	 = 'atlas-phys-susy-conveners@cern.ch'
info.publication = 'JHEP 06 (2018) 108'

T2tt = {
'name' 		 : ['T2tt','T2ttoff','T2bbffff'],
'info' 		 :{'figure' 		: ['Fig.20', 'Fig.21'],
			   'figureUrl' 		: ['https://atlas.web.cern.ch/Atlas/GROUPS/PHYSICS/PAPERS/SUSY-2016-16/fig_20.png', 'https://atlas.web.cern.ch/Atlas/GROUPS/PHYSICS/PAPERS/SUSY-2016-16/fig_21.png'],
			   'dataUrl' 		: ['https://www.hepdata.net/record/ins1639856?version=4&table=Table60','https://www.hepdata.net/record/ins1639856?version=4&table=Table61']},
'sources'	 :{'expExcl'		: ['orig/HEPData-ins1639856-v4-Table_16.csv','orig/HEPData-ins1639856-v4-Table_19.csv'],
			   'obsExcl'		: ['orig/HEPData-ins1639856-v4-Table_17.csv','orig/HEPData-ins1639856-v4-Table_20.csv'],
			   'upLimit'		: ['orig/HEPData-ins1639856-v4-Table_60.csv','orig/HEPData-ins1639856-v4-Table_61.csv']},
'constraint' : ['[[[t]],[[t]]]','[[[b, W]],[[b, W]]]','27./8.*[[[b, l, nu]],[[b, jet, jet]]]'],
'massPlane'  : [2*[[x, y]],2*[[x, y]]]}

T6bbWW = {
'name' 		 : 'T6bbWW',
'info' 		 :{'figure' 		: 'Fig.23', 
			   'figureUrl' 		: 'https://atlas.web.cern.ch/Atlas/GROUPS/PHYSICS/PAPERS/SUSY-2016-16/fig_23.png', 
			   'dataUrl' 		: 'https://www.hepdata.net/record/ins1639856?version=4&table=Table70'},
'sources'	 :{'expExcl'		: 'orig/HEPData-ins1639856-v4-Table_28.csv',
			   'obsExcl'		: 'orig/HEPData-ins1639856-v4-Table_29.csv',
			   'upLimit'		: 'orig/HEPData-ins1639856-v4-Table_70.csv'},
'constraint' : '[[[b],[W]],[[b],[W]]]',
'massPlane'  : 2*[[x, x - 10, y]]}


DATA = [T2tt]

#+++++++ dataset block ++++++++++++++
dataset = DataSetInput('data')
dataset.setInfo(dataType = 'upperLimit', dataId = None)


for TX in DATA:
	#+++++++ next txName block ++++++++++++++
	newTx 							= dataset.addTxName(TX['name'][0])
	newTx.checked 					= 'False'
	newTx.constraint 				= TX['constraint'][0]
	newTx.conditionDescription 		= None
	newTx.condition 				= None
	newTx.source 					= 'ATLAS'
	#+++++++ next txName block ++++++++++++++
	newTxOff1 						= dataset.addTxName(TX['name'][1])
	newTxOff1.checked 				= 'False'
	newTxOff1.constraint 			= TX['constraint'][1]
	newTxOff1.conditionDescription 	= None
	newTxOff1.condition 			= None
	newTxOff1.massConstraint 		= [['80 <= dm < 169.0'], ['80 <= dm < 169.0']]
	newTxOff1.source 				= 'ATLAS'
	#+++++++ next mass plane block ++++++++++++++
	newTxOff2						= dataset.addTxName(TX['name'][2])
	newTxOff2.checked				= 'False'
	newTxOff2.constraint			= TX['constraint'][2]
	newTxOff2.conditionDescription 	= None
	newTxOff2.condition				= None
	newTxOff2.source				= 'ATLAS'
	newTxOff2.massConstraint		= [['dm < 80'], ['dm < 80']]

	#for i in range(len(TX['info']['figure'])):
	i = 0

	#+++++++ next mass plane block ++++++++++++++
	newPlane 						= newTx.addMassPlane(TX['massPlane'][i])
	newPlane.figure 				= TX['info']['figure'][i]
	newPlane.figureUrl 				= TX['info']['figureUrl'][i]
	newPlane.dataUrl 				= TX['info']['dataUrl'][i]
	newPlane.setSources(dataLabels 	= ['expExclusion', 'obsExclusion', 'upperLimits'],
					dataFiles 		= [TX['sources']['expExcl'][i], TX['sources']['obsExcl'][i], TX['sources']['upLimit'][i]],
					units			= [ None, None, 'pb' ],
				 	coordinates 	= [ {x: 0, y: 1, 'value': None}, {x: 0, y: 1, 'value': None},  {x : 1, y: 0, 'value' :2} ],
	             	dataFormats 	= ['csv', 'csv', 'csv'])

	newTxOff1.addMassPlane(newPlane)
	newTxOff2.addMassPlane(newPlane)

DATA = [T6bbWW]

for TX in DATA:
	#+++++++ next txName block ++++++++++++++
	newTx 							= dataset.addTxName(TX['name'])
	newTx.checked 					= 'False'
	newTx.constraint 				= TX['constraint']
	newTx.conditionDescription 		= None
	newTx.condition 				= None
	newTx.source 					= 'ATLAS'
	#+++++++ next mass plane block ++++++++++++++
	newPlane 						= newTx.addMassPlane(TX['massPlane'])
	newPlane.figure 				= TX['info']['figure']
	newPlane.figureUrl 				= TX['info']['figureUrl']
	newPlane.dataUrl 				= TX['info']['dataUrl']
	newPlane.setSources(dataLabels 	= ['expExclusion', 'obsExclusion', 'upperLimits'],
					dataFiles 		= [TX['sources']['expExcl'], TX['sources']['obsExcl'], TX['sources']['upLimit']],
					units			= [ None, None, 'pb' ],
				 	coordinates 	= [ {x: 0, y: 1, 'value': None}, {x: 0, y: 1, 'value': None},  {x : 1, y: 0, 'value' :2} ],
                 	dataFormats 	= ['csv', 'csv', 'csv'])

databaseCreator.create()
