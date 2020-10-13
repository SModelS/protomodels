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
info 			= MetaInfoInput('ATLAS-SUSY-2018-16')
info.url 		= 'https://atlas.web.cern.ch/Atlas/GROUPS/PHYSICS/PAPERS/SUSY-2018-16/'
info.sqrts 		= 13
info.lumi 		= 139
info.prettyName 	= 'EW production in models with compressed mass spectra'
info.private 		= False
info.arxiv 		= 'arXiv:1911.12606'
info.contact            = 'atlas-phys-susy-conveners@cern.ch'

SR   = {'obsN' 	: [472,131,341], 
		'expN'  : [443,137,331], 
		'bgErr' : [31,6,11], 
		'SR' 	: ['SR-S','SR-S-high','SR-S-low']}

#+++++++ next txName block ++++++++++++++
TSlepSlep = {
'name' 		 : 'TSlepSlep',
'info' 		 :{'figure' 		: 'Fig.16a', 
			   'figureUrl' 		: 'https://atlas.web.cern.ch/Atlas/GROUPS/PHYSICS/PAPERS/SUSY-2018-16/fig_16a.png', 
			   'dataUrl' 		: 'https://www.hepdata.net/record/ins1767649?version=1&table=Figure 16a Observed'},
'sources'	 :{'expExcl'		: 'orig/Figure16aExpected.csv',
			   'obsExcl'		: 'orig/Figure16aObserved.csv',
			   'effMap'			: 'orig/EffMap_TSlepSlep_'},
'constraint' : "[[['mu+']],[['mu-']]]+[[['e+']],[['e-']]]",
'massConstr' : None,
'massPlane'  : 2*[[x, (x-y)]]}

TX = TSlepSlep
for i in range(len(SR['SR'])):
	#+++++++ dataset block ++++++++++++++
	dataset = DataSetInput(SR['SR'][i])
	dataset.setInfo(dataType = 'efficiencyMap', dataId = SR['SR'][i], observedN = SR['obsN'][i], expectedBG = SR['expN'][i], bgError = SR['bgErr'][i])
	#+++++++ next txName block ++++++++++++++
	newTx							= dataset.addTxName(TX['name'])
	newTx.checked					= 'no'
	newTx.constraint				= TX['constraint']
	newTx.conditionDescription 		= None
	newTx.condition					= None
	newTx.source					= 'ATLAS'
	newTx.massConstraint			= TX['massConstr']
	#+++++++ next mass plane block ++++++++++++++
	newPlane 						= newTx.addMassPlane(TX['massPlane'])
	newPlane.figure 				= TX['info']['figure']
	newPlane.figureUrl 				= TX['info']['figureUrl']
	newPlane.dataUrl 				= TX['info']['dataUrl']
	newPlane.setSources(dataLabels 	= ['expExclusion', 'obsExclusion', 'efficiencyMap'],
                        dataFiles 	= [TX['sources']['expExcl'], TX['sources']['obsExcl'], TX['sources']['effMap'] + SR['SR'][i] + '.txt'],
                        dataFormats 	= ['csv', 'csv', 'txt'])


databaseCreator.create()
