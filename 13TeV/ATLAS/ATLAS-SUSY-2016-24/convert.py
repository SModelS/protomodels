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
info = MetaInfoInput('ATLAS-SUSY-2016-24')
info.url = 'https://atlas.web.cern.ch/Atlas/GROUPS/PHYSICS/PAPERS/SUSY-2016-24/'
info.sqrts = 13
info.lumi = 36.1
info.prettyName = '2-3 leptons + Etmiss, EWino'
info.private = False
info.arxiv =  'https://arxiv.org/abs/1803.02762'
info.contact = 'atlas-phys-susy-conveners@cern.ch'
info.publication = 'https://link.springer.com/article/10.1140/epjc/s10052-018-6423-7'

TChipChimSlepSlep = {
'name' 		 : 'TChipChimSlepSlep',
'valTarball' : 'TChipChimSlepSlepAll.tar.gz',
'info' 		 :{'figure' 		: 'figaux_31a', 
			   'figureUrl' 		: 'https://atlas.web.cern.ch/Atlas/GROUPS/PHYSICS/PAPERS/SUSY-2016-24/figaux_31a.png', 
			   'dataUrl' 		: 'https://www.hepdata.net/record/ins1658902?version=1&table=Table78'},
'sources'	 :{'expExcl'		: 'orig/HEPData-ins1658902-v1-Table_13.csv',
			   'obsExcl'		: 'orig/HEPData-ins1658902-v1-Table_14.csv',
			   'upLimit'		: 'orig/HEPData-ins1658902-v1-Table_78.csv'},
'constraint' : "2.25*([[['mu+'],['nu']],[['mu-'],['nu']]] + [[['nu'],['mu+']],[['nu'],['mu-']]] + [[['mu+'],['nu']],[['nu'],['mu-']]] + [[['e+'],['nu']],[['e-'],['nu']]] + [[['nu'],['e+']],[['nu'],['e-']]] + [[['e+'],['nu']],[['nu'],['e-']]] + [[['nu'],['mu+']],[['nu'],['e-']]] + [[['nu'],['mu+']],[['e-'],['nu']]] + [[['mu+'],['nu']],[['e-'],['nu']]] + [[['mu+'],['nu']],[['nu'],['e-']]] + [[['nu'],['mu-']],[['nu'],['e+']]] + [[['nu'],['mu-']],[['e+'],['nu']]] + [[['mu-'],['nu']],[['e+'],['nu']]] + [[['mu-'],['nu']],[['nu'],['e+']]])",
'condDesc'   : "[[['mu+'],['nu']],[['mu-'],['nu']]] + [[['nu'],['mu+']],[['nu'],['mu-']]] + [[['mu+'],['nu']],[['nu'],['mu-']]] >= [[['e+'],['nu']],[['e-'],['nu']]] + [[['nu'],['e+']],[['nu'],['e-']]] + [[['e+'],['nu']],[['nu'],['e-']]], 2*([[['mu+'],['nu']],[['mu-'],['nu']]] + [[['nu'],['mu+']],[['nu'],['mu-']]] + [[['mu+'],['nu']],[['nu'],['mu-']]]) >= [[['nu'],['mu+']],[['nu'],['e-']]] + [[['nu'],['mu+']],[['e-'],['nu']]] + [[['mu+'],['nu']],[['e-'],['nu']]] + [[['mu+'],['nu']],[['nu'],['e-']]] + [[['nu'],['mu-']],[['nu'],['e+']]] + [[['nu'],['mu-']],[['e+'],['nu']]] + [[['mu-'],['nu']],[['e+'],['nu']]] + [[['mu-'],['nu']],[['nu'],['e+']]]",
'condition'  : "Cgtr([[['mu+'],['nu']],[['mu-'],['nu']]] + [[['nu'],['mu+']],[['nu'],['mu-']]] + [[['mu+'],['nu']],[['nu'],['mu-']]], [[['e+'],['nu']],[['e-'],['nu']]] + [[['nu'],['e+']],[['nu'],['e-']]] + [[['e+'],['nu']],[['nu'],['e-']]]); Cgtr(2.*( [[['mu+'],['nu']],[['mu-'],['nu']]] + [[['nu'],['mu+']],[['nu'],['mu-']]] + [[['mu+'],['nu']],[['nu'],['mu-']]] ), [[['nu'],['mu+']],[['nu'],['e-']]] + [[['nu'],['mu+']],[['e-'],['nu']]] + [[['mu+'],['nu']],[['e-'],['nu']]] + [[['mu+'],['nu']],[['nu'],['e-']]] + [[['nu'],['mu-']],[['nu'],['e+']]] + [[['nu'],['mu-']],[['e+'],['nu']]] + [[['mu-'],['nu']],[['e+'],['nu']]] + [[['mu-'],['nu']],[['nu'],['e+']]])",
'massPlane'  : 2*[[x, 0.5*(x+y), y]]}

TSlepSlep = {
'name' 		 : 'TSlepSlep',
'valTarball' : 'TSlepSlepAll.tar.gz',
'info' 		 :{'figure' 		: 'figaux_31b', 
			   'figureUrl' 		: 'https://atlas.web.cern.ch/Atlas/GROUPS/PHYSICS/PAPERS/SUSY-2016-24/figaux_31b.png', 
			   'dataUrl' 		: 'https://www.hepdata.net/record/ins1658902?version=1&table=Table79'},
'sources'	 :{'expExcl'		: 'orig/HEPData-ins1658902-v1-Table_15.csv',
			   'obsExcl'		: 'orig/HEPData-ins1658902-v1-Table_16.csv',
			   'upLimit'		: 'orig/HEPData-ins1658902-v1-Table_79.csv'},
'constraint' : "2.25*([[['e+']],[['e-']]] + [[['mu+']],[['mu-']]])",
'condDesc'	 : "[[['mu+']],[['mu-']]] > [[['e+']],[['e-']]]",
'condition'	 : "Cgtr([[['mu+']],[['mu-']]], [[['e+']],[['e-']]])",
'massPlane'  : 2*[[x, y]]}

TChiChipmSlepSlep = {
'name' 		 : 'TChiChipmSlepSlep',
'valTarball' : 'TChiChipmSlepLNoTau.tar.gz',
'info' 		 :{'figure' 		: 'figaux_32a',
			   'figureUrl' 		: 'https://atlas.web.cern.ch/Atlas/GROUPS/PHYSICS/PAPERS/SUSY-2016-24/figaux_32a.png', 
			   'dataUrl' 		: 'https://www.hepdata.net/record/ins1658902?version=1&table=Table80'},
'sources'	 :{'expExcl'		: 'orig/HEPData-ins1658902-v1-Table_17.csv',
			   'obsExcl'		: 'orig/HEPData-ins1658902-v1-Table_18.csv',
			   'upLimit'		: 'orig/HEPData-ins1658902-v1-Table_80.csv'},
'constraint' : "2.25*([[['e+'],['e-']],[['l'],['nu']]] + [[['e-'],['e+']],[['l'],['nu']]] + [[['e+'],['e-']],[['nu'],['l']]] + [[['e-'],['e+']],[['nu'],['l']]] + [[['mu+'],['mu-']],[['l'],['nu']]] + [[['mu-'],['mu+']],[['l'],['nu']]] + [[['mu+'],['mu-']],[['nu'],['l']]] + [[['mu-'],['mu+']],[['nu'],['l']]])",
'condDesc'	 : "[[['mu+'],['mu-']],[['l'],['nu']]] + [[['mu-'],['mu+']],[['l'],['nu']]] + [[['mu+'],['mu-']],[['nu'],['l']]] + [[['mu-'],['mu+']],[['nu'],['l']]] >= [[['e+'],['e-']],[['l'],['nu']]] + [[['e-'],['e+']],[['l'],['nu']]] + [[['e+'],['e-']],[['nu'],['l']]] + [[['e-'],['e+']],[['nu'],['l']]]",
'condition'	 : "Cgtr([[['mu+'],['mu-']],[['l'],['nu']]] + [[['mu-'],['mu+']],[['l'],['nu']]] + [[['mu+'],['mu-']],[['nu'],['l']]] + [[['mu-'],['mu+']],[['nu'],['l']]], [[['e+'],['e-']],[['l'],['nu']]] + [[['e-'],['e+']],[['l'],['nu']]] + [[['e+'],['e-']],[['nu'],['l']]] + [[['e-'],['e+']],[['nu'],['l']]])",
'massPlane'  : 2*[[x, 0.5*(x+y), y]]}

TChiWZ = {
'name' 		 : 'TChiWZ',
'info' 		 :{'figure' 		: 'figaux_32b', 
			   'figureUrl' 		: 'https://atlas.web.cern.ch/Atlas/GROUPS/PHYSICS/PAPERS/SUSY-2016-24/figaux_32b.png', 
			   'dataUrl' 		: 'https://www.hepdata.net/record/ins1658902?version=1&table=Table81'},
'sources'	 :{'expExcl'		: 'orig/HEPData-ins1658902-v1-Table_19.csv',
			   'obsExcl'		: 'orig/HEPData-ins1658902-v1-Table_20.csv',
			   'upLimit'		: 'orig/HEPData-ins1658902-v1-Table_81.csv'},
'constraint' : "[[['W']],[['Z']]]",
'condDesc'	 :	None,
'condition'	 :	None,
'massPlane'  : 2*[[x, y]]}

DATA = [TChipChimSlepSlep, TSlepSlep, TChiChipmSlepSlep, TChiWZ]

#+++++++ dataset block ++++++++++++++
dataset = DataSetInput('data')
dataset.setInfo(dataType = 'upperLimit', dataId = None)

for TX in DATA:
	#+++++++ next txName block ++++++++++++++
	newTx 							= dataset.addTxName(TX['name'])
	newTx.checked 					= 'False'
	newTx.constraint 				= TX['constraint']
	newTx.conditionDescription 		= TX['condDesc']
	newTx.condition 				= TX['condition']
	newTx.source 					= 'ATLAS'
	if 'valTarball' in TX:
		print('validating ' + TX['valTarball'] + ' with ' + TX['name'])
		newTx.validationTarball = TX['valTarball']
	#+++++++ next mass plane block ++++++++++++++
	newPlane 						= newTx.addMassPlane(TX['massPlane'])
	newPlane.figure 				= TX['info']['figure']
	newPlane.figureUrl 				= TX['info']['figureUrl']
	newPlane.dataUrl 				= TX['info']['dataUrl']
	newPlane.setSources(dataLabels 	= ['expExclusion', 'obsExclusion', 'upperLimits'],
					dataFiles 		= [TX['sources']['expExcl'], TX['sources']['obsExcl'], TX['sources']['upLimit']],
					units			= [ None, None, 'fb' ],
				 	coordinates 	= [ {x: 0, y: 1, 'value': None}, {x: 0, y: 1, 'value': None},  {x : 1, y: 0, 'value' :2} ],
                 	dataFormats 	= ['csv', 'csv', 'csv'])


databaseCreator.create()

