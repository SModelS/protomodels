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
info 				= MetaInfoInput('ATLAS-SUSY-2018-32')
info.url 			= 'https://atlas.web.cern.ch/Atlas/GROUPS/PHYSICS/PAPERS/SUSY-2018-32/'
info.sqrts 			= 13
info.lumi 			= 139
info.prettyName 	= '2 OS leptons + Etmiss'
info.private 		= False
info.arxiv 			= 'arXiv:1908.08215'
info.contact 		= 'atlas-phys-susy-conveners@cern.ch'
info.publication 	= 'https://epjc.epj.org/articles/epjc/abs/2020/02/10052_2019_Article_7594/10052_2019_Article_7594.html'


TChiWW = {
'name' 		 : 'TChiWW',
'info' 		 :{'figure' 		: 'Fig.7a', 
			   'figureUrl' 		: 'https://atlas.web.cern.ch/Atlas/GROUPS/PHYSICS/PAPERS/SUSY-2018-32/figaux_03a.png', 
			   'dataUrl' 		: 'https://www.hepdata.net/record/ins1750597?version=1&table=xsec upper limits 1'},
'sources'	 :{'expExcl'		: 'orig/HEPData-ins1750597-v1-Exclusion_contour_(exp)_1.csv',
			   'obsExcl'		: 'orig/HEPData-ins1750597-v1-Exclusion_contour_(obs)_1.csv',
			   'upLimit'		: 'orig/HEPData-ins1750597-v1-xsec_upper_limits_1.csv'},
'constraint' : "[[['W']],[['W']]]",
'condDesc'	 : None,
'condition'	 : None,
'massPlane'  : 2*[[x, y]]}


TChipChimSlepSlep = {
'name' 		 : 'TChipChimSlepSlep',
'valTarball' : 'TChipChimSlepSlepAll.tar.gz',
'info' 		 :{'figure' 		: 'Fig.7b', 
			   'figureUrl' 		: 'https://atlas.web.cern.ch/Atlas/GROUPS/PHYSICS/PAPERS/SUSY-2018-32/fig_07b.png', 
			   'dataUrl' 		: 'https://www.hepdata.net/record/ins1750597?version=1&table=xsec upper limits 2'},
'sources'	 :{'expExcl'		: 'orig/HEPData-ins1750597-v1-Exclusion_contour_(exp)_2.csv',
			   'obsExcl'		: 'orig/HEPData-ins1750597-v1-Exclusion_contour_(exp)_2.csv',
			   'upLimit'		: 'orig/HEPData-ins1750597-v1-xsec_upper_limits_2.csv'},
'constraint' : "2.25*([[['mu+'],['nu']],[['mu-'],['nu']]] + [[['nu'],['mu+']],[['nu'],['mu-']]] + [[['mu+'],['nu']],[['nu'],['mu-']]] + [[['e+'],['nu']],[['e-'],['nu']]] + [[['nu'],['e+']],[['nu'],['e-']]] + [[['e+'],['nu']],[['nu'],['e-']]] + [[['nu'],['mu+']],[['nu'],['e-']]] + [[['nu'],['mu+']],[['e-'],['nu']]] + [[['mu+'],['nu']],[['e-'],['nu']]] + [[['mu+'],['nu']],[['nu'],['e-']]] + [[['nu'],['mu-']],[['nu'],['e+']]] + [[['nu'],['mu-']],[['e+'],['nu']]] + [[['mu-'],['nu']],[['e+'],['nu']]] + [[['mu-'],['nu']],[['nu'],['e+']]])",
'condDesc'   : "[[['mu+'],['nu']],[['mu-'],['nu']]] + [[['nu'],['mu+']],[['nu'],['mu-']]] + [[['mu+'],['nu']],[['nu'],['mu-']]] >= [[['e+'],['nu']],[['e-'],['nu']]] + [[['nu'],['e+']],[['nu'],['e-']]] + [[['e+'],['nu']],[['nu'],['e-']]], 2*([[['mu+'],['nu']],[['mu-'],['nu']]] + [[['nu'],['mu+']],[['nu'],['mu-']]] + [[['mu+'],['nu']],[['nu'],['mu-']]]) >= [[['nu'],['mu+']],[['nu'],['e-']]] + [[['nu'],['mu+']],[['e-'],['nu']]] + [[['mu+'],['nu']],[['e-'],['nu']]] + [[['mu+'],['nu']],[['nu'],['e-']]] + [[['nu'],['mu-']],[['nu'],['e+']]] + [[['nu'],['mu-']],[['e+'],['nu']]] + [[['mu-'],['nu']],[['e+'],['nu']]] + [[['mu-'],['nu']],[['nu'],['e+']]]",
'condition'  : "Cgtr([[['mu+'],['nu']],[['mu-'],['nu']]] + [[['nu'],['mu+']],[['nu'],['mu-']]] + [[['mu+'],['nu']],[['nu'],['mu-']]], [[['e+'],['nu']],[['e-'],['nu']]] + [[['nu'],['e+']],[['nu'],['e-']]] + [[['e+'],['nu']],[['nu'],['e-']]]); Cgtr(2.*( [[['mu+'],['nu']],[['mu-'],['nu']]] + [[['nu'],['mu+']],[['nu'],['mu-']]] + [[['mu+'],['nu']],[['nu'],['mu-']]] ), [[['nu'],['mu+']],[['nu'],['e-']]] + [[['nu'],['mu+']],[['e-'],['nu']]] + [[['mu+'],['nu']],[['e-'],['nu']]] + [[['mu+'],['nu']],[['nu'],['e-']]] + [[['nu'],['mu-']],[['nu'],['e+']]] + [[['nu'],['mu-']],[['e+'],['nu']]] + [[['mu-'],['nu']],[['e+'],['nu']]] + [[['mu-'],['nu']],[['nu'],['e+']]])",
'massPlane'  : 2*[[x, 0.5*(x+y), y]]}

TSlepSlep = {
'name' 		 : 'TSlepSlep',
#'valTarball' : 'TSlepSlepAll.tar.gz',
'info' 		 :{'figure' 		: 'Fig.7c', 
			   'figureUrl' 		: 'https://atlas.web.cern.ch/Atlas/GROUPS/PHYSICS/PAPERS/SUSY-2018-32/fig_07c.png', 
			   'dataUrl' 		: 'https://www.hepdata.net/record/ins1750597?version=1&table=xsec upper limits 3'},
'sources'	 :{'expExcl'		: 'orig/HEPData-ins1750597-v1-Exclusion_contour_(exp)_3.csv',
			   'obsExcl'		: 'orig/HEPData-ins1750597-v1-Exclusion_contour_(obs)_3.csv',
			   'upLimit'		: 'orig/HEPData-ins1750597-v1-xsec_upper_limits_3.csv'},
'constraint' : "[[['e+']],[['e-']]] + [[['mu+']],[['mu-']]]",
'condDesc'	 : "[[['mu+']],[['mu-']]] > [[['e+']],[['e-']]]",
'condition'	 : "Cgtr([[['mu+']],[['mu-']]], [[['e+']],[['e-']]])",
'massPlane'  : 2*[[x, y]]}

DATA = [TChiWW, TChipChimSlepSlep, TSlepSlep]

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
                 	dataFormats 	= ['csv', 'csv', 'csv'])

databaseCreator.create()
