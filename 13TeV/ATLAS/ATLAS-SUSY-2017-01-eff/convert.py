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
info 				= MetaInfoInput('ATLAS-SUSY-2017-01')
info.url 			= 'https://atlas.web.cern.ch/Atlas/GROUPS/PHYSICS/PAPERS/SUSY-2017-01/'
info.sqrts 			= 13
info.lumi 			= 36.1
info.prettyName 	= 'EWK WH(bb) + Etmiss'
info.private 		= False
info.arxiv 			= 'arXiv:1812.09432'
info.contact 		= 'atlas-phys-susy-conveners@cern.ch'
info.publication 	= 'https://journals.aps.org/prd/abstract/10.1103/PhysRevD.100.012006'
info.comment	 	= 'Each SR has their observed and expected upper limits in datainfo.txt replaced with the upper limits provided by ATLAS (arxiv.org/pdf/1812.09432 Table 17).'
	


SR   = {'obsN' 	: [7, 1, 5, 7, 6], 
		'expN'  : [8, 2.5, 4.6, 2.8, 5.7], 
		'bgErr' : [4, 1.3, 1.2, 1, 2.3], 
		'SR' 	: ['SRHad-Low', 'SRHad-High', 'SR1Lbb-High', 'SR1Lbb-Medium', 'SR1Lbb-Low']}

DATA = {'expExcl' : ['orig/Expectedlimit0lbb.csv', 'orig/Expectedlimit0lbb.csv', 'orig/Expectedlimit1lbb.csv', 'orig/Expectedlimit1lbb.csv', 'orig/Expectedlimit1lbb.csv'],
		'obsExcl' : ['orig/Observedlimit0lbb.csv', 'orig/Observedlimit0lbb.csv', 'orig/Observedlimit1lbb.csv', 'orig/Observedlimit1lbb.csv', 'orig/Observedlimit1lbb.csv'], 
		'effMap'  : ['orig/EffMap_TChiWH_SRHad-Low.txt', 'orig/EffMap_TChiWH_SRHad-High.txt', 'orig/EffMap_TChiWH_SR1Lbb-High.txt', 'orig/EffMap_TChiWH_SR1Lbb-Medium.txt', 'orig/EffMap_TChiWH_SR1Lbb-Low.txt'],
		'fig'	  : ['Aux Fig. 12a,b','Aux Fig. 12c,d','Aux Fig. 15e,f','Aux Fig. 15c,d','Aux Fig. 15a,b'],
		'figUrl'  : ['12b.png', '12d.png', '15f.png', '15d.png', '15b.png'],
		'dataUrl' : ['Efficiency 0lbb SRHad-Low','Efficiency 0lbb SRHad-High','Efficiency 1lbb SR1Lbb-High','Efficiency 1lbb SR1Lbb-Medium','Efficiency 1lbb SR1Lbb-Low']}

for i in range(len(SR['obsN'])):
	#+++++++ dataset block ++++++++++++++
	dataset = DataSetInput(SR['SR'][i])
	dataset.setInfo(dataType = 'efficiencyMap', dataId = SR['SR'][i], observedN = SR['obsN'][i], expectedBG = SR['expN'][i], bgError = SR['bgErr'][i])
	#+++++++ next txName block ++++++++++++++
	newTx 							= dataset.addTxName('TChiWH')
	newTx.checked 					= 'No'
	newTx.constraint 				= "[[['W']],[['higgs']]]"
	newTx.conditionDescription 		= None
	newTx.condition 				= None
	newTx.source 					= 'ATLAS'
	#+++++++ next mass plane block ++++++++++++++
	newPlane 						= newTx.addMassPlane(2*[[x, y]])
	newPlane.figure 				= DATA['fig'][i]
	newPlane.figureUrl 				= 'https://atlas.web.cern.ch/Atlas/GROUPS/PHYSICS/PAPERS/SUSY-2017-01/fig_' + DATA['figUrl'][i]
	newPlane.dataUrl 				= 'https://www.hepdata.net/record/ins1711261?version=1&table=' + DATA['dataUrl'][i]
	newPlane.setSources(dataLabels 	= ['expExclusion', 'obsExclusion', 'efficiencyMap'],
					dataFiles 		= [DATA['expExcl'][i], DATA['obsExcl'][i], DATA['effMap'][i]],
					dataFormats		= ['csv', 'csv', 'txt'])

databaseCreator.create()


### manually replace ULs in dataInfo.txt
### taken from https://arxiv.org/pdf/1812.09432.pdf Table 17

from shutil import move
from os import remove

UL_obs = {'SRHad-Low' 		: 0.26, 
		  'SRHad-High' 		: 0.10, 
		  'SR1Lbb-High' 	: 0.18, 
		  'SR1Lbb-Medium' 	: 0.28,
		  'SR1Lbb-Low' 		: 0.23}

UL_exp = {'SRHad-Low' 		: 9.5, 
		  'SRHad-High' 		: 4.3, 
		  'SR1Lbb-High' 	: 6.1, 
		  'SR1Lbb-Medium' 	: 5.6,
		  'SR1Lbb-Low' 		: 8.0}

for sr, ul in UL_exp.items(): UL_exp[sr] = round(ul / 36.1, 3)

name = 'dataInfo.txt'
phrase_obs = 'upperLimit: '
phrase_exp = 'expectedUpperLimit: '
pos_obs = 5
pos_exp	= 6

for sr in SR['SR']:

	path = sr + '/' + name
	
	with open(path) as file:
		out = file.read()
		lines = out.split('\n')
		old_obs = lines[pos_obs].split(phrase_obs)[1]
		old_exp = lines[pos_exp].split(phrase_exp)[1]
		lines[pos_obs] = '{}{}*fb\t# official number taken from ATLAS (calculated value: {})'.format(phrase_obs, UL_obs[sr], old_obs)
		lines[pos_exp] = '{}{}*fb\t# official number taken from ATLAS (calculated value: {})'.format(phrase_exp, UL_exp[sr], old_exp)
		
		with open(name, 'w') as temp:
			for line in lines:
				temp.write(line + '\n')

	remove(path)
	move(name, path)

