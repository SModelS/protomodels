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
info.prettyName = '2+ leptons (e,mu) + jets + Etmiss'
info.private = False
info.arxiv =  'https://arxiv.org/abs/1803.02762'
info.contact = 'atlas-phys-susy-conveners@cern.ch'
info.publication = 'https://link.springer.com/article/10.1140/epjc/s10052-018-6423-7'

SF   = {'obsN' 	: [153, 9], 
		'expN'  : [133, 9.8], 
		'bgErr' : [22, 2.9], 
		'SR' 	: ['SR2-SF-loose', 'SR2-SF-tight'],
		'Fig'	: ['22a; 22b','22c; 22d'],
		'URL'	: [[['25','26'],['27','28']],
				  [['37','38'],['39','40']]]}
DF   = {'obsN' 	: [78, 11, 6, 2], 
		'expN'  : [68, 11.5, 2.1, 0.6], 
		'bgErr' : [7, 3.1, 1.9, 0.6], 
		'SR' 	: ['SR2-DF-100', 'SR2-DF-150', 'SR2-DF-200', 'SR2-DF-300'],
		'Fig'	: ['15a; 15b','15c; 15d','16a; 16b','16c; 16d'],
		'URL'	: [['29','30'],['31','32'],['33','34'],['35','36']]}
SR2  = {'obsN' 	: [11, 2, 0], 
		'expN'  : [4.2, 4.1, 1.6], 
		'bgErr' : [3.4, 2.6, 1.6], 
		'SR' 	: ['SR2-low', 'SR2-int', 'SR2-high'],
		'Fig'	: ['23a; 23b','23c; 23d','23e; 23f'],
		'URL'	: [['41','42'],['43','44'],['45','46']]}
WZ   = {'obsN' 	: [21, 1, 2, 1, 3, 4], 
		'expN'  : [21.7, 2.7, 1.6, 2.2, 1.8, 1.3], 
		'bgErr' : [2.9, 0.5, 0.3, 0.5, 0.3, 0.3], 
		'Fig'	: ['26a; 26b','26c; 26d','26e; 26f','27a; 27b','27c; 27d','27e; 27f'],
		'SR' 	: ['WZ-0Ja', 'WZ-0Jb', 'WZ-0Jc', 'WZ-1Ja', 'WZ-1Jb', 'WZ-1Jc'],
		'URL'	: [['57','58'],['59','60'],['61','62'],['63','64'],['65','66'],['67','68']]}
SLEP = {'obsN' 	: [4, 3, 9, 0, 0], 
		'expN'  : [2.2, 2.8, 5.4, 1.4, 1.1], 
		'bgErr' : [0.8, 0.4, 0.9, 0.4, 0.2], 
		'SR' 	: ['slep-a', 'slep-b', 'slep-c', 'slep-d', 'slep-e'],
		'Fig'	: ['24a; 24b','24c; 24d','24e; 24f','25a; 25b','25c; 25d'],
		'URL'	: [['47','48'],['49','50'],['51','52'],['53','54'],['55','56']]}

TChipChimSlepSlep = {
'name' 		 : 'TChipChimSlepSlep',
'valTarball' : 'TChipChimSlepSlepAll.tar.gz',
'sources'	 :{'expExcl': 'orig/HEPData-ins1658902-v1-Table_13.csv',
			   'obsExcl': 'orig/HEPData-ins1658902-v1-Table_14.csv'},
'constraint' : "[[['mu+'],['nu']],[['mu-'],['nu']]] + [[['nu'],['mu+']],[['nu'],['mu-']]] + [[['mu+'],['nu']],[['nu'],['mu-']]] + [[['e+'],['nu']],[['e-'],['nu']]] + [[['nu'],['e+']],[['nu'],['e-']]] + [[['e+'],['nu']],[['nu'],['e-']]] + [[['nu'],['mu+']],[['nu'],['e-']]] + [[['nu'],['mu+']],[['e-'],['nu']]] + [[['mu+'],['nu']],[['e-'],['nu']]] + [[['mu+'],['nu']],[['nu'],['e-']]] + [[['nu'],['mu-']],[['nu'],['e+']]] + [[['nu'],['mu-']],[['e+'],['nu']]] + [[['mu-'],['nu']],[['e+'],['nu']]] + [[['mu-'],['nu']],[['nu'],['e+']]]",
'condDesc'   : "[[['mu+'],['nu']],[['mu-'],['nu']]] + [[['nu'],['mu+']],[['nu'],['mu-']]] + [[['mu+'],['nu']],[['nu'],['mu-']]] >= [[['e+'],['nu']],[['e-'],['nu']]] + [[['nu'],['e+']],[['nu'],['e-']]] + [[['e+'],['nu']],[['nu'],['e-']]], 2*([[['mu+'],['nu']],[['mu-'],['nu']]] + [[['nu'],['mu+']],[['nu'],['mu-']]] + [[['mu+'],['nu']],[['nu'],['mu-']]]) >= [[['nu'],['mu+']],[['nu'],['e-']]] + [[['nu'],['mu+']],[['e-'],['nu']]] + [[['mu+'],['nu']],[['e-'],['nu']]] + [[['mu+'],['nu']],[['nu'],['e-']]] + [[['nu'],['mu-']],[['nu'],['e+']]] + [[['nu'],['mu-']],[['e+'],['nu']]] + [[['mu-'],['nu']],[['e+'],['nu']]] + [[['mu-'],['nu']],[['nu'],['e+']]]",
'condition'  : "Cgtr([[['mu+'],['nu']],[['mu-'],['nu']]] + [[['nu'],['mu+']],[['nu'],['mu-']]] + [[['mu+'],['nu']],[['nu'],['mu-']]], [[['e+'],['nu']],[['e-'],['nu']]] + [[['nu'],['e+']],[['nu'],['e-']]] + [[['e+'],['nu']],[['nu'],['e-']]]); Cgtr(2.*( [[['mu+'],['nu']],[['mu-'],['nu']]] + [[['nu'],['mu+']],[['nu'],['mu-']]] + [[['mu+'],['nu']],[['nu'],['mu-']]] ), [[['nu'],['mu+']],[['nu'],['e-']]] + [[['nu'],['mu+']],[['e-'],['nu']]] + [[['mu+'],['nu']],[['e-'],['nu']]] + [[['mu+'],['nu']],[['nu'],['e-']]] + [[['nu'],['mu-']],[['nu'],['e+']]] + [[['nu'],['mu-']],[['e+'],['nu']]] + [[['mu-'],['nu']],[['e+'],['nu']]] + [[['mu-'],['nu']],[['nu'],['e+']]])",
'massPlane'  : 2*[[x, 0.5*(x+y), y]],
'scales'	 : [1., 1., 2.25]}

TSlepSlep = {
'name' 		 : 'TSlepSlep',
'sources'	 :{'expExcl': 'orig/HEPData-ins1658902-v1-Table_15.csv',
			   'obsExcl': 'orig/HEPData-ins1658902-v1-Table_16.csv'},
'constraint' : "[[['e+']],[['e-']]]+[[['mu+']],[['mu-']]]",
'condDesc'	 : "[[['mu+']],[['mu-']]] > [[['e+']],[['e-']]]",
'condition'	 : "Cgtr([[['mu+']],[['mu-']]],[[['e+']],[['e-']]])",
'massPlane'  : 2*[[x, y]],
'scales'	 : [1., 1., 1.5]}

TSelSel = {
'name' 		 : 'TSelSel',
'sources'	 :{'expExcl': 'orig/HEPData-ins1658902-v1-Table_15.csv',
			   'obsExcl': 'orig/HEPData-ins1658902-v1-Table_16.csv'},
'constraint' : "[[['e+']],[['e-']]]",
'condDesc'	 : None,
'condition'	 : None,
'massPlane'  : 2*[[x, y]],
'scales'	 : [1., 1., 1.5]}

TSmuSmu = {
'name' 		 : 'TSmuSmu',
'sources'	 :{'expExcl': 'orig/HEPData-ins1658902-v1-Table_15.csv',
			   'obsExcl': 'orig/HEPData-ins1658902-v1-Table_16.csv'},
'constraint' : "[[['mu+']],[['mu-']]]",
'condDesc'	 : None,
'condition'	 : None,
'massPlane'  : 2*[[x, y]],
'scales'	 : [1., 1., 1.5]}

TChiChipmSlepSlep = {
'name' 		 : 'TChiChipmSlepSlep',
'valTarball' : 'TChiChipmSlepLNoTau.tar.gz',
'sources'	 :{'expExcl': 'orig/HEPData-ins1658902-v1-Table_17.csv',
			   'obsExcl': 'orig/HEPData-ins1658902-v1-Table_18.csv'},
'constraint' : "[[['e+'],['e-']],[['l'],['nu']]] + [[['e-'],['e+']],[['l'],['nu']]] + [[['e+'],['e-']],[['nu'],['l']]] + [[['e-'],['e+']],[['nu'],['l']]] + [[['mu+'],['mu-']],[['l'],['nu']]] + [[['mu-'],['mu+']],[['l'],['nu']]] + [[['mu+'],['mu-']],[['nu'],['l']]] + [[['mu-'],['mu+']],[['nu'],['l']]]",
'condDesc'	 : "[[['mu+'],['mu-']],[['l'],['nu']]] + [[['mu-'],['mu+']],[['l'],['nu']]] + [[['mu+'],['mu-']],[['nu'],['l']]] + [[['mu-'],['mu+']],[['nu'],['l']]] > [[['e+'],['e-']],[['l'],['nu']]] + [[['e-'],['e+']],[['l'],['nu']]] + [[['e+'],['e-']],[['nu'],['l']]] + [[['e-'],['e+']],[['nu'],['l']]]",
'condition'	 : "Cgtr([[['mu+'],['mu-']],[['l'],['nu']]] + [[['mu-'],['mu+']],[['l'],['nu']]] + [[['mu+'],['mu-']],[['nu'],['l']]] + [[['mu-'],['mu+']],[['nu'],['l']]], [[['e+'],['e-']],[['l'],['nu']]] + [[['e-'],['e+']],[['l'],['nu']]] + [[['e+'],['e-']],[['nu'],['l']]] + [[['e-'],['e+']],[['nu'],['l']]])",
'massPlane'  : 2*[[x, 0.5*(x+y), y]],
'scales'	 : [1., 1., 2.25]}

TChiWZ = {
'name' 		 : 'TChiWZ',
'sources'	 :{'expExcl': 'orig/HEPData-ins1658902-v1-Table_19.csv',
			   'obsExcl': 'orig/HEPData-ins1658902-v1-Table_20.csv'},
'constraint' : "[[['W']],[['Z']]]",
'condDesc'	 :	None,
'condition'	 :	None,
'massPlane'  : 2*[[x, y]],
'scales'	 : [1., 1., 1.]}

DATA = [([TChipChimSlepSlep, TSelSel, TSmuSmu], SF),
		(TChipChimSlepSlep, DF),
		(TChiChipmSlepSlep, SLEP),
		(TChiWZ, SR2),
		(TChiWZ, WZ),]

for TX, SR in DATA:
	for i in range(len(SR['obsN'])):
		#+++++++ dataset block ++++++++++++++
		dataset = DataSetInput(SR['SR'][i])
		dataset.setInfo(dataType = 'efficiencyMap', dataId = SR['SR'][i], observedN = SR['obsN'][i], expectedBG = SR['expN'][i], bgError = SR['bgErr'][i])
		if not isinstance(TX, list):
			#+++++++ next txName block ++++++++++++++
			newTx 							= dataset.addTxName(TX['name'])
			newTx.checked 					= 'No'
			newTx.constraint 				= TX['constraint']
			newTx.conditionDescription 		= TX['condDesc']
			newTx.condition 				= TX['condition']
			newTx.source 					= 'ATLAS'
			if 'valTarball' in TX:
				print('validating ' + TX['valTarball'] + ' with ' + TX['name'])
				newTx.validationTarball = TX['valTarball']
			#+++++++ next mass plane block ++++++++++++++
			newPlane 						= newTx.addMassPlane(TX['massPlane'])
			newPlane.figure 				= 'Aux. Fig ' + SR['Fig'][i]
			newPlane.figureUrl 				= 'https://atlas.web.cern.ch/Atlas/GROUPS/PHYSICS/PAPERS/SUSY-2016-24/figaux_' + SR['Fig'][i][0:3] + '.png' + '; ' + 'https://atlas.web.cern.ch/Atlas/GROUPS/PHYSICS/PAPERS/SUSY-2016-24/figaux_' + SR['Fig'][i][5:8] + '.png'
			newPlane.dataUrl 				= 'https://www.hepdata.net/record/ins1658902?version=1&table=Table' + SR['URL'][i][0] + '; ' + 'https://www.hepdata.net/record/ins1658902?version=1&table=Table' + SR['URL'][i][1]
			newPlane.setSources(dataLabels 	= ['expExclusion', 'obsExclusion', 'efficiencyMap'],
							dataFiles 		= [TX['sources']['expExcl'], TX['sources']['obsExcl'], 'orig/EffMap_' + TX['name'] + '_' + SR['SR'][i] + '.txt'],
							dataFormats		= ['csv', 'csv', 'txt'],
							scales		 	= TX['scales'])
		else:
			for n in range(len(TX)):
				#+++++++ next txName block ++++++++++++++
				newTx 							= dataset.addTxName(TX[n]['name'])
				newTx.checked 					= 'No'
				newTx.constraint 				= TX[n]['constraint']
				newTx.conditionDescription 		= TX[n]['condDesc']
				newTx.condition 				= TX[n]['condition']
				newTx.source 					= 'ATLAS'
				if 'valTarball' in TX[n]:
					print('using ' + TX[n]['valTarball'] + ' for ' + TX[n]['name'])
					newTx.validationTarball = TX[n]['valTarball']
				#+++++++ next mass plane block ++++++++++++++
				if TX[n]['name'] == 'TSelSel' or TX[n]['name'] == 'TSmuSmu': efftxname = 'TSlepSlep'
				else: efftxname = TX[n]['name']
				newPlane1 						= newTx.addMassPlane(TX[n]['massPlane'])
				newPlane1.figure 				= 'Aux. Fig ' + SR['Fig'][i]
				newPlane1.figureUrl 			= 'https://atlas.web.cern.ch/Atlas/GROUPS/PHYSICS/PAPERS/SUSY-2016-24/figaux_' + SR['Fig'][i][0:3] + '.png' + '; ' + 'https://atlas.web.cern.ch/Atlas/GROUPS/PHYSICS/PAPERS/SUSY-2016-24/figaux_' + SR['Fig'][i][5:8] + '.png'
				newPlane1.dataUrl 				= None#'https://www.hepdata.net/record/ins1658902?version=1&table=Table' + SR['URL'][n][i][0] + '; ' + 'https://www.hepdata.net/record/ins1658902?version=1&table=Table' + SR['URL'][n][i][1]
				newPlane1.setSources(dataLabels = ['expExclusion', 'obsExclusion', 'efficiencyMap'],
								dataFiles 		= [TX[n]['sources']['expExcl'], TX[n]['sources']['obsExcl'], 'orig/EffMap_' + efftxname + '_' + SR['SR'][i] + '.txt'],
								dataFormats		= ['csv', 'csv', 'txt'],
								scales		 	= TX[n]['scales'])

databaseCreator.create()





from shutil import move
from os import remove
### manually replace ULs in dataInfo.txt
### taken from https://arxiv.org/pdf/1803.02762.pdf table 15, p29

def replaceULs(UL_obs, UL_exp):

	
	name = 'dataInfo.txt'
	phrase_obs = 'upperLimit: '
	phrase_exp = 'expectedUpperLimit: '
	pos_obs = 5
	pos_exp	= 6

	for key in UL_obs:

		sr = key

		path = sr + '/' + name
		
		with open(path) as file:
			out = file.read()
			lines = out.split('\n')
			old_obs = lines[pos_obs].split(phrase_obs)[1]
			old_exp = lines[pos_exp].split(phrase_exp)[1]

			#ratio_obs = round(UL_obs[sr] / float(old_obs[:-3]), 3)
			#ratio_exp = round(UL_exp[sr] / float(old_exp[:-3]), 3)
			#print(ratio_obs, ratio_exp)

			lines[pos_obs] = '{}{}*fb\t# official number taken from ATLAS (calculated value: {})'.format(phrase_obs, UL_obs[sr], old_obs)
			lines[pos_exp] = '{}{}*fb\t# official number taken from ATLAS (calculated value: {})'.format(phrase_exp, UL_exp[sr], old_exp)
			
			with open(name, 'w') as temp:
				for line in lines:
					temp.write(line + '\n')

		remove(path)
		move(name, path)


UL_obs = {'SR2-SF-loose' : 2.02,
		  'SR2-SF-tight' : 0.29,
		  'SR2-DF-100'   : 0.88,
		  'SR2-DF-150'	 : 0.32,
	 	  'SR2-DF-200'	 : 0.33,
		  'SR2-DF-300'	 : 0.18,
		  'SR2-low'		 : 0.43,
		  'SR2-int'		 : 0.13, 
		  'SR2-high'	 : 0.09,
		  'WZ-0Ja'		 : 0.35, 
		  'WZ-0Jb'		 : 0.10, 
		  'WZ-0Jc'		 : 0.13, 
		  'WZ-1Ja'		 : 0.09, 
		  'WZ-1Jb'		 : 0.16, 
		  'WZ-1Jc'		 : 0.20,
		  'slep-a'		 : 0.19,
		  'slep-b'		 : 0.14,
		  'slep-c'		 : 0.29, 
		  'slep-d'		 : 0.08, 
		  'slep-e'		 : 0.09}


# S_obs/S_exp

UL_exp = {'SR2-SF-loose' : 73./53.,
		  'SR2-SF-tight' : 10.5/12.,
		  'SR2-DF-100'   : 32./27.,
		  'SR2-DF-150'	 : 11.4/12.,
	 	  'SR2-DF-200'	 : 12.0/10.3,
		  'SR2-DF-300'	 : 6.6/5.6,
		  'SR2-low'		 : 15.7/12.,
		  'SR2-int'		 : 4.5/5.6, 
		  'SR2-high'	 : 3.1/3.1,
		  'WZ-0Ja'		 : 12.8/14., 
		  'WZ-0Jb'		 : 3.7/4.6, 
		  'WZ-0Jc'		 : 4.8/4.1, 
		  'WZ-1Ja'		 : 3.2/4.5, 
		  'WZ-1Jb'		 : 5.6/4.3, 
		  'WZ-1Jc'		 : 7.2/4.2,
		  'slep-a'		 : 6.8/4.7,
		  'slep-b'		 : 5.2/5.1,
		  'slep-c'		 : 10.5/6.8, 
		  'slep-d'		 : 3.0/3.6, 
		  'slep-e'		 : 3.3/3.6}


# <sigma*eps> / (S_obs/S_exp)

for sr, ul in UL_exp.items(): UL_exp[sr] = round(UL_obs[sr] / ul, 3)

replaceULs(UL_obs, UL_exp)





