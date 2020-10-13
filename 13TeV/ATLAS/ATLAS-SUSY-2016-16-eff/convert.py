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
info.lumi 		 = 36
info.prettyName  = '1L stop'
info.private 	 = False
info.arxiv 		 = 'https://arxiv.org/abs/1711.11520'
info.contact 	 = 'atlas-phys-susy-conveners@cern.ch'
info.publication = 'JHEP 06 (2018) 108'
info.comment 	 = 'Information about SR-mixing in official data missing for T6bbWW' 




SR   = {'obsN' 	: [8,50,19,115,34,68,70], 
		'expN'  : [3.8,36.3,18.3,115,30.3,71,60.5], 
		'bgErr' : [1,6.6,2.2,31,5.9,16,6.1], 
		'SR' 	: ['tN_high','tN_med','tN_diag_high','tN_diag_med','tN_diag_low','bWN','bffN']}

SR2  = {'obsN'  : [5,13,65,22,4,25,33,19,2],
        'expN'  : [7.4,13.8,48.3,21.3,5.8,25.1,24.7,13.7,1.8],
        'bgErr' : [2.1,3.6,8.2,5.0,1.6,3.8,3.1,2.1,0.3],
        'SR'    : ['DM_high','DM_low','DM_low_loose','bC2x_diag','bC2x_med','bCbv','bCsoft_diag','bCsoft_med','bCsoft_high']}

SR3  = {'obsN'  : [83,52,69,33],
        'expN'  : [61,38.2,50,24.7],
        'bgErr' : [9.7,6.9,8.7,3.1],
        'SR'    : ['bCsoft_diag+tN_med','bCsoft_high+tN_med','bCsoft_med+tN_med', 'bCsoft_diag']}

T2tt = {
'name' 		 : 'T2tt',
'sources'	 :{'expExcl'		: 'orig/HEPData-ins1639856-v4-Table_16.csv',
			   'obsExcl'		: 'orig/HEPData-ins1639856-v4-Table_17.csv',
			   'effMap'			: 'orig/EffMap_T2tt_'},
'constraint' : '[[[t]],[[t]]]',
'massConstr' : None,
'massPlane'  : 2*[[x, y]]}

T2ttoff = {
'name' 		 : 'T2ttoff',
'sources'	 :{'expExcl'		: 'orig/HEPData-ins1639856-v4-Table_16.csv',
			   'obsExcl'		: 'orig/HEPData-ins1639856-v4-Table_17.csv',
			   'effMap'			: 'orig/EffMap_T2tt_'},
'constraint' : '[[[b, W]],[[b, W]]]',
'massConstr' : [['80 <= dm < 169.0'], ['80 <= dm < 169.0']],
'massPlane'  : 2*[[x, y]]}

T2bbffff = {
'name' 		 : 'T2bbffff',
'sources'	 :{'expExcl'		: 'orig/HEPData-ins1639856-v4-Table_16.csv',
			   'obsExcl'		: 'orig/HEPData-ins1639856-v4-Table_17.csv',
			   'effMap'			: 'orig/EffMap_T2tt_'},
'constraint' : '[[[b, l, nu]],[[b, jet,jet]]]',
'massConstr' : [['dm < 80'], ['dm < 80']],
'massPlane'  : 2*[[x, y]]}

for i in range(len(SR['SR'])):
	if i < 5: TX = T2tt
	elif i == 5: TX = T2ttoff
	else: TX = T2bbffff
	#+++++++ dataset block ++++++++++++++
	dataset = DataSetInput(SR['SR'][i])
	dataset.setInfo(dataType = 'efficiencyMap', dataId = SR['SR'][i], observedN = SR['obsN'][i], expectedBG = SR['expN'][i], bgError = SR['bgErr'][i])
	#+++++++ next txName block ++++++++++++++
	newTx							= dataset.addTxName(TX['name'])
	newTx.checked					= 'False'
	newTx.constraint				= TX['constraint']
	newTx.conditionDescription 		= None
	newTx.condition					= None
	newTx.source					= 'ATLAS'
	newTx.massConstraint			= TX['massConstr']
	#+++++++ next mass plane block ++++++++++++++
	newPlane 						= newTx.addMassPlane(TX['massPlane'])
	newPlane.figure 				= 'Fig14a;Fig14b'
	newPlane.figureUrl 				= 'https://atlas.web.cern.ch/Atlas/GROUPS/PHYSICS/PAPERS/SUSY-2016-16/figaux_14a.png;https://atlas.web.cern.ch/Atlas/GROUPS/PHYSICS/PAPERS/SUSY-2016-16/figaux_14b.png'
	newPlane.dataUrl 				= 'https://www.hepdata.net/record/ins1639856?version=4&table=Table87'
	if TX != T2bbffff:
		newPlane.setSources(dataLabels 	= ['expExclusion', 'obsExclusion', 'efficiencyMap'],
						dataFiles 		= [TX['sources']['expExcl'], TX['sources']['obsExcl'], TX['sources']['effMap'] + SR['SR'][i] + '.txt'],
						#units			= [ None, None, 'pb' ],
					 	coordinates 	= [ {x: 0, y: 1, 'value': None}, {x: 0, y: 1, 'value': None},  {x : 1, y: 0, 'value' :2} ],
		             	dataFormats 	= ['csv', 'csv', 'txt'])
	else:
		newPlane.setSources(dataLabels 	= ['expExclusion', 'obsExclusion', 'efficiencyMap'],
						dataFiles 		= [TX['sources']['expExcl'], TX['sources']['obsExcl'], TX['sources']['effMap'] + SR['SR'][i] + '.txt'],
						#units			= [ None, None, 'pb' ],
					 	coordinates 	= [ {x: 0, y: 1, 'value': None}, {x: 0, y: 1, 'value': None},  {x : 1, y: 0, 'value' :2} ],
		             	dataFormats 	= ['csv', 'csv', 'txt'],
						scales = [1, 1, 27./8.]) #(9./4.)*(3./2.)

databaseCreator.create()


### manually replace ULs in dataInfo.txt
### taken from https://atlas.web.cern.ch/Atlas/GROUPS/PHYSICS/PAPERS/SUSY-2016-16/tab_23.png

from shutil import move
from os import remove

UL_obs = {'tN_high' 		: 10, 
		  'tN_med' 			: 31, 
		  'tN_diag_high' 	: 11, 
		  'tN_diag_med' 	: 58,
		  'tN_diag_low' 	: 17,
		  'bWN'				: 31,
		  'bffN'			: 28}

UL_exp = {'tN_high' 		: 5.8, 
		  'tN_med' 			: 19, 
		  'tN_diag_high' 	: 11, 
		  'tN_diag_med' 	: 58,
		  'tN_diag_low' 	: 19,
		  'bWN'				: 33,
		  'bffN'			: 21}

for sr, ul in UL_obs.items(): UL_obs[sr] = round(ul / 36., 3)
for sr, ul in UL_exp.items(): UL_exp[sr] = round(ul / 36., 3)

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

