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
info.prettyName  = '0L stop + MET'
info.private 	 = False
info.arxiv 		 = 'https://arxiv.org/abs/1709.04183'
info.contact 	 = 'atlas-phys-susy-conveners@cern.ch'
info.publication = 'https://link.springer.com/article/10.1007/JHEP12(2017)085'
info.comment	 = 'Calculated obs and exp ULs replaced with official ATLAS data from https://atlas.web.cern.ch/Atlas/GROUPS/PHYSICS/PAPERS/SUSY-2016-15/tab_14.png'


obsN 	= [27, 11]
expN 	= [25.1, 8.5]
bgErr 	= [6.2, 1.5]
SR	 	= ['SRDlow', 'SRDhigh']
fig		= ['11a',	'11a']

for i in range(len(obsN)):
	#+++++++ dataset block ++++++++++++++
	dataset = DataSetInput(SR[i])
	dataset.setInfo(dataType = 'efficiencyMap', dataId = SR[i], observedN = obsN[i], expectedBG = expN[i], bgError = bgErr[i])
	#+++++++ next txName block ++++++++++++++
	T2bb 									= dataset.addTxName('T2bb')
	T2bb.checked 							= 'N/A'
	T2bb.constraint 						= '[[[b]],[[b]]]'
	T2bb.conditionDescription 				= None
	T2bb.condition 							= None
	T2bb.source 							= 'ATLAS'
	#+++++++ next mass plane block ++++++++++++++
	T2bb_1 									= T2bb.addMassPlane(2*[[x, y]])
	T2bb_1.figure 							= 'Fig.' + fig[i]
	T2bb_1.figureUrl 						= 'https://atlas.web.cern.ch/Atlas/GROUPS/PHYSICS/PAPERS/SUSY-2016-15/figaux_' + fig[i] + '.png'
	T2bb_1.dataUrl 							= 'https://www.hepdata.net/record/ins1623207?version=7&table=Eff' + SR[i]
	T2bb_1.setSources(dataLabels 			= ['efficiencyMap'],
								dataFiles 	= ['orig/EffMap_T2bb_' + SR[i] + '.txt'],
								coordinates = [{x : 1, y: 0, 'value' :2}],																	 
								dataFormats	= ['txt'])


obsN 	= [18, 	 	   9,	   11,	   206,		53,	 	 38,	  20,	  22,	  22,	  1,	  0,	  3]
expN 	= [18.7, 	 9.3,	  8.6, 	   179.,	52.4,	 39.3,	  20.6,	  27.6,	  18.9,	  7.7,	  0.91,	  3.64]
bgErr 	= [2.7, 	 2.2,	  2.1, 	   26.,		7.4,	 7.6,	  6.5,	  4.9,	  3.4,	  1.2,	  0.73,	  0.79]
SR 		= ['SRAT0', 'SRATW', 'SRATT', 'SRBT0', 'SRBTW', 'SRBTT', 'SRC1', 'SRC2', 'SRC3', 'SRC4', 'SRC5', 'SRE']
fig		= ['09c',	'09b',	 '09a',	  '09c',   '09b',	'09a',	 '10a',	 '10b',	 '10c',	 '10d',	 '10e',	 '11b']

for i in range(len(obsN)):
	#+++++++ dataset block ++++++++++++++
	dataset = DataSetInput(SR[i])
	dataset.setInfo(dataType = 'efficiencyMap', dataId = SR[i], observedN = obsN[i], expectedBG = expN[i], bgError = bgErr[i])
	#+++++++ next txName block ++++++++++++++
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
	'''
	#+++++++ next mass plane block ++++++++++++++
	T2bbffff						= dataset.addTxName('T2bbffff')
	T2bbffff.checked				= 'False'
	T2bbffff.constraint				= "[[['b', 'jet','jet']],[['b', 'jet','jet']]]"
	T2bbffff.conditionDescription 	= None
	T2bbffff.condition				= None
	T2bbffff.source					= 'ATLAS'
	T2bbffff.massConstraint			= [['dm < 80'], ['dm < 80']]
	'''
	#+++++++ next mass plane block ++++++++++++++
	T2tt_1 							= T2tt.addMassPlane(2*[[x, y]])
	T2tt_1.figure 					= 'Fig.' + fig[i]
	T2tt_1.figureUrl 				= 'https://atlas.web.cern.ch/Atlas/GROUPS/PHYSICS/PAPERS/SUSY-2016-15/figaux_' + fig[i] + '.png'
	T2tt_1.dataUrl 					= 'https://www.hepdata.net/record/ins1623207?version=7&table=Eff' + SR[i]
	T2tt_1.setSources(dataLabels 	= ['expExclusion', 'obsExclusion', 'efficiencyMap'],
						dataFiles 	= ['orig/ExpectedexclusioncontourdirectTT.csv', 'orig/ObservedexclusioncontourdirectTT.csv', 'orig/EffMap_T2tt_' + SR[i] + '.txt'],
						coordinates = [ {x: 0, y: 1, 'value': None}, {x: 0, y: 1, 'value': None},  {x : 1, y: 0, 'value' :2} ],				 
						dataFormats	= ['csv', 'csv', 'txt'])
	'''
	T2bbffff_1 						 = T2bbffff.addMassPlane(2*[[x, y]])
	T2bbffff_1.figure 				 = 'Fig.' + fig[i]
	T2bbffff_1.figureUrl 			 = 'https://atlas.web.cern.ch/Atlas/GROUPS/PHYSICS/PAPERS/SUSY-2016-15/figaux_' + fig[i] + '.png'
	T2bbffff_1.dataUrl 				 = 'https://www.hepdata.net/record/ins1623207?version=7&table=Eff' + SR[i]
	T2bbffff_1.setSources(dataLabels = ['expExclusion', 'obsExclusion', 'efficiencyMap'],
						dataFiles 	 = ['orig/ExpectedexclusioncontourdirectTT.csv', 'orig/ObservedexclusioncontourdirectTT.csv', 'orig/EffMap_T2tt_' + SR[i] + '.txt'],
						coordinates  = [ {x: 0, y: 1, 'value': None}, {x: 0, y: 1, 'value': None},  {x : 1, y: 0, 'value' :2} ],
						scales		 = [1., 1., 9./4.],				 
						dataFormats	 = ['csv', 'csv', 'txt'])
	'''
	T2ttoff.addMassPlane(T2tt_1)
	#T2bbffff.addMassPlane(T2tt_1)





databaseCreator.create()


from shutil import move
from os import remove
### manually replace ULs in dataInfo.txt
### taken from https://atlas.web.cern.ch/Atlas/GROUPS/PHYSICS/PAPERS/SUSY-2016-15/tab_14.png

def replaceULs(UL_obs, UL_exp):

	for sr, ul in UL_exp.items(): UL_exp[sr] = round(ul / 36.1, 3)

	name = 'dataInfo.txt'
	phrase_obs = 'upperLimit: '
	phrase_exp = 'expectedUpperLimit: '
	pos_obs = 5
	pos_exp	= 6

	for sr in SR:

		path = sr + '/' + name
		
		with open(path) as file:
			out = file.read()
			lines = out.split('\n')
			old_obs = lines[pos_obs].split(phrase_obs)[1]
			old_exp = lines[pos_exp].split(phrase_exp)[1]
			lines[pos_obs] = '{}{}*fb\t# official number taken from ATLAS (table 14) (calculated value: {})'.format(phrase_obs, UL_obs[sr], old_obs)
			lines[pos_exp] = '{}{}*fb\t# official number taken from ATLAS (table 14) (calculated value: {})'.format(phrase_exp, UL_exp[sr], old_exp)
			
			with open(name, 'w') as temp:
				for line in lines:
					temp.write(line + '\n')

		remove(path)
		move(name, path)




UL_obs = {'SRAT0' : 0.31, 
		  'SRATW' : 0.27, 
		  'SRATT' : 0.30, 
		  'SRBT0' : 2.19,
		  'SRBTW' : 0.60,
		  'SRBTT' : 0.54,
		  'SRC1'  : 0.42,
		  'SRC2'  : 0.31,
		  'SRC3'  : 0.42,
		  'SRC4'  : 0.10,
		  'SRC5'  : 0.09,
		  'SRE'   : 0.17}

UL_exp = {'SRAT0' : 11.5, 
		  'SRATW' : 9.6, 
		  'SRATT' : 8.7, 
		  'SRBT0' : 58.0,
		  'SRBTW' : 21.0,
		  'SRBTT' : 20.0,
		  'SRC1'  : 15.8,
		  'SRC2'  : 13.9,
		  'SRC3'  : 12.3,
		  'SRC4'  : 6.7,
		  'SRC5'  : 3.0,
		  'SRE'   : 6.4}

replaceULs(UL_obs, UL_exp)

SR = ['SRDlow', 'SRDhigh']
UL_obs = {'SRDlow' : 0.5, 
		  'SRDhigh' : 0.3}

UL_exp = {'SRDlow' : 16.4, 
		  'SRDhigh' : 8.0}

replaceULs(UL_obs, UL_exp)

