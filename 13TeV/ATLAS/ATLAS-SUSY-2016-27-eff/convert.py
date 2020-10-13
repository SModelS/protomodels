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
info 				= MetaInfoInput('ATLAS-SUSY-2016-27')
info.url 			= "https://atlas.web.cern.ch/Atlas/GROUPS/PHYSICS/PAPERS/SUSY-2016-27/"
info.sqrts 			= 13
info.lumi			= 36.1
info.prettyName 	= "jets + photon + Etmiss"
info.private		= False
info.arxiv			= 'https://arxiv.org/abs/1802.03158'
info.contact		= 'atlas-phys-susy-conveners@cern.ch'
info.publication	= 'https://journals.aps.org/prd/abstract/10.1103/PhysRevD.97.092006'

obsN 	= [0, 0]
expN 	= [0.5, 0.48]
bgErr 	= [0.28, 0.275]
SR 		= ['SRyy-SL', 'SRyy-SH']
TP		= ['T5gg', 'T6gg']
DATA    = [['https://www.hepdata.net/record/ins1654357?version=1&table=Acceptance/Efficiency 1', 'https://www.hepdata.net/record/ins1654357?version=1&table=Acceptance/Efficiency 2'],
		   ['https://www.hepdata.net/record/ins1654357?version=1&table=Acceptance/Efficiency 3', 'https://www.hepdata.net/record/ins1654357?version=1&table=Acceptance/Efficiency 4']]


for i in range(len(obsN)):

	#+++++++ dataset block ++++++++++++++
	dataset = DataSetInput(SR[i])
	dataset.setInfo(dataType = 'efficiencyMap', dataId = SR[i] + '-q', observedN = obsN[i], expectedBG = expN[i], bgError = bgErr[i])
	#+++++++ next txName block ++++++++++++++
	T6gg 					 	= dataset.addTxName(TP[1])
	T6gg.checked 			 	= 'No'
	T6gg.constraint 			 	= "[[['jet'],['photon']],[['jet'],['photon']]]"
	T6gg.conditionDescription 	= None
	T6gg.condition 			 	= None
	T6gg.source 				 	= "ATLAS"
	#+++++++ next mass plane block ++++++++++++++
	T6gg_1 						= T6gg.addMassPlane(2*[[x, y, 1.]])
	T6gg_1.dataUrl   			= DATA[1][i]
	T6gg_1.setSources(dataLabels	= ['expExclusion', 'obsExclusion', 'efficiencyMap'],
						dataFiles 	= ['orig/HEPData-ins1654357-v1-Exclusion_contour_(expected)_3.csv', 'orig/HEPData-ins1654357-v1-Exclusion_contour_(observed)_4.csv', 'orig/EffMap_' + TP[1].replace("gg","Gamma") + '_' + SR[i] + '.txt'],																	 
						dataFormats	= ['csv', 'csv', 'txt'])

	T5gg 						= dataset.addTxName(TP[0])
	T5gg.checked 				= 'No'
	T5gg.constraint 		 		= "[[['jet','jet'],['photon']],[['jet','jet'],['photon']]]"
	T5gg.conditionDescription	= None
	T5gg.condition 		 		= None
	T5gg.source 					= "ATLAS"
	#+++++++ next mass plane block ++++++++++++++
	T5gg_1						= T5gg.addMassPlane(2*[[x, y, 1.]])
	T5gg_1.dataUrl   			= DATA[0][i]
	T5gg_1.setSources(dataLabels	= ['expExclusion', 'obsExclusion', 'efficiencyMap'],
						dataFiles 	= ['orig/HEPData-ins1654357-v1-Exclusion_contour_(expected)_1.csv', 'orig/HEPData-ins1654357-v1-Exclusion_contour_(observed)_2.csv', 'orig/EffMap_' + TP[0].replace("gg","Gamma") + '_' + SR[i] + '.txt'],																	 
						dataFormats	= ['csv', 'csv', 'txt'])


obsN 	= [6, 1]
expN 	= [3.7, 2.05]
bgErr 	= [1.1, 0.64]
SR 		= ['SRyy-WL', 'SRyy-WH']
TP		= 'TChipChimgg'
DATA    = ['https://www.hepdata.net/record/ins1654357?version=1&table=Acceptance/Efficiency 5', 'https://www.hepdata.net/record/ins1654357?version=1&table=Acceptance/Efficiency 6']

for i in range(len(obsN)):

	#+++++++ dataset block ++++++++++++++
	dataset = DataSetInput(SR[i])
	dataset.setInfo(dataType = 'efficiencyMap', dataId = SR[i], observedN = obsN[i], expectedBG = expN[i], bgError = bgErr[i])
	#+++++++ next txName block ++++++++++++++
	TChipChimgg 							= dataset.addTxName(TP)
	TChipChimgg.checked 					= 'No'
	TChipChimgg.constraint 				= "[[['W'],['photon']],[['Z'],['photon']]]+[[['W'],['photon']],[['W'],['photon']]]+[[['W'],['photon']],[['higgs'],['photon']]]"
	TChipChimgg.conditionDescription 	= None
	TChipChimgg.condition 				= None
	TChipChimgg.source 			    	= "ATLAS"
	#+++++++ next mass plane block ++++++++++++++
	TChipChimgg_1 		   				= TChipChimgg.addMassPlane(2*[[x, y, 1.]])
	TChipChimgg_1.dataUrl   				= DATA[i]
	TChipChimgg_1.setSources(dataLabels	= ['expExclusion', 'obsExclusion', 'efficiencyMap'],
								dataFiles 	= ['orig/HEPData-ins1654357-v1-Exclusion_contour_(expected)_5.csv', 'orig/HEPData-ins1654357-v1-Exclusion_contour_(observed)_6.csv', 'orig/EffMap_' + TP.replace("gg","Gamma") + '_' + SR[i] + '.txt'],																	 
								dataFormats	= ['csv', 'csv', 'txt'])


obsN 	= [4, 3]
expN 	= [1.33, 1.14]
bgErr 	= [0.38, 0.485]
SR 		= ['SRyj-L', 'SRyj-H']
TP		= 'T5Zg'
DATA    = ['https://www.hepdata.net/record/ins1654357?version=1&table=Acceptance/Efficiency 7', 'https://www.hepdata.net/record/ins1654357?version=1&table=Acceptance/Efficiency 8']

'''

obsN 	= [4, 3, 0, 0]
expN 	= [1.33, 1.14, 0.5, 0.48]
bgErr 	= [0.38, 0.485, 0.28, 0.275]
SR 		= ['SRyj-L', 'SRyj-H', 'SRyy-SL', 'SRyy-SH']
TP		= 'T5Zg'
DATA    = ['https://www.hepdata.net/record/ins1654357?version=1&table=Acceptance/Efficiency 7', 'https://www.hepdata.net/record/ins1654357?version=1&table=Acceptance/Efficiency 8', 'https://www.hepdata.net/record/ins1654357?version=1&table=Acceptance/Efficiency 1', 'https://www.hepdata.net/record/ins1654357?version=1&table=Acceptance/Efficiency 2']

'''

for i in range(len(obsN)):

	#+++++++ dataset block ++++++++++++++
	dataset = DataSetInput(SR[i])
	dataset.setInfo(dataType = 'efficiencyMap', dataId = SR[i] + '-q', observedN = obsN[i], expectedBG = expN[i], bgError = bgErr[i])
	#+++++++ next txName block ++++++++++++++
	T5Zg 					  		= dataset.addTxName(TP)
	T5Zg.checked 			  		= 'No'
	T5Zg.constraint 		  		= "[[['jet','jet'],['Z']],[['jet','jet'],['photon']]]+[[['jet','jet'],['photon']],[['jet','jet'],['photon']]]+[[['jet','jet'],['Z']],[['jet','jet'],['Z']]]"#"[[['jet','jet'],['Z']],[['jet','jet'],['photon']]]"
	T5Zg.conditionDescription 		= None
	T5Zg.condition 			  		= "Csim([[['jet','jet'],['Z']],[['jet','jet'],['photon']]]+2.*[[['jet','jet'],['photon']],[['jet','jet'],['photon']]]+2.*[[['jet','jet'],['Z']],[['jet','jet'],['Z']]])"
	T5Zg.source 			  		= "ATLAS"
	#+++++++ next mass plane block ++++++++++++++
	T5Zg_1 							= T5Zg.addMassPlane(2*[[x, y, 1.]])
	T5Zg_1.dataUrl   				= DATA[i]
	T5Zg_1.setSources(dataLabels	= ['expExclusion', 'obsExclusion', 'efficiencyMap'],
						dataFiles 	= ['orig/HEPData-ins1654357-v1-Exclusion_contour_(expected)_7.csv', 'orig/HEPData-ins1654357-v1-Exclusion_contour_(observed)_8.csv', 'orig/EffMap_' + TP.replace("g","Gamma") + '_' + SR[i] + '.txt'],																	 
						dataFormats	= ['csv', 'csv', 'txt'])



databaseCreator.create()
