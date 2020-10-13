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

#+++++++ dataset block ++++++++++++++
dataset = DataSetInput('data')
dataset.setInfo(dataType = 'upperLimit', dataId = None)

lsp_masses = [1.]
#planes = []

#+++++++ next txName block ++++++++++++++
T5gg 						= dataset.addTxName('T5gg')
T5gg.checked 				= 'No'
T5gg.constraint 		 		= "[[['jet','jet'],['photon']],[['jet','jet'],['photon']]]"
T5gg.conditionDescription	= None
T5gg.condition 		 		= None
T5gg.source 					= "ATLAS"
#+++++++ next mass plane block ++++++++++++++
for lsp in lsp_masses:
	p						= T5gg.addMassPlane(2*[[x, y, lsp]])
	p.figure    			= 'Fig.8'
	p.figureUrl 			= 'https://atlas.web.cern.ch/Atlas/GROUPS/PHYSICS/PAPERS/SUSY-2016-27/fig_08.png'
	p.dataUrl   			= "https://www.hepdata.net/record/ins1654357?version=1&table=Cross section UL 1"
	p.setSources(dataLabels	= ['expExclusion', 'obsExclusion', 'upperLimits'],
				dataFiles	= ['orig/HEPData-ins1654357-v1-Exclusion_contour_(expected)_1.csv', 'orig/HEPData-ins1654357-v1-Exclusion_contour_(observed)_2.csv', 'orig/HEPData-ins1654357-v1-Cross_section_UL_1.csv'],
				units 		= [ None, None, 'fb' ],
				dataFormats	= ['csv', 'csv', 'csv'])
#	planes.append(p)

#+++++++ next txName block ++++++++++++++
T6gg 					 	= dataset.addTxName('T6gg')
T6gg.checked 			 	= 'No'
T6gg.constraint 			 	= "[[['jet'],['photon']],[['jet'],['photon']]]"#+[[['b'],['photon']],[['b'],['photon']]]"
T6gg.conditionDescription 	= None
T6gg.condition 			 	= None
T6gg.source 				 	= "ATLAS"
#+++++++ next mass plane block ++++++++++++++
for lsp in lsp_masses:
	p 						= T6gg.addMassPlane(2*[[x, y, lsp]])
	p.figure    			= 'Fig.9'
	p.figureUrl 			= 'https://atlas.web.cern.ch/Atlas/GROUPS/PHYSICS/PAPERS/SUSY-2016-27/fig_09.png'
	p.dataUrl   			= "https://www.hepdata.net/record/ins1654357?version=1&table=Cross section UL 2"
	p.setSources(dataLabels = ['expExclusion', 'obsExclusion', 'upperLimits'],
				dataFiles 	= ['orig/HEPData-ins1654357-v1-Exclusion_contour_(expected)_3.csv', 'orig/HEPData-ins1654357-v1-Exclusion_contour_(observed)_4.csv', 'orig/HEPData-ins1654357-v1-Cross_section_UL_2.csv'],
					units 	= [ None, None, 'fb' ],
				dataFormats = ['csv', 'csv', 'csv'])
#	planes.append(p)

#+++++++ next txName block ++++++++++++++
TChipChimgg 							= dataset.addTxName('TChipChimgg')
TChipChimgg.checked 					= 'No'
TChipChimgg.constraint 				= "[[['W'],['photon']],[['Z'],['photon']]]+[[['W'],['photon']],[['W'],['photon']]]+[[['W'],['photon']],[['higgs'],['photon']]]"
TChipChimgg.conditionDescription 	= None
TChipChimgg.condition 				= None
TChipChimgg.source 			    	= "ATLAS"
#+++++++ next mass plane block ++++++++++++++
for lsp in lsp_masses:
	p 		   				= TChipChimgg.addMassPlane(2*[[x, y, lsp]])
	p.figure    			= 'Fig.10'
	p.figureUrl 			= 'https://atlas.web.cern.ch/Atlas/GROUPS/PHYSICS/PAPERS/SUSY-2016-27/fig_10.png'
	p.dataUrl   			= "https://www.hepdata.net/record/ins1654357?version=1&table=Cross section UL 3"
	p.setSources(dataLabels = ['expExclusion', 'obsExclusion', 'upperLimits'],
				dataFiles	= ['orig/HEPData-ins1654357-v1-Exclusion_contour_(expected)_5.csv', 'orig/HEPData-ins1654357-v1-Exclusion_contour_(observed)_6.csv', 'orig/HEPData-ins1654357-v1-Cross_section_UL_3.csv'],
					units	= [ None, None, 'fb' ],
				dataFormats	= ['csv', 'csv', 'csv'])
#	planes.append(p)

#+++++++ next txName block ++++++++++++++
T5Zg 					  		= dataset.addTxName('T5Zg')
T5Zg.checked 			  		= 'No'
T5Zg.constraint 		  		= "[[['jet','jet'],['Z']],[['jet','jet'],['photon']]]+[[['jet','jet'],['photon']],[['jet','jet'],['photon']]]+[[['jet','jet'],['Z']],[['jet','jet'],['Z']]]"
T5Zg.conditionDescription 		= None
T5Zg.condition 			  		= "Csim([[['jet','jet'],['Z']],[['jet','jet'],['photon']]]+2.*[[['jet','jet'],['photon']],[['jet','jet'],['photon']]]+2.*[[['jet','jet'],['Z']],[['jet','jet'],['Z']]])"
T5Zg.source 			  		= "ATLAS"
#+++++++ next mass plane block ++++++++++++++
for lsp in lsp_masses:
	p 						= T5Zg.addMassPlane(2*[[x, y, lsp]])
	p.figure    			= 'Fig.11'
	p.figureUrl 			= 'https://atlas.web.cern.ch/Atlas/GROUPS/PHYSICS/PAPERS/SUSY-2016-27/fig_11.png'
	p.dataUrl   			= "https://www.hepdata.net/record/ins1654357?version=1&table=Cross section UL 4"
	p.setSources(dataLabels	= ['expExclusion', 'obsExclusion', 'upperLimits'],
				dataFiles	= ['orig/HEPData-ins1654357-v1-Exclusion_contour_(expected)_7.csv', 'orig/HEPData-ins1654357-v1-Exclusion_contour_(observed)_8.csv', 'orig/HEPData-ins1654357-v1-Cross_section_UL_4.csv'],
					units	= [ None, None, 'fb' ],
				dataFormats	= ['csv', 'csv', 'csv'])
#	planes.append(p)


databaseCreator.create()
