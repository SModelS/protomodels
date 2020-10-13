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
info 				= MetaInfoInput('ATLAS-SUSY-2017-02')
info.url 			= 'https://atlas.web.cern.ch/Atlas/GROUPS/PHYSICS/PAPERS/SUSY-2017-02/'
info.sqrts 			= 13
info.lumi 			= 36.1
info.prettyName 	= '0L + jets + Etmiss'
info.private 		= False
info.arxiv 			= 'https://arxiv.org/abs/1806.04030'
info.contact 		= 'atlas-phys-susy-conveners@cern.ch'
info.publication 	= 'https://journals.aps.org/prd/abstract/10.1103/PhysRevD.98.092002'
info.comment		= 'Low-mass results (figaux 2) have been added to the upper limit map of the high-mass analysis (figaux 3) for higgsino masses <= 300GeV.'

#+++++++ dataset block ++++++++++++++
dataset = DataSetInput('data')
dataset.setInfo(dataType = 'upperLimit', dataId = None)

lsp_masses = [1.]

#+++++++ next txName block ++++++++++++++
TChiH 						= dataset.addTxName('TChiH')
TChiH.checked 				= 'No'
TChiH.constraint 			= "[[['Z']],[['Z']]]+[[['higgs']],[['higgs']]]+[[['higgs']],[['Z']]]"
TChiH.conditionDescription = None
TChiH.condition 			= None
TChiH.source 				= "ATLAS"
#+++++++ next mass plane block ++++++++++++++
for lsp in lsp_masses:
	plane 						= TChiH.addMassPlane(2*[[x, y]])
	plane.figure 				= 'figaux_04'
	plane.figureUrl 			= 'https://atlas.web.cern.ch/Atlas/GROUPS/PHYSICS/PAPERS/SUSY-2017-02/figaux_04.png'
	plane.dataUrl 				= 'https://www.hepdata.net/record/ins1677389?version=1&table=Table3'
	plane.setSources(dataLabels = ['upperLimits', 'expectedUpperLimits'],
					units 		= ['fb', 'fb'],
					dataFiles 	= ['orig/HEPData-ins1677389-v1-Table_3_addedG_obs.csv', 'orig/HEPData-ins1677389-v1-Table_3_addedG_exp.csv'],
          			coordinates = [{ x: 0, y: 1, "value": 2 }, { x: 0, y: 1, "value": 2 }],
					dataFormats	= ['csv', 'csv'])

'''

for lsp in lsp_masses:
	plane 						= TChiH.addMassPlane(2*[[x, lsp]])
	plane.figure 				= 'figaux_04'
	plane.figureUrl 			= 'https://atlas.web.cern.ch/Atlas/GROUPS/PHYSICS/PAPERS/SUSY-2017-02/figaux_04.png'
	plane.dataUrl 				= 'https://www.hepdata.net/record/ins1677389?version=1&table=Table3'
	plane.setSources(dataLabels = ['upperLimits', 'expectedUpperLimits'],
					units 		= ['fb', 'fb'],
					dataFiles 	= ['orig/HEPData-ins1677389-v1-Table_3_addedG_obs.csv', 'orig/HEPData-ins1677389-v1-Table_3_addedG_exp.csv'],
          			coordinates = [{ x: 0, "value": 2 }, { x: 0, "value": 2 }],
					dataFormats	= ['csv', 'csv'])
'''

databaseCreator.create()
