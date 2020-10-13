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

#+++++++ dataset block ++++++++++++++
dataset = DataSetInput('data')
dataset.setInfo(dataType = 'upperLimit', dataId = None)
#+++++++ next txName block ++++++++++++++
TChiWZoffqq 					= dataset.addTxName('TChiWZoffqq')
TChiWZoffqq.checked 				= 'no'
TChiWZoffqq.constraint 				= "23*([[['q','q']],[['mu+','mu-']]]+[[['q','q']],[['e+','e-']]])"
TChiWZoffqq.conditionDescription 			= None
TChiWZoffqq.condition 				=  "Cgtr([[['q','q']],[['mu+','mu-']]],[[['q','q']],[['e+','e-']]])"
TChiWZoffqq.source 				= "ATLAS"
TChiWZoffqq.massConstraint			=[['dm < 76'],[' dm < 87']]

#+++++++ next mass plane block ++++++++++++++
TChiWZoffqq_1 					= TChiWZoffqq.addMassPlane(2*[[x, (x-y)]])
TChiWZoffqq_1.figure 				= 'Fig.41ab'
TChiWZoffqq_1.figureUrl 				= 'https://atlas.web.cern.ch/Atlas/GROUPS/PHYSICS/PAPERS/SUSY-2018-16/figaux_41a.png'
TChiWZoffqq_1.dataUrl 				= 'https://www.hepdata.net/record/ins1767649?version=1&table=Figure 41ab'
TChiWZoffqq_1.exclusionDataUrl			= 'https://www.hepdata.net/record/ins1767649?version=1&table=Figure 14c Observed'
TChiWZoffqq_1.setSources(dataLabels 	= ['expExclusion', 'obsExclusion','upperLimits','expectedUpperLimits' ],
					units 		= [None, None,'pb','pb'],
					dataFiles 	= ['orig/Figure14cExpected.csv','orig/Figure14cObserved.csv','orig/Figure41ab-obs.csv', 'orig/Figure41ab-exp.csv'],
					dataFormats	= ['csv','csv','csv','csv'])

#+++++++ next txName block ++++++++++++++
TSlepSlep 					= dataset.addTxName('TSlepSlep')
TSlepSlep.checked 				= 'no'
TSlepSlep.constraint 				= "[[['mu+']],[['mu-']]]+[[['e+']],[['e-']]]"
TSlepSlep.conditionDescription 			= None
TSlepSlep.condition 				= None
TSlepSlep.source 				= "ATLAS"

#+++++++ next mass plane block ++++++++++++++
TSlepSlep_1 					= TSlepSlep.addMassPlane(2*[[x, (x-y)]])
TSlepSlep_1.figure 				= 'Fig.44ab'
TSlepSlep_1.figureUrl 				= 'https://atlas.web.cern.ch/Atlas/GROUPS/PHYSICS/PAPERS/SUSY-2018-16/figaux_44a.png'
TSlepSlep_1.dataUrl 				= 'https://www.hepdata.net/record/ins1767649?version=1&table=Figure 44ab'
TSlepSlep_1.exclusionDataUrl			= 'https://www.hepdata.net/record/ins1767649?version=1&table=Figure 16a Observed'
TSlepSlep_1.setSources(dataLabels 	= ['expExclusion', 'obsExclusion','upperLimits','expectedUpperLimits' ],
					units 		= [None, None,'pb','pb'],
					dataFiles 	= ['orig/Figure16aExpected.csv','orig/Figure16aObserved.csv','orig/Figure44ab-obs.csv', 'orig/Figure44ab-exp.csv'],
					dataFormats	= ['csv','csv','csv','csv'])



databaseCreator.create()
