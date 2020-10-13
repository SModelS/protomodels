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
info.comment		= 'Only 0lbb part of analysis used (for mass > 350GeV), 1lbb will follow with database update. Full-luminosity result available in ATLAS-SUSY-2019-08.'

#+++++++ dataset block ++++++++++++++
dataset = DataSetInput('data')
dataset.setInfo(dataType = 'upperLimit', dataId = None)
#+++++++ next txName block ++++++++++++++
TChiWH 							= dataset.addTxName('TChiWH')
TChiWH.checked 					= 'NO'
TChiWH.constraint 				= "[[['W']],[['higgs']]]"
TChiWH.conditionDescription 	= None
TChiWH.condition 				= None
TChiWH.source 					= "ATLAS"
#+++++++ next mass plane block ++++++++++++++
TChiWH_1 						= TChiWH.addMassPlane(2*[[x, y]])
TChiWH_1.figure 				= 'figaux_13'
TChiWH_1.figureUrl 				= 'https://atlas.web.cern.ch/Atlas/GROUPS/PHYSICS/PAPERS/SUSY-2017-01/figaux_13.png'
TChiWH_1.dataUrl 				= 'https://www.hepdata.net/record/ins1711261?version=1&table=Upper limit 0lbb'
TChiWH_1.setSources(dataLabels 	= ['expExclusion', 'obsExclusion', 'upperLimits'],
					units 		= [None, None, 'pb'],
					dataFiles 	= ['orig/Expectedlimit0lbb.csv', 'orig/Observedlimit0lbb.csv', 'orig/Upperlimit0lbb.csv'],
					dataFormats	= ['csv', 'csv', 'csv'])

databaseCreator.create()
