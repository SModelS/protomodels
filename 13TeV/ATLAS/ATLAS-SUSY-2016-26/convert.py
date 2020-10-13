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
info = MetaInfoInput('ATLAS-SUSY-2016-26')
info.url = "https://atlas.web.cern.ch/Atlas/GROUPS/PHYSICS/PAPERS/SUSY-2016-26/"
info.sqrts = 13
info.lumi = 36.1
info.prettyName = ">=2 c jets + Etmiss"
info.private = False
info.arxiv =  'https://arxiv.org/abs/1805.01649'
info.contact = 'ATLAS collaboration'
info.publication ='JHEP 09 (2018) 050'

#+++++++ dataset block ++++++++++++++
dataset = DataSetInput('data')
dataset.setInfo(dataType = 'upperLimit', dataId = None)

#+++++++ next txName block ++++++++++++++
T2cc = dataset.addTxName('T2cc')
T2cc.checked = ''
T2cc.constraint ="[[['c']],[['c']]]"
T2cc.conditionDescription = None
T2cc.condition = None
T2cc.source = "ATLAS"
#+++++++ next mass plane block ++++++++++++++
T2cc_1 = T2cc.addMassPlane(2*[[x, y]])
T2cc_1.figure = 'Fig.6'
T2cc_1.figureUrl ='https://atlas.web.cern.ch/Atlas/GROUPS/PHYSICS/PAPERS/SUSY-2016-26/fig_06.png'
T2cc_1.dataUrl = "https://www.hepdata.net/record/ins1672099?version=3&table=Table%2030"
T2cc_1.setSources(dataLabels= ['expExclusion', 'obsExclusion', 'upperLimits'],
                 dataFiles= ['orig/T2cc_Exp_Excl.csv', 'orig/T2cc_Obs_Excl.csv', 'orig/T2cc_Obs_UL.csv'],
                 dataFormats= ['csv', 'csv', 'csv'])



databaseCreator.create()
