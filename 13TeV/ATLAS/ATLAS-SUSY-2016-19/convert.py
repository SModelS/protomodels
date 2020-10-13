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
info = MetaInfoInput('ATLAS-SUSY-2016-19')
info.url = "https://atlas.web.cern.ch/Atlas/GROUPS/PHYSICS/PAPERS/SUSY-2016-19/"
info.sqrts = 13
info.lumi = 36.1
info.prettyName = "stops to staus"
info.private = False
info.arxiv =  'https://arxiv.org/abs/1803.10178'
info.contact = 'ATLAS collaboration'
info.publication ='in process of being published'

#+++++++ dataset block ++++++++++++++
dataset = DataSetInput('data')
dataset.setInfo(dataType = 'upperLimit', dataId = None)

#+++++++ next txName block ++++++++++++++
T4bnutaubnutau = dataset.addTxName('T4bnutaubnutau')
T4bnutaubnutau.checked = ''
T4bnutaubnutau.constraint ="[[[b,nu],[ta]],[[b,nu],[ta]]]"
T4bnutaubnutau.conditionDescription = None
T4bnutaubnutau.condition = None
T4bnutaubnutau.source = "ATLAS"
#+++++++ next mass plane block ++++++++++++++
T4bnutaubnutau_1 = T4bnutaubnutau.addMassPlane(2*[[x, y, 0.]])
T4bnutaubnutau_1.figure = 'Fig.7'
T4bnutaubnutau_1.figureUrl ='https://atlas.web.cern.ch/Atlas/GROUPS/PHYSICS/PAPERS/SUSY-2016-19/fig_07.png'
T4bnutaubnutau_1.setSources(dataLabels= ['expExclusion', 'obsExclusion', 'upperLimits'],
                 dataFiles= ['orig/Exp_Excl.csv', 'orig/Obs_Excl.csv', 'orig/Obs_UL_fixed.csv'],
                 units = [ None, None, 'fb' ],
                 dataFormats= ['csv', 'csv', 'csv'])

databaseCreator.create()
