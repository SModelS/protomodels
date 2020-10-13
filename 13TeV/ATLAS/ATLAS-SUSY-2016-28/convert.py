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
info = MetaInfoInput('ATLAS-SUSY-2016-28')
info.url = "https://atlas.web.cern.ch/Atlas/GROUPS/PHYSICS/PAPERS/SUSY-2016-28/"
info.sqrts = 13
info.lumi = 36.1
info.prettyName = "2 b-jets + Etmiss"
info.private = False
info.arxiv =  'https://arxiv.org/abs/1708.09266'
info.contact = 'atlas-phys-susy-conveners@cern.ch'
info.publication ='https://link.springer.com/article/10.1007%2FJHEP11%282017%29195'


#+++++++ dataset block ++++++++++++++
dataset = DataSetInput('data')
dataset.setInfo(dataType = 'upperLimit', dataId = None)

#+++++++ next txName block ++++++++++++++
T2bb = dataset.addTxName('T2bb')
T2bb.checked = ''
T2bb.constraint ="[[['b']],[['b']]]"
T2bb.conditionDescription = None
T2bb.condition = None
T2bb.source = "ATLAS"
#+++++++ next mass plane block ++++++++++++++
T2bb_1 = T2bb.addMassPlane(2*[[x, y]])
T2bb_1.figure = 'Fig.7a'
T2bb_1.figureUrl ='https://atlas.web.cern.ch/Atlas/GROUPS/PHYSICS/PAPERS/SUSY-2016-28/figaux_08a.png'
T2bb_1.dataUrl = 'https://www.hepdata.net/record/79165'
T2bb_1.exclusionDataUrl = 'https://www.hepdata.net/record/79165?version=1&table=Limitoncrosssection1'
T2bb_1.setSources(dataLabels= ['expExclusion', 'obsExclusion', 'upperLimits'],
                 dataFiles= ['orig/Contour17.csv' ,'orig/Contour18.csv', 'orig/Limitoncrosssection1.csv'],
                 dataFormats= ['csv', 'csv', 'csv'])



databaseCreator.create()
