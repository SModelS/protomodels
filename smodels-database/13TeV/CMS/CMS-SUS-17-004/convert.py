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
argparser.add_argument ('-utilsPath', '-utilsPath', 
help = 'path to the package smodels_utils',\
type = str )
argparser.add_argument ('-smodelsPath', '-smodelsPath', 
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
info = MetaInfoInput('CMS-SUS-17-004')
info.url = 'http://cms-results.web.cern.ch/cms-results/public-results/publications/SUS-17-004/index.html'
info.sqrts = 13
info.lumi = 35.9
info.prettyName = 'EW-ino combination'
info.private = False
info.arxiv = 'https://arxiv.org/abs/1801.03957'
info.contact = 'cms-phys-conveners-sus@cern.ch'
info.publication = 'JHEP 1803 (2018) 160, http://dx.doi:10.1007/JHEP03(2018)160'
info.comment = 'Uses Fig. 8a for WZ topology instead of Fig. 7'
info.implementedBy = 'Sabine'

#+++++++ dataset block ++++++++++++++
dataset = DataSetInput('data')
dataset.setInfo(dataType = 'upperLimit', dataId = None)

#+++++txName block +++++++++++++++++

TChiWZ=dataset.addTxName('TChiWZ')
TChiWZ.checked=''
TChiWZ.constraint="[[['W']],[['Z']]]"
TChiWZ.condition=None
TChiWZ.conditionDescription = None
#TChiWZ.massConstraint=[['dm > 169.']]*2  #Use only on-shell region to avoid interpolating in the excluded band
TChiWZ.source="CMS"


#+++++txName block +++++++++++++++++
TChiWZoff=dataset.addTxName('TChiWZoff')
TChiWZoff.checked=''
TChiWZoff.constraint="71.*([[['mu+','mu-']],[['l','nu']]] + [[['e+','e-']],[['l','nu']]])"
TChiWZoff.condition = "cGtr([[['mu+','mu-']],[['l','nu']]],[[['e+','e-']],[['l','nu']]])"
TChiWZoff.massConstraint = [['dm < 86.0'], ['dm < 76.0']]
TChiWZoff.conditionDescription=None
TChiWZoff.source="CMS"


#++++++next mass plane block+++++++++

TChiWZ_1 = TChiWZ.addMassPlane(2*[[x,y]])
TChiWZ_1.figure='Fig. 8'
TChiWZ_1.figureUrl='http://cms-results.web.cern.ch/cms-results/public-results/publications/SUS-17-004/CMS-SUS-17-004_Figure_008-a.png'
TChiWZ_1.dataUrl='http://cms-results.web.cern.ch/cms-results/public-results/publications/SUS-17-004/CMS-SUS-17-004_Figure_008-a.root'
TChiWZ_1.histoDataUrl='http://cms-results.web.cern.ch/cms-results/public-results/publications/SUS-17-004/CMS-SUS-17-004_Figure_008-a.root'
TChiWZ_1.exclusionDataUrl='http://cms-results.web.cern.ch/cms-results/public-results/publications/SUS-17-004/CMS-SUS-17-004_Figure_008-a.root'
TChiWZ_1.setSources(dataLabels=['expExclusion','expExclusionM1','expExclusionP1','obsExclusion','obsExclusionM1','obsExclusionP1','upperLimits'],
                    dataFiles=['orig/CMS-SUS-17-004_Figure_008-a.root','orig/CMS-SUS-17-004_Figure_008-a.root',
                    'orig/CMS-SUS-17-004_Figure_008-a.root','orig/CMS-SUS-17-004_Figure_008-a.root','orig/CMS-SUS-17-004_Figure_008-a.root',
                    'orig/CMS-SUS-17-004_Figure_008-a.root','orig/CMS-SUS-17-004_Figure_008-a.root'],
indices=[4, 6, 5, 7, 9, 8, 2],dataFormats=['canvas','canvas','canvas','canvas','canvas','canvas','canvas'],objectNames=['TChiWZ;1','TChiWZ;1','TChiWZ;1',
'TChiWZ;1','TChiWZ;1','TChiWZ;1','TChiWZ;1'],units=[None,None,None,None,None,None,'pb'])
TChiWZoff.addMassPlane(TChiWZ_1)
 
 
 
 
 
 #++++++next txName block+++++++++++++++

TChiWH=dataset.addTxName('TChiWH')
TChiWH.checked=''
TChiWH.constraint="[[['W']],[['higgs']]]"
TChiWH.condition=None
TChiWH.conditionDescription=None
TChiWH.source="CMS"


#++++++next mass plane block++++++++
TChiWH_1 = TChiWH.addMassPlane(2*[[x,y]])
TChiWH_1.figure='Fig. 8-b'
TChiWH_1.figureUrl='http://cms-results.web.cern.ch/cms-results/public-results/publications/SUS-17-004/CMS-SUS-17-004_Figure_008-b.png'
TChiWH_1.dataUrl='http://cms-results.web.cern.ch/cms-results/public-results/publications/SUS-17-004/CMS-SUS-17-004_Figure_008-b.root'
TChiWH_1.histoDataUrl='http://cms-results.web.cern.ch/cms-results/public-results/publications/SUS-17-004/CMS-SUS-17-004_Figure_008-b.root'
TChiWH_1.exclusionDataUrl='http://cms-results.web.cern.ch/cms-results/public-results/publications/SUS-17-004/CMS-SUS-17-004_Figure_008-b.root'
TChiWH_1.setSources(dataLabels=['expExclusion','expExclusionM1','expExclusionP1','obsExclusion','obsExclusionM1','obsExclusionP1','upperLimits'],
                    dataFiles=['orig/CMS-SUS-17-004_Figure_008-b.root','orig/CMS-SUS-17-004_Figure_008-b.root','orig/CMS-SUS-17-004_Figure_008-b.root','orig/CMS-SUS-17-004_Figure_008-b.root',
                    'orig/CMS-SUS-17-004_Figure_008-b.root','orig/CMS-SUS-17-004_Figure_008-b.root','orig/CMS-SUS-17-004_Figure_008-b.root'],
                    dataFormats=['canvas','canvas','canvas','canvas','canvas','canvas','canvas'],objectNames=['TChiWH;1','TChiWH;1','TChiWH;1','TChiWH;1','TChiWH;1','TChiWH;1','TChiWH;1'],
                    indices= [4, 6, 5, 7, 9, 8, 2],units=[None,None,None,None,None,None,'pb'])
 
 
 
 
 
 
 
 

databaseCreator.create()



