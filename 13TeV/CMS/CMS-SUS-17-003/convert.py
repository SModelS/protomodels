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
help = 'path to smodels',\
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

databaseCreator.ncpus = 1


#+++++++ global info block ++++++++++++++
info = MetaInfoInput('CMS-SUS-17-003')
info.url = 'https://cms-results.web.cern.ch/cms-results/public-results/publications/SUS-17-003/'
info.sqrts = 13
info.lumi = 35.9
info.prettyName = '2 taus + Etmiss'
info.private = False
info.arxiv = 'https://arxiv.org/abs/1807.02048'
info.contact = 'cms-phys-conveners-sus@cern.ch'
info.publication = 'JHEP 1811 (2018) 151'
info.comment = 'the direct stau pair production followed by each stau decaying to a tau lepton and neutralino, has not been included because corresponding higher luminosity ATLAS result is already implemented'


#+++++++ dataset block ++++++++++++++
dataset = DataSetInput('data')
dataset.setInfo(dataType = 'upperLimit', dataId = None)



#+++++next txName block+++++++++++++++
TChiChipmStauStau=dataset.addTxName('TChiChipmStauStau')
TChiChipmStauStau.checked=''
TChiChipmStauStau.constraint="[[['ta+'],['ta-']],[['nu'],['ta']]]+[[['ta-'],['ta+']],[['nu'],['ta']]]"
TChiChipmStauStau.condition=None
TChiChipmStauStau.conditionDescription=None
TChiChipmStauStau.source = "CMS"
TChiChipmStauStau.massConstraint=None

#++++++next txName block+++++++++++++++

TChiChipmStauStau_1=TChiChipmStauStau.addMassPlane(2*[[x,0.5*x+0.5*y,y]])
TChiChipmStauStau_1.figure='Fig. 14'
TChiChipmStauStau_1.figureUrl='https://cms-results.web.cern.ch/cms-results/public-results/publications/SUS-17-003/CMS-SUS-17-003_Figure_014.png'
TChiChipmStauStau_1.dataUrl='https://cms-results.web.cern.ch/cms-results/public-results/publications/SUS-17-003/CMS-SUS-17-003_Figure_014.root'
TChiChipmStauStau_1.histoDataUrl='https://cms-results.web.cern.ch/cms-results/public-results/publications/SUS-17-003/CMS-SUS-17-003_Figure_014.root'
TChiChipmStauStau_1.exclusionDataUrl='https://cms-results.web.cern.ch/cms-results/public-results/publications/SUS-17-003/CMS-SUS-17-003_Figure_014.root'
TChiChipmStauStau_1.setSources(dataLabels=['expExclusionP1','expExclusionM1','expExclusion','obsExclusionM1','obsExclusionP1','obsExclusion','upperLimits'],
                    dataFiles=['orig/CMS-SUS-17-003_Figure_014.root','orig/CMS-SUS-17-003_Figure_014.root','orig/CMS-SUS-17-003_Figure_014.root','orig/CMS-SUS-17-003_Figure_014.root','orig/CMS-SUS-17-003_Figure_014.root','orig/CMS-SUS-17-003_Figure_014.root','orig/CMS-SUS-17-003_Figure_014.root'], dataFormats=['root','root','root','root','root','root','root'],
                        objectNames=['graph_smoothed_ExpP;1','graph_smoothed_ExpM;1','graph_smoothed_Exp;1','graph_smoothed_ObsM;1','graph_smoothed_ObsP;1','graph_smoothed_Obs;1','hXsec_exp_corr;1'],
                    indices= [1,2,3, 6, 7, 8, 9],units=[None,None,None,None,None,None,'pb'])
 
 

#+++++next txName block+++++++++++++++
TChipChimStauSnu=dataset.addTxName('TChipChimStauSnu')
TChipChimStauSnu.checked=''
#TChipChimStauSnu.constraint="[[[nu],[ta+]],[[nu],[ta-]]]+[[[ta+],[nu]],[[ta-],[nu]]]"
TChipChimStauSnu.constraint="[[[nu],[ta+]],[[ta-],[nu]]]+[[[nu],[ta-]],[[ta+],[nu]]]+[[[nu],[ta+]],[[nu],[ta-]]]+[[[ta+],[nu]],[[ta-],[nu]]]"
TChipChimStauSnu.condition=None
TChipChimStauSnu.conditionDescription=None
TChipChimStauSnu.source = "CMS"
TChipChimStauSnu.massConstraint=None

#++++++next txName block+++++++++++++++

TChipChimStauSnu_1=TChipChimStauSnu.addMassPlane(2*[[x,0.5*x+0.5*y,y]])
TChipChimStauSnu_1.figure='Fig. 15'
TChipChimStauSnu_1.figureUrl='https://cms-results.web.cern.ch/cms-results/public-results/publications/SUS-17-003/CMS-SUS-17-003_Figure_015.png'
TChipChimStauSnu_1.dataUrl='https://cms-results.web.cern.ch/cms-results/public-results/publications/SUS-17-003/CMS-SUS-17-003_Figure_015.root'
TChipChimStauSnu_1.histoDataUrl='https://cms-results.web.cern.ch/cms-results/public-results/publications/SUS-17-003/CMS-SUS-17-003_Figure_015.root'
TChipChimStauSnu_1.exclusionDataUrl='https://cms-results.web.cern.ch/cms-results/public-results/publications/SUS-17-003/CMS-SUS-17-003_Figure_015.root'
TChipChimStauSnu_1.setSources(dataLabels=['expExclusionP1','expExclusionM1','expExclusion','obsExclusionM1','obsExclusionP1','obsExclusion','upperLimits'],
                    dataFiles=['orig/CMS-SUS-17-003_Figure_015.root','orig/CMS-SUS-17-003_Figure_015.root','orig/CMS-SUS-17-003_Figure_015.root','orig/CMS-SUS-17-003_Figure_015.root','orig/CMS-SUS-17-003_Figure_015.root','orig/CMS-SUS-17-003_Figure_015.root','orig/CMS-SUS-17-003_Figure_015.root'], dataFormats=['root','root','root','root','root','root','root'],
                        objectNames=['graph_smoothed_ExpP;1','graph_smoothed_ExpM;1','graph_smoothed_Exp;1','graph_smoothed_ObsM;1','graph_smoothed_ObsP;1','graph_smoothed_Obs;1','hXsec_exp_corr;1'],
                    indices= [5,6,4,3,2,1,9],units=[None,None,None,None,None,None,'pb'])

	
	
databaseCreator.create()
