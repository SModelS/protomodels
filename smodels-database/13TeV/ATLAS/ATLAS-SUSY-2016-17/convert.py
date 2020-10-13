#!/usr/bin/env python3

"""
.. module:: convert
   :synopsis: used to create info.txt and the <txname>.txt files.

"""
import sys
import os
import argparse
import types

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

databaseCreator.ncpus = 1

#+++++++ global info block ++++++++++++++
info = MetaInfoInput('ATLAS-SUSY-2016-17')
info.url = 'http://atlas.web.cern.ch/Atlas/GROUPS/PHYSICS/PAPERS/SUSY-2016-17/'
info.sqrts = 13
info.lumi = 36.1
info.prettyName = '2 opposite sign leptons + Etmiss'
info.private = False
info.arxiv = 'https://arxiv.org/abs/1708.03247'
info.contact = 'atlas-phys-susy-conveners@cern.ch'
info.publication = 'https://link.springer.com/article/10.1140/epjc/s10052-017-5445-x'
info.comment = ''


#+++++++ dataset block ++++++++++++++
dataset = DataSetInput('data')
dataset.setInfo(dataType = 'upperLimit', dataId = None)
 
#+++++++txName block++++++++++++++++++++
T2tt=dataset.addTxName('T2tt')
T2tt.checked=''
T2tt.constraint="[[['t']],[['t']]]"
T2tt.condition=None
T2tt.conditionDescription = None
T2tt.source="ATLAS"
T2tt.massConstraint=[['dm>=169.0'],['dm>=169.0']]

T2ttoff=dataset.addTxName('T2ttoff')
T2ttoff.checked=''
T2ttoff.constraint="[[['b','W']],[['b','W']]]"
T2ttoff.condition=None
T2ttoff.source="ATLAS"
T2ttoff.massConstraint=[['80 < dm < 169.0'], ['80 < dm < 169.0']]
 
T2bbWWoff=dataset.addTxName('T2bbWWoff')
T2bbWWoff.checked=''
T2bbWWoff.constraint="[[['b','l','nu']],[['b','l','nu']]]"
T2bbWWoff.condition=None
T2bbWWoff.source="CMS"
T2bbWWoff.massConstraint=[['dm<80'], ['dm < 80']]


#++++++next mass plane block+++++++++

T2tt_1 = T2tt.addMassPlane(2*[[x,y]])
T2tt_1.figure='Fig. 8'
T2tt_1.figureUrl='https://atlas.web.cern.ch/Atlas/GROUPS/PHYSICS/PAPERS/SUSY-2016-17/figaux_08.png'
T2tt_1.dataUrl='http://dx.doi.org/10.17182/hepdata.78219.v1/t59'
T2tt_1.histoDataUrl='http://dx.doi.org/10.17182/hepdata.78219.v1/t59'
T2tt_1.exclusionDataUrl='http://dx.doi.org/10.17182/hepdata.78219.v1/t33'
T2tt_1.setSources(dataLabels=['expExclusion', 'obsExclusion', 'upperLimits'],
                 dataFiles= [ 'orig/exclusion_exp_T2tt.csv', 'orig/exclusion_obs_T2tt.csv', 'orig/CrosssectionUpperLimits1.csv'],
                 dataFormats= ['csv', 'csv', 'csv'], objectNames= [ None, None, '$\sigma$ [FB]'
 ], units = [ None , None, 'fb' ] )
T2bbWWoff_1 = T2bbWWoff.addMassPlane(2*[[x,y]])
T2bbWWoff_1.figure='Fig. 8'
T2bbWWoff_1.figureUrl='hbbWWoffps://atlas.web.cern.ch/Atlas/GROUPS/PHYSICS/PAPERS/SUSY-2016-17/figaux_08.png'
T2bbWWoff_1.dataUrl='hbbWWoffp://dx.doi.org/10.17182/hepdata.78219.v1/t59'
T2bbWWoff_1.histoDataUrl='hbbWWoffp://dx.doi.org/10.17182/hepdata.78219.v1/t59'
T2bbWWoff_1.exclusionDataUrl='hbbWWoffp://dx.doi.org/10.17182/hepdata.78219.v1/t33'
T2bbWWoff_1.setSources(dataLabels=['expExclusion', 'obsExclusion', 'upperLimits'],
                 dataFiles= [ 'orig/exclusion_exp_T2bbWWoff.csv', 'orig/exclusion_obs_T2bbWWoff.csv', 'orig/CrosssectionUpperLimits1.csv'],
                 dataFormats= ['csv', 'csv', 'csv'], objectNames= [ None, None, '$\sigma$ [FB]'
 ], units = [ None , None, 'fb' ] )
T2ttoff.addMassPlane(T2bbWWoff_1)
T2bbWWoff.addMassPlane(T2bbWWoff_1)

#+++++++txName block++++++++++++++++++++
T6bbWW=dataset.addTxName('T6bbWW')
T6bbWW.checked=''
T6bbWW.constraint="[[['b'],['W']],[['b'],['W']]]"
T6bbWW.condition=None
T6bbWW.conditionDescription = None
T6bbWW.source="ATLAS"
T6bbWW.massConstraint=[['dm>0.0','dm>=76'],['dm>0.0','dm>=76.0']]


T6bbWWoff=dataset.addTxName('T6bbWWoff')
T6bbWWoff.checked=''
T6bbWWoff.constraint="[[['b'],['l','nu']],[['b'],['l','nu']]]"
T6bbWWoff.condition=None
T6bbWWoff.source="ATLAS"
T6bbWWoff.massConstraint=[['dm>0.0','dm<76.0'],['dm>0.0','dm<76.0']]

#++++++next mass plane block+++++++++

T6bbWW_1 = T6bbWW.addMassPlane(2*[[x,x-10.0,y]])
T6bbWW_1.figure='Fig. 9'
T6bbWW_1.figureUrl='https://atlas.web.cern.ch/Atlas/GROUPS/PHYSICS/PAPERS/SUSY-2016-17/figaux_09.png'
T6bbWW_1.dataUrl='http://dx.doi.org/10.17182/hepdata.78219.v1/t60'
T6bbWW_1.histoDataUrl='http://dx.doi.org/10.17182/hepdata.78219.v1/t60'
T6bbWW_1.exclusionDataUrl='http://dx.doi.org/10.17182/hepdata.78219.v1/t34'
T6bbWW_1.setSources(dataLabels=['expExclusion', 'obsExclusion', 'upperLimits'],
                 dataFiles= [ 'orig/Exclusioncontour(exp)2.csv', 'orig/Exclusioncontour(obs)2.csv', 'orig/CrosssectionUpperLimits2.csv'],
                 dataFormats= ['csv', 'csv', 'csv'], objectNames= [ None, None, 'upperLimits[fb]'], units = [ None, None, 'fb' ] )
T6bbWWoff.addMassPlane(T6bbWW_1)
 

databaseCreator.create()
