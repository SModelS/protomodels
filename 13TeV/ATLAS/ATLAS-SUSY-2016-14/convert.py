#!/usr/bin/env python

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
type = types.StringType)
argparser.add_argument ('-smodelsPath', '--smodelsPath', 
help = 'path to the package smodels_utils',\
type = types.StringType)
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
info = MetaInfoInput('ATLAS-SUSY-2016-14')
info.url = 'http://atlas.web.cern.ch/Atlas/GROUPS/PHYSICS/PAPERS/SUSY-2016-14/'
info.sqrts = 13
info.lumi = 36.1
info.prettyName = '2 same-sign or 3 leptons + jets + Etmiss'
info.private = False
info.arxiv = 'https://arxiv.org/abs/1706.03731'
info.contact = 'atlas-phys-susy-conveners@cern.ch'
info.publication = 'https://link.springer.com/article/10.1007/JHEP09(2017)084'
info.comment = 'Moriond 2017. Omitted RPV SUSY and long cascade topologies.'


#+++++++ dataset block ++++++++++++++
dataset = DataSetInput('data')
dataset.setInfo(dataType = 'upperLimit', dataId = None)

#+++++++txName block++++++++++++++++++++
T1tttt=dataset.addTxName('T1tttt')
T1tttt.checked=''
T1tttt.constraint="[[['t','t']],[['t','t']]]"
T1tttt.condition=None
T1tttt.conditionDescription = None
T1tttt.source="ATLAS"
T1tttt.massConstraint=[['dm>=338.0'],['dm>=338.0']]

T1ttttoff=dataset.addTxName('T1ttttoff')
T1ttttoff.checked=''
T1ttttoff.constraint="[[['b','W','b','W']],[['b','W','b','W']]]"
T1ttttoff.condition=None
T1ttttoff.source="ATLAS"
T1ttttoff.massConstraint=[['160 < dm < 338.0'], ['160 < dm < 338.0']]

T1ttofftt=dataset.addTxName('T1ttofftt')
T1ttofftt.checked=''
T1ttofftt.constraint="[[['t','b','W']],[['t','b','W']]]"
T1ttofftt.condition=None
T1ttofftt.source="CMS"
T1ttofftt.massConstraint=[['dm>169'], ['dm>169']]

#++++++next mass plane block+++++++++

T1tttt_1 = T1tttt.addMassPlane(2*[[x,y]])
T1tttt_1.figure='Fig. 12-a'
T1tttt_1.figureUrl='http://atlas.web.cern.ch/Atlas/GROUPS/PHYSICS/PAPERS/SUSY-2016-14/figaux_12a.png'
T1tttt_1.dataUrl='http://dx.doi.org/10.17182/hepdata.77719.v1/t17'
T1tttt_1.histoDataUrl='http://dx.doi.org/10.17182/hepdata.77719.v1/t17'
T1tttt_1.exclusionDataUrl='http://dx.doi.org/10.17182/hepdata.77719.v1/t1'
T1tttt_1.setSources(dataLabels=['expExclusion', 'obsExclusion', 'upperLimits'],
                 dataFiles= [ 'orig/Exclusioncontour1(exp.).csv', 'orig/Exclusioncontour1(obs.).csv', 'orig/X-sectionU.L.,acc.,eff.1.csv'],
                 dataFormats= ['csv', 'csv', 'csv'], objectNames= [ None, None, 'SIG 95%CL [FB]' ], units = [ None, None, 'fb' ] )
T1ttttoff.addMassPlane(T1tttt_1)
T1ttofftt.addMassPlane(T1tttt_1)

#+++++++txName block++++++++++++++++++++
T6ttWW=dataset.addTxName('T6ttWW')
T6ttWW.checked=''
T6ttWW.constraint="[[['t'],['W']],[['t'],['W']]]"
T6ttWW.condition=None
T6ttWW.conditionDescription = None
T6ttWW.source="ATLAS"
T6ttWW.massConstraint=[['dm>=169.0',' dm>=76'],['dm>=169.0',' dm>=76']]


T6ttWWoff=dataset.addTxName('T6ttWWoff')
T6ttWWoff.checked=''
T6ttWWoff.constraint="[[['b','W'],['W']],[['b','W'],['W']]]"
T6ttWWoff.condition=None
T6ttWWoff.source="ATLAS"
T6ttWWoff.massConstraint=[['80<dm<169.0',' dm>76'],['80<dm<169.0',' dm>76']]

#++++++next mass plane block+++++++++

T6ttWW_1 = T6ttWW.addMassPlane(2*[[x,y+100.0,y]])
T6ttWW_1.figure='Fig. 12-d'
T6ttWW_1.figureUrl='http://atlas.web.cern.ch/Atlas/GROUPS/PHYSICS/PAPERS/SUSY-2016-14/figaux_12d.png'
T6ttWW_1.dataUrl='http://dx.doi.org/10.17182/hepdata.77719.v1/t20'
T6ttWW_1.histoDataUrl='http://dx.doi.org/10.17182/hepdata.77719.v1/t20'
T6ttWW_1.exclusionDataUrl='http://dx.doi.org/10.17182/hepdata.77719.v1/t8'
T6ttWW_1.setSources(dataLabels=['expExclusion', 'obsExclusion', 'upperLimits'],
                 dataFiles= [ 'orig/Exclusioncontour4(exp.).csv', 'orig/Exclusioncontour4(obs.).csv', 'orig/X-sectionU.L.,acc.,eff.4.csv'],
                 dataFormats= ['csv', 'csv', 'csv'], objectNames= [ None, None, 'SIG 95%CL [FB]' ], units = [ None, None, 'fb' ] )
T6ttWWoff.addMassPlane(T6ttWW_1)


databaseCreator.create()
