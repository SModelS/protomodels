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
info = MetaInfoInput('CMS-SUS-16-009')
info.url = 'https://cms-results.web.cern.ch/cms-results/public-results/publications/SUS-16-009/'
info.sqrts = 13
info.lumi = 2.3
info.prettyName = 'multijets + Etmiss, top tagging'
info.private = False
info.arxiv = ''
info.contact = 'cms-phys-conveners-sus@cern.ch'
info.publication = 'Phys. Rev. D 96 (2017), 012004'
info.comment = ''



#+++++++ dataset block ++++++++++++++
dataset = DataSetInput('data')
dataset.setInfo(dataType = 'upperLimit', dataId = None)

#+++++++ next txName block ++++++++++++++
T2tt = dataset.addTxName('T2tt')
T2tt.checked = ''
T2tt.constraint = "[[['t']],[['t']]]"
T2tt.conditionDescription = None
T2tt.condition = None
T2tt.source = "CMS"

#+++++++ next txName block ++++++++++++++
T2ttoff=dataset.addTxName('T2ttoff')
T2ttoff.checked=''
T2ttoff.constraint="[[['b','W']],[['b','W']]]"
T2ttoff.condition=None
T2ttoff.conditionDescription=None
T2ttoff.massConstraint=[['80<dm<169'],['80<dm<169']]
T2ttoff.source="CMS"
#+++++++ next mass plane block ++++++++++++++
T2tt_1 = T2tt.addMassPlane(2*[[x, y]])
T2tt_1.figure = 'Fig.12-a'
T2tt_1.figureUrl = 'http://cms-results.web.cern.ch/cms-results/public-results/publications/SUS-16-009/CMS-SUS-16-009_Figure_012-a.png'
T2tt_1.dataUrl = 'http://cms-results.web.cern.ch/cms-results/public-results/publications/SUS-16-009/CMS-SUS-16-009_Figure_012-a.root'
T2tt_1.histoDataUrl = 'http://cms-results.web.cern.ch/cms-results/public-results/publications/SUS-16-009/CMS-SUS-16-009_Figure_012-a.root'
T2tt_1.dataUrl = 'http://cms-results.web.cern.ch/cms-results/public-results/publications/SUS-16-009/CMS-SUS-16-009_Figure_012-a.root'
T2tt_1.exclusionDataUrl = 'http://cms-results.web.cern.ch/cms-results/public-results/publications/SUS-15-002/CMS-SUS-16-009_Figure_012-a.root'
T2tt_1.setSources(dataLabels= ['upperLimits','expectedUpperLimits','obsExclusion','expExclusion'],
                 dataFiles= ['orig/CMS-SUS-16-009_Figure_012-a.root', 'orig/CMS-SUS-16-009_Figure_012-a.root','orig/CMS-SUS-16-009_Figure_012-aObs.csv','orig/CMS-SUS-16-009_Figure_012-aExp.csv'],
                 dataFormats= ['root', 'root','csv','csv'],objectNames= ['combined_obsLimit_BR100pct', 'combined_expLimit_BR100pct', None, None ])
                 
T2ttoff.addMassPlane(T2tt_1)                 
                 

                 
                 
#+++++++ next txName block ++++++++++++++
T1tttt = dataset.addTxName('T1tttt')
T1tttt.checked = ''
T1tttt.constraint =  "[[['t','t']],[['t','t']]]"
T1tttt.conditionDescription = None
T1tttt.condition = None
T1tttt.source = "CMS"
T1tttt.massConstraint = None
T1ttttoff = dataset.addTxName('T1ttttoff')
T1ttttoff.checked = ''
T1ttttoff.constraint = "[[['b','W','b','W']],[['b','W','b','W']]]"
T1ttttoff.conditionDescription = None
T1ttttoff.condition = None
T1ttttoff.massConstraint = [['dm <= 338.0'], ['dm <= 338.0']]
T1ttttoff.source = "CMS"
#+++++++ next mass plane block ++++++++++++++
T1tttt_1 = T1tttt.addMassPlane(2*[[x, y]])
T1tttt_1.figure = 'Fig.13-a'
T1tttt_1.figureUrl = 'http://cms-results.web.cern.ch/cms-results/public-results/publications/SUS-16-009/CMS-SUS-16-009_Figure_013-a.png'
T1tttt_1.dataUrl = 'http://cms-results.web.cern.ch/cms-results/public-results/publications/SUS-16-009/CMS-SUS-16-009_Figure_013-a.root'
T1tttt_1.histoDataUrl = 'http://cms-results.web.cern.ch/cms-results/public-results/publications/SUS-16-009/CMS-SUS-16-009_Figure_013-a.root'
T1tttt_1.dataUrl = 'http://cms-results.web.cern.ch/cms-results/public-results/publications/SUS-16-009/CMS-SUS-16-009_Figure_013-a.root'
T1tttt_1.exclusionDataUrl = 'http://cms-results.web.cern.ch/cms-results/public-results/publications/SUS-16-009/CMS-SUS-16-009_Figure_013-a.root'
T1tttt_1.setSources(dataLabels= ['upperLimits','expectedUpperLimits','obsExclusion','expExclusion'],
                 dataFiles= ['orig/CMS-SUS-16-009_Figure_013-a.root', 'orig/CMS-SUS-16-009_Figure_013-a.root','orig/CMS-SUS-16-009_Figure_013-aObs.csv','orig/CMS-SUS-16-009_Figure_013-aExp.csv'],
                 dataFormats= ['root', 'root','csv','csv'],objectNames= ['combined_obsLimit_BR100pct', 'combined_expLimit_BR100pct', None, None ])
T1ttttoff.addMassPlane(T1tttt_1)

#+++++++ next txName block ++++++++++++++
T5tctc = dataset.addTxName('T5tctc')
T5tctc.checked =''
T5tctc.constraint = "[[['t'],['c']],[['t'],['c']]]"
T5tctc.conditionDescription = None
T5tctc.condition = None
T5tctc.source = "CMS"
#+++++++ next mass plane block ++++++++++++++
T5tctc_1 = T5tctc.addMassPlane(2*[[x,y+20, y]])
T5tctc_1.figure = 'Fig.13-b'
T5tctc_1.figureUrl = 'http://cms-results.web.cern.ch/cms-results/public-results/publications/SUS-16-009/CMS-SUS-16-009_Figure_013-b.png'
T5tctc_1.dataUrl = 'http://cms-results.web.cern.ch/cms-results/public-results/publications/SUS-16-009/CMS-SUS-16-009_Figure_013-b.root'
T5tctc_1.histoDataUrl = 'http://cms-results.web.cern.ch/cms-results/public-results/publications/SUS-16-009/CMS-SUS-16-009_Figure_013-b.root'
T5tctc_1.dataUrl = 'http://cms-results.web.cern.ch/cms-results/public-results/publications/SUS-15-002/CMS-SUS-16-009_Figure_013-b.root'
T5tctc_1.exclusionDataUrl = 'http://cms-results.web.cern.ch/cms-results/public-results/publications/SUS-16-009/CMS-SUS-16-009_Figure_013-b.root'
T5tctc_1.setSources(dataLabels= ['upperLimits','expectedUpperLimits','obsExclusion','expExclusion'],
                 dataFiles= ['orig/CMS-SUS-16-009_Figure_013-b.root', 'orig/CMS-SUS-16-009_Figure_013-b.root','orig/CMS-SUS-16-009_Figure_013-bObs.csv','orig/CMS-SUS-16-009_Figure_013-bExp.csv'],
                 dataFormats= ['root', 'root','csv','csv'],objectNames= ['combined_obsLimit_BR100pct', 'combined_expLimit_BR100pct', None, None ])



databaseCreator.create()
