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
argparser.add_argument ('-no', '--noUpdate',
help = 'do not update the lastUpdate field.',\
action= "store_true" )
argparser.add_argument ('-r', '--resetValidation',
help = 'reset the validation flag',\
action= "store_true" )
args = argparser.parse_args()

if args.noUpdate:
    os.environ["SMODELS_NOUPDATE"]="1"

if args.resetValidation:
    os.environ["SMODELS_RESETVALIDATION"]="1"

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
info = MetaInfoInput('CMS-SUS-17-006')
info.url = 'https://cms-results.web.cern.ch/cms-results/public-results/publications/SUS-17-006/'
info.sqrts = 13
info.lumi = 35.9
info.prettyName = 'jets + boosted H(bb) + Etmiss'
info.private = False
info.arxiv = 'https://arxiv.org/abs/1712.08501'
info.contact = 'cms-phys-conveners-sus@cern.ch'
info.publication = 'Phys. Rev. Lett.  120, no. 24, 241801 (2018)'
info.comment = ''



#+++++++ dataset block ++++++++++++++
dataset = DataSetInput('data')
dataset.setInfo(dataType = 'upperLimit', dataId = None)


lsp_masses = [0., 1.]
#+++++txName block +++++++++++++++++
T5HH=dataset.addTxName('T5HH')
T5HH.checked=''
T5HH.constraint="[[['jet','jet'],['higgs']],[['jet','jet'],['higgs']]]"
T5HH.condition=None
T5HH.conditionDescription = None
T5HH.source="CMS"
#T5HH.massConstraint=None




#+++++++ next mass plane block ++++++++++++++
for lsp in lsp_masses:
	plane 						= T5HH.addMassPlane(2*[[x,x-50,lsp]])	#2*[[x, 1.]])#[[x],[x]])#2*[[x]]
	plane.figure 				= 'Fig. 3'
	plane.figureUrl 			= 'http://cms-results.web.cern.ch/cms-results/public-results/publications/SUS-17-006/CMS-SUS-17-006_Figure_003.png'
	plane.dataUrl 				= 'http://cms-results.web.cern.ch/cms-results/public-results/publications/SUS-17-006/CMS-SUS-17-006_Figure_003.root'
	plane.setSources(dataLabels=['expectedUpperLimits','upperLimits'],
                    dataFiles=['orig/CMS-SUS-17-006_Figure_003.root','orig/CMS-SUS-17-006_Figure_003.root'],
                    dataFormats=['root','root',],objectNames=['ExpectedLimitT5qqqqHH;1','ObservedLimitT5qqqqHH;1'],
                    units=['pb','pb'])


lsp_masses = [0., 1.]
 #+++++txName block +++++++++++++++++
T5HZ=dataset.addTxName('T5HZ')
T5HZ.checked=''
#T5HZ.constraint="[[['jet','jet'],['higgs']],[['jet','jet'],['Z']]]"
T5HZ.constraint="[[['jet','jet'],['higgs']],[['jet','jet'],['Z']]]+[[['jet','jet'],['higgs']],[['jet','jet'],['higgs']]]+[[['jet','jet'],['Z']],[['jet','jet'],['Z']]]"
T5HZ.condition="Csim([[['jet','jet'],['higgs']],[['jet','jet'],['Z']]],2.0*([[['jet','jet'],['higgs']],[['jet','jet'],['higgs']]]),2.0*([[['jet','jet'],['Z']],[['jet','jet'],['Z']]]))"
T5HZ.conditionDescription = None
T5HZ.source="CMS"
#T5HH.massConstraint=None



#+++++++ next mass plane block ++++++++++++++
for lsp in lsp_masses:
	plane 						= T5HZ.addMassPlane(2*[[x,x-50,lsp]])	#2*[[x, 1.]])#[[x],[x]])#2*[[x]]
	plane.figure 				= 'Fig. 3'
	plane.figureUrl 			= 'http://cms-results.web.cern.ch/cms-results/public-results/publications/SUS-17-006/CMS-SUS-17-006_Figure_003.png'
	plane.dataUrl 				= 'http://cms-results.web.cern.ch/cms-results/public-results/publications/SUS-17-006/CMS-SUS-17-006_Figure_003.root'
	plane.setSources(dataLabels=['expectedUpperLimits','upperLimits'],
                    dataFiles=['orig/CMS-SUS-17-006_Figure_003.root','orig/CMS-SUS-17-006_Figure_003.root'],
                    dataFormats=['root','root',],objectNames=['ExpectedLimitT5qqqqZH;1','ObservedLimitT5qqqqZH;1'],
                    units=['pb','pb'])


databaseCreator.create()
