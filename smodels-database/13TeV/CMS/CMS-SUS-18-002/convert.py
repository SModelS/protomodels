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
info = MetaInfoInput('CMS-SUS-18-002')
info.url = 'https://cms-results.web.cern.ch/cms-results/public-results/publications/SUS-18-002/'
info.sqrts = 13
info.lumi = 35.9
info.prettyName = 'photon, jets, b-jets+ Etmiss, top tagging'
info.private = False
info.arxiv = 'https://arxiv.org/abs/1901.06726'
info.contact = 'cms-phys-conveners-sus@cern.ch'
info.publication = 'Eur. Phys. J. C 79 (2019) no.5,  444'
info.comment = ''


 


#+++++++ dataset block ++++++++++++++
dataset = DataSetInput('data')
dataset.setInfo(dataType = 'upperLimit', dataId = None)


lsp_masses = [1.]

#+++++++ next txName block ++++++++++++++
T5Hg 					  		= dataset.addTxName('T5Hg')
T5Hg.checked 			  		= 'No'
T5Hg.constraint 		  		= "[[['jet','jet'],['higgs']],[['jet','jet'],['photon']]]+[[['jet','jet'],['photon']],[['jet','jet'],['photon']]]+[[['jet','jet'],['higgs']],[['jet','jet'],['higgs']]]"
T5Hg.conditionDescription 		= None
T5Hg.condition 			  		= "Csim([[['jet','jet'],['higgs']],[['jet','jet'],['photon']]]+2.*[[['jet','jet'],['photon']],[['jet','jet'],['photon']]]+2.*[[['jet','jet'],['higgs']],[['jet','jet'],['higgs']]])"
T5Hg.source 			  		= "CMS"
#+++++++ next mass plane block ++++++++++++++
for lsp in lsp_masses:
	p 						= T5Hg.addMassPlane(2*[[x, y, lsp]])
	p.figure    			= 'Fig.5a'
	p.figureUrl 			= 'http://cms-results.web.cern.ch/cms-results/public-results/publications/SUS-18-002/CMS-SUS-18-002_Figure_005-a.png'
	p.dataUrl   			= "http://cms-results.web.cern.ch/cms-results/public-results/publications/SUS-18-002/CMS-SUS-18-002_Figure_005-a.root"
	p.setSources(dataLabels	= ['expExclusionP1','expExclusionM1','obsExclusionP1','obsExclusionM1','expExclusion', 'obsExclusion', 'upperLimits','expectedUpperLimits'],
				dataFiles	= ['orig/CMS-SUS-18-002_Figure_005-a.root', 'orig/CMS-SUS-18-002_Figure_005-a.root', 'orig/CMS-SUS-18-002_Figure_005-a.root', 'orig/CMS-SUS-18-002_Figure_005-a.root','orig/CMS-SUS-18-002_Figure_005-a.root', 'orig/CMS-SUS-18-002_Figure_005-a.root', 'orig/CMS-SUS-18-002_Figure_005-a.root', 'orig/CMS-SUS-18-002_Figure_005-a.root'],
					units	= [ None, None,None, None,None, None, 'fb','fb' ],
				dataFormats	= ['root', 'root', 'root', 'root','root', 'root', 'root', 'root'],objectNames= ['exp1up;1','exp1dn;1','obs_XsecUp;1','obs_XsecDn;1','exp;1', 'obs;1','obs_XsecLimit;1','exp_XsecLimit;1'],
                    indices= [5,4,7,8,3,6,1,2])






#+++++++ next txName block ++++++++++++++
T5bbbbZg 					  		= dataset.addTxName('T5bbbbZg')
T5bbbbZg.checked 			  		= 'No'
T5bbbbZg.constraint 		  		= "[[['b','b'],['Z']],[['b','b'],['photon']]]+[[['b','b'],['photon']],[['b','b'],['photon']]]+[[['b','b'],['Z']],[['b','b'],['Z']]]"
T5bbbbZg.conditionDescription 		= None
T5bbbbZg.condition 			  		= "Csim([[['jet','jet'],['Z']],[['jet','jet'],['photon']]]+2.*[[['jet','jet'],['photon']],[['jet','jet'],['photon']]]+2.*[[['jet','jet'],['Z']],[['jet','jet'],['Z']]])"
T5bbbbZg.source 			  		= "CMS"
#+++++++ next mass plane block ++++++++++++++
for lsp in lsp_masses:
	p 						= T5bbbbZg.addMassPlane(2*[[x, y, lsp]])
	p.figure    			= 'Fig.5b'
	p.figureUrl 			= 'http://cms-results.web.cern.ch/cms-results/public-results/publications/SUS-18-002/CMS-SUS-18-002_Figure_005-b.png'
	p.dataUrl   			= "http://cms-results.web.cern.ch/cms-results/public-results/publications/SUS-18-002/CMS-SUS-18-002_Figure_005-b.root"
	p.setSources(dataLabels	= ['expExclusionP1','expExclusionM1','obsExclusionP1','obsExclusionM1','expExclusion', 'obsExclusion', 'upperLimits','expectedUpperLimits'],
				dataFiles	= ['orig/CMS-SUS-18-002_Figure_005-b.root', 'orig/CMS-SUS-18-002_Figure_005-b.root', 'orig/CMS-SUS-18-002_Figure_005-b.root', 'orig/CMS-SUS-18-002_Figure_005-b.root','orig/CMS-SUS-18-002_Figure_005-b.root', 'orig/CMS-SUS-18-002_Figure_005-b.root', 'orig/CMS-SUS-18-002_Figure_005-b.root', 'orig/CMS-SUS-18-002_Figure_005-b.root'],
					units	= [ None, None,None, None,None, None, 'fb','fb' ],
				dataFormats	= ['root', 'root', 'root', 'root','root', 'root', 'root', 'root'],objectNames= ['exp1up;1','exp1dn;1','obs_XsecUp;1','obs_XsecDn;1','exp;1', 'obs;1','obs_XsecLimit;1','exp_XsecLimit;1'],
                    indices= [5,4,7,8,3,6,1,2])


#+++++++ next txName block ++++++++++++++
T5ttttZg 					  		= dataset.addTxName('T5ttttZg')
T5ttttZg.checked 			  		= 'No'
T5ttttZg.constraint 		  		= "[[['t','t'],['Z']],[['t','t'],['photon']]]+[[['t','t'],['photon']],[['t','t'],['photon']]]+[[['t','t'],['Z']],[['t','t'],['Z']]]"
T5ttttZg.conditionDescription 		= None
T5ttttZg.condition 			  		= "Csim([[['t','t'],['Z']],[['t','t'],['photon']]]+2.*[[['t','t'],['photon']],[['t','t'],['photon']]]+2.*[[['t','t'],['Z']],[['t','t'],['Z']]])"
T5ttttZg.source 			  		= "CMS"
#+++++++ next mass plane block ++++++++++++++
for lsp in lsp_masses:
	p 						= T5ttttZg.addMassPlane(2*[[x, y, lsp]])
	p.figure    			= 'Fig.5c'
	p.figureUrl 			= 'http://cms-results.web.cern.ch/cms-results/public-results/publications/SUS-18-002/CMS-SUS-18-002_Figure_005-c.png'
	p.dataUrl   			= "http://cms-results.web.cern.ch/cms-results/public-results/publications/SUS-18-002/CMS-SUS-18-002_Figure_005-c.root"
	p.setSources(dataLabels	= ['expExclusionP1','expExclusionM1','obsExclusionP1','obsExclusionM1','expExclusion', 'obsExclusion', 'upperLimits','expectedUpperLimits'],
				dataFiles	= ['orig/CMS-SUS-18-002_Figure_005-c.root', 'orig/CMS-SUS-18-002_Figure_005-c.root', 'orig/CMS-SUS-18-002_Figure_005-c.root','orig/CMS-SUS-18-002_Figure_005-c.root','orig/CMS-SUS-18-002_Figure_005-c.root', 'orig/CMS-SUS-18-002_Figure_005-c.root', 'orig/CMS-SUS-18-002_Figure_005-c.root','orig/CMS-SUS-18-002_Figure_005-c.root'],
					units	= [ None, None,None, None,None, None, 'fb','fb' ],
				dataFormats	= ['root', 'root', 'root', 'root','root', 'root', 'root', 'root'],objectNames= ['exp1up;1','exp1dn;1','obs_XsecUp;1','obs_XsecDn;1','exp;1', 'obs;1','obs_XsecLimit;1','exp_XsecLimit;1'],
                    indices= [5,4,7,8,3,6,1,2])


#+++++++ next txName block ++++++++++++++
T6ttZg 					  		= dataset.addTxName('T6ttZg')
T6ttZg.checked 			  		= 'No'
T6ttZg.constraint 		  		= "[[['t'],['Z']],[['t'],['photon']]]+[[['t'],['photon']],[['t'],['photon']]]+[[['t'],['Z']],[['t'],['Z']]]"
T6ttZg.conditionDescription 		= None
T6ttZg.condition 			  		= "Csim([[['t'],['Z']],[['t'],['photon']]]+2.*[[['t'],['photon']],[['t'],['photon']]]+2.*[[['t'],['Z']],[['t'],['Z']]])"
T6ttZg.source 			  		= "CMS"
#+++++++ next mass plane block ++++++++++++++
for lsp in lsp_masses:
	p 						= T6ttZg.addMassPlane(2*[[x, y, lsp]])
	p.figure    			= 'Fig.5d'
	p.figureUrl 			= 'http://cms-results.web.cern.ch/cms-results/public-results/publications/SUS-18-002/CMS-SUS-18-002_Figure_005-d.png'
	p.dataUrl   			= "http://cms-results.web.cern.ch/cms-results/public-results/publications/SUS-18-002/CMS-SUS-18-002_Figure_005-d.root"
	p.setSources(dataLabels	= ['expExclusionP1','expExclusionM1','obsExclusionP1','obsExclusionM1','expExclusion', 'obsExclusion', 'upperLimits','expectedUpperLimits'],
				dataFiles	= ['orig/CMS-SUS-18-002_Figure_005-d.root', 'orig/CMS-SUS-18-002_Figure_005-d.root', 'orig/CMS-SUS-18-002_Figure_005-d.root','orig/CMS-SUS-18-002_Figure_005-d.root','orig/CMS-SUS-18-002_Figure_005-d.root', 'orig/CMS-SUS-18-002_Figure_005-d.root', 'orig/CMS-SUS-18-002_Figure_005-d.root','orig/CMS-SUS-18-002_Figure_005-d.root'],
					units	= [ None, None,None, None,None, None, 'fb', 'fb' ],
				dataFormats	= ['root', 'root', 'root', 'root','root', 'root', 'root', 'root'],objectNames= ['exp1up;1','exp1dn;1','obs_XsecUp;1','obs_XsecDn;1','exp;1', 'obs;1','obs_XsecLimit;1','exp_XsecLimit;1'],
                    indices= [5,4,7,8,3,6,1,2])
databaseCreator.create()
