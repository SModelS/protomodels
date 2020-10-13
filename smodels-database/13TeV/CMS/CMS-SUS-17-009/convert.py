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
info = MetaInfoInput('CMS-SUS-17-009')
info.url = 'https://cms-results.web.cern.ch/cms-results/public-results/publications/SUS-17-009/'
info.sqrts = 13
info.lumi = 35.9
info.prettyName = 'SFOS leptons + Etmiss'
info.private = False
info.arxiv = 'https://arxiv.org/abs/1806.05264'
info.contact ='cms-phys-conveners-sus@cern.ch'
info.publication = 'Phys. Lett.B 790 (2019) 140'
info.comment = ''

  

#+++++++ dataset block ++++++++++++++
dataset = DataSetInput('data')
dataset.setInfo(dataType = 'upperLimit', dataId = None)



#+++++++ next txName block ++++++++++++++
TSlepSlep = dataset.addTxName('TSlepSlep')
TSlepSlep.checked =""
TSlepSlep.constraint ="[[['e+']],[['e-']]]+[[['mu+']],[['mu-']]]"
TSlepSlep.conditionDescription ="[[['mu+']],[['mu-']]] > [[['e+']],[['e-']]]"
TSlepSlep.condition ="Cgtr([[['mu+']],[['mu-']]],[[['e+']],[['e-']]])"
TSlepSlep.source = "CMS"
#+++++++ next mass plane block ++++++++++++++
TSlepSlep_1 = TSlepSlep.addMassPlane(2*[[x, y]])
TSlepSlep_1.figure = "Fig. 5a"
TSlepSlep_1.figureUrl = "https://cms-results.web.cern.ch/cms-results/public-results/publications/SUS-17-009/CMS-SUS-17-009_Figure_005-a.png"
TSlepSlep_1.dataUrl = 'https://cms-results.web.cern.ch/cms-results/public-results/publications/SUS-17-009/CMS-SUS-17-009_Figure_005-a.root'
TSlepSlep_1.setSources(dataLabels= ['obsExclusion', 'expExclusion','expExclusionP1','expExclusionM1','obsExclusionP1','obsExclusionM1', 'upperLimits','expectedUpperLimits'],
                 dataFiles= ['orig/CMS-SUS-17-009_Figure_005-a.root', 'orig/CMS-SUS-17-009_Figure_005-a.root','orig/CMS-SUS-17-009_Figure_005-a.root', 'orig/CMS-SUS-17-009_Figure_005-a.root','orig/CMS-SUS-17-009_Figure_005-a.root', 'orig/CMS-SUS-17-009_Figure_005-a.root','orig/CMS-SUS-17-009_Figure_005-a.root', 'orig/CMS-SUS-17-009_Figure_005-a.root'],
                 dataFormats= ['root', 'root','root', 'root','root', 'root','root', 'root'],objectNames= ['gr_obs_smoothed;1','gr_exp_smoothed;1','gr_ep1s_smoothed;1','gr_em1s_smoothed;1','gr_op1s_smoothed;1','gr_om1s_smoothed;1','obs_xs0;1','exp_xs0;1'],units= [None,None,None,None,None,None,'fb','fb'],
                    indices= [1,2,3,4,7,8,10,16])


 
 
 
 #+++++++ next txName block ++++++++++++++
TSelSel = dataset.addTxName('TSelSel')
TSelSel.checked =""
TSelSel.constraint ="[[['e+']],[['e-']]]"
TSelSel.conditionDescription =None
TSelSel.condition =None
TSelSel.source = "CMS"
#+++++++ next mass plane block ++++++++++++++
TSelSel_1 = TSelSel.addMassPlane(2*[[x, y]])
TSelSel_1.constraint ="[[['e+']],[['e-']]]"
TSelSel_1.conditionDescription =None
TSelSel_1.condition =None
TSelSel_1.figure = "Fig. 6a"
TSelSel_1.figureUrl = "https://cms-results.web.cern.ch/cms-results/public-results/publications/SUS-17-009/CMS-SUS-17-009_Figure_006-a.png"
TSelSel_1.dataUrl = 'https://cms-results.web.cern.ch/cms-results/public-results/publications/SUS-17-009/CMS-SUS-17-009_Figure_006-a.root'
TSelSel_1.setSources(dataLabels= ['obsExclusion', 'expExclusion','expExclusionP1','expExclusionM1','obsExclusionP1','obsExclusionM1','upperLimits','expectedUpperLimits'],
                 dataFiles= ['orig/CMS-SUS-17-009_Figure_006-a.root', 'orig/CMS-SUS-17-009_Figure_006-a.root','orig/CMS-SUS-17-009_Figure_006-a.root', 'orig/CMS-SUS-17-009_Figure_006-a.root','orig/CMS-SUS-17-009_Figure_006-a.root', 'orig/CMS-SUS-17-009_Figure_006-a.root','orig/CMS-SUS-17-009_Figure_006-a.root', 'orig/CMS-SUS-17-009_Figure_006-a.root'],
                 dataFormats= ['root', 'root','root', 'root','root', 'root','root', 'root'],objectNames= ['gr_obs_smoothed;1','gr_exp_smoothed;1','gr_ep1s_smoothed;1','gr_em1s_smoothed;1','gr_op1s_smoothed;1','gr_om1s_smoothed;1','obs_xs0;1','exp_xs0;1'],units= [None,None,None,None,None,None,'fb','fb'],
                    indices= [1,2,3,4,7,8,10,16])




#+++++++ next txName block ++++++++++++++
TSmuSmu = dataset.addTxName('TSmuSmu')
TSmuSmu.checked =""
TSmuSmu.constraint ="[[['mu+']],[['mu-']]]"
TSmuSmu.conditionDescription =None
TSmuSmu.condition =None
TSmuSmu.source = "CMS" 
 #+++++++ next mass plane block ++++++++++++++
TSmuSmu_1 = TSmuSmu.addMassPlane(2*[[x, y]])
TSmuSmu_1.constraint ="[[['mu+']],[['mu-']]]"
TSmuSmu_1.conditionDescription =None
TSmuSmu_1.condition =None
TSmuSmu_1.figure = "Fig. 7a"
TSmuSmu_1.figureUrl = "https://cms-results.web.cern.ch/cms-results/public-results/publications/SUS-17-009/CMS-SUS-17-009_Figure_007-a.png"
TSmuSmu_1.dataUrl = 'https://cms-results.web.cern.ch/cms-results/public-results/publications/SUS-17-009/CMS-SUS-17-009_Figure_007-a.root'
TSmuSmu_1.setSources(dataLabels= ['obsExclusion', 'expExclusion','expExclusionP1','expExclusionM1','obsExclusionP1','obsExclusionM1','upperLimits','expectedUpperLimits'],
                 dataFiles= ['orig/CMS-SUS-17-009_Figure_007-a.root', 'orig/CMS-SUS-17-009_Figure_007-a.root','orig/CMS-SUS-17-009_Figure_007-a.root', 'orig/CMS-SUS-17-009_Figure_007-a.root','orig/CMS-SUS-17-009_Figure_007-a.root', 'orig/CMS-SUS-17-009_Figure_007-a.root','orig/CMS-SUS-17-009_Figure_007-a.root', 'orig/CMS-SUS-17-009_Figure_007-a.root'],
                 dataFormats= ['root', 'root','root', 'root','root', 'root','root', 'root'],objectNames= ['gr_obs_smoothed;1','gr_exp_smoothed;1','gr_ep1s_smoothed;1','gr_em1s_smoothed;1','gr_op1s_smoothed;1','gr_om1s_smoothed;1','obs_xs0;1','exp_xs0;1'],units= [None,None,None,None,None,None,'fb','fb'],
                    indices= [1,2,3,4,7,8,10,16]) 
 
 

 
databaseCreator.create()
