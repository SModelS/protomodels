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
info = MetaInfoInput('CMS-SUS-17-010')
info.url = 'http://cms-results.web.cern.ch/cms-results/public-results/publications/SUS-17-010'
info.sqrts = 13
info.lumi = 35.9 
info.prettyName = '2L stop'
info.private = False
info.arxiv       = 'https://arxiv.org/abs/1807.07799'
info.contact     = 'cms-phys-conveners-sus@cern.ch'
info.publication = 'JHEP 1811 (2018) 079'
info.comment     = ''


#+++++++ dataset block ++++++++++++++
dataset = DataSetInput('data')
dataset.setInfo(dataType = 'upperLimit', dataId = None)

#+++++++ next txName block ++++++++++++++
T2tt = dataset.addTxName('T2tt')
T2tt.checked =''
T2tt.constraint = "[[['t']],[['t']]]"
T2tt.conditionDescription = None
T2tt.condition = None
T2tt.source = "CMS"
T2tt.massConstraint = None
T2ttoff = dataset.addTxName('T2ttoff')
T2ttoff.checked =''
T2ttoff.constraint ="[[['b','W']],[['b','W']]]"
T2ttoff.conditionDescription = None
T2ttoff.condition = None
T2ttoff.massConstraint=[['80<dm<169'],['80<dm<169']]
T2ttoff.source = "CMS"
#+++++++ next mass plane block ++++++++++++++
T2tt_1 = T2tt.addMassPlane(2*[[x, y]])
T2tt_1.figure = 'Figure 9a'
T2tt_1.figureUrl = 'http://cms-results.web.cern.ch/cms-results/public-results/publications/SUS-17-010/CMS-SUS-17-010_Figure_009-a.png'
T2tt_1.dataUrl = 'http://cms-results.web.cern.ch/cms-results/public-results/publications/SUS-17-010/'
T2tt_1.histoDataUrl = 'http://cms-results.web.cern.ch/cms-results/public-results/publications/SUS-17-010/CMS-SUS-17-010_Figure_009-a.root'
T2tt_1.exclusionDataUrl = 'http://cms-results.web.cern.ch/cms-results/public-results/publications/SUS-17-010/CMS-SUS-17-010_Figure_009-a.root'
T2tt_1.setSources(dataLabels= ['expExclusion','expExclusionM1','expExclusionP1','obsExclusionP1','obsExclusion','obsExclusionM1','expectedUpperLimits','upperLimits'],
                 dataFiles= ['orig/CMS-SUS-17-010_Figure_009-a.root', 'orig/CMS-SUS-17-010_Figure_009-a.root', 'orig/CMS-SUS-17-010_Figure_009-a.root', 'orig/CMS-SUS-17-010_Figure_009-a.root','orig/CMS-SUS-17-010_Figure_009-a.root', 'orig/CMS-SUS-17-010_Figure_009-a.root', 'orig/CMS-SUS-17-010_Figure_009-a.root', 'orig/CMS-SUS-17-010_Figure_009-a.root'],
                 dataFormats= ['root', 'root', 'root', 'root','root', 'root', 'root', 'root'],objectNames= ['gr_Exp', 'gr_ExpM', 'gr_ExpP', 'gr_ObsP','gr_Obs', 'gr_ObsM', 'hXsec_exp', 'hXsec_obs'],
                    indices= [1, 2, 3, 4, 5, 6, 7,8],units=[None,None,None,None,None,None,'pb','pb'])
T2ttoff.addMassPlane(T2tt_1)







#+++++++ next txName block ++++++++++++++
T6bbWW = dataset.addTxName('T6bbWW')
T6bbWW.checked = ''
T6bbWW.constraint ="[[['b'],['W']],[['b'],['W']]]"
T6bbWW.conditionDescription = None
T6bbWW.condition = None
T6bbWW.source = "CMS"
T6bbWW.massConstraint = None
T6bbWWoff = dataset.addTxName('T6bbWWoff')
T6bbWWoff.checked =''
T6bbWWoff.constraint ="2.3*([[['b'],['jet','jet']],[['b'],['jet','jet']]])"
T6bbWWoff.conditionDescription = None
T6bbWWoff.condition = None
T6bbWWoff.massConstraint = [['dm >= 0.0', 'dm <= 76.0'], ['dm >= 0.0', 'dm <= 76.0']]
T6bbWWoff.source = "CMS"
#+++++++ next mass plane block ++++++++++++++
T6bbWW_1 = T6bbWW.addMassPlane(2*[[x, (0.5*x+0.5*y), y]])
T6bbWW_1.figure = 'Figure 9b'
T6bbWW_1.figureUrl = 'http://cms-results.web.cern.ch/cms-results/public-results/publications/SUS-17-010/CMS-SUS-17-010_Figure_009-b.png'
T6bbWW_1.dataUrl ='http://cms-results.web.cern.ch/cms-results/public-results/publications/SUS-17-010/'
T6bbWW_1.histoDataUrl ='http://cms-results.web.cern.ch/cms-results/public-results/publications/SUS-17-010/CMS-SUS-17-010_Figure_009-b.root'
T6bbWW_1.exclusionDataUrl ='http://cms-results.web.cern.ch/cms-results/public-results/publications/SUS-17-010/CMS-SUS-17-010_Figure_009-b.root'
T6bbWW_1.setSources(dataLabels= ['expExclusion','expExclusionM1','expExclusionP1','obsExclusionP1','obsExclusion','obsExclusionM1','expectedUpperLimits','upperLimits'],
                 dataFiles= ['orig/CMS-SUS-17-010_Figure_009-b.root', 'orig/CMS-SUS-17-010_Figure_009-b.root','orig/CMS-SUS-17-010_Figure_009-b.root', 'orig/CMS-SUS-17-010_Figure_009-b.root','orig/CMS-SUS-17-010_Figure_009-b.root', 'orig/CMS-SUS-17-010_Figure_009-b.root','orig/CMS-SUS-17-010_Figure_009-b.root', 'orig/CMS-SUS-17-010_Figure_009-b.root'],
                 dataFormats= ['root', 'root','root', 'root','root', 'root','root', 'root'],objectNames= ['gr_Exp','gr_ExpM','gr_ExpP','gr_ObsP','gr_Obs','gr_ObsM','hXsec_exp','hXsec_obs'])
T6bbWWoff.addMassPlane(T6bbWW_1)







#+++++++ next txName block ++++++++++++++
TChipChimSlepSnu = dataset.addTxName('TChipChimSlepSnu')
TChipChimSlepSnu.checked =''
TChipChimSlepSnu.constraint = "[[['L'],[nu]],[[nu],['L']]]+[[['L'],[nu]],[[nu],['L']]]+[[['L'],[nu]],[['L'],[nu]]]+[[[nu],['L']],[[nu],['L']]]"

TChipChimSlepSnu.constraint ="[[['L-'],['nu']],[['nu'],['L+']]] + [[['L+'],['nu']],[['nu'],['L-']]] + [[['L+'],['nu']],[['L-'],['nu']]] + [[['nu'],['L+']],[['nu'],['L-']]]"
TChipChimSlepSnu.conditionDescription ="[[['L-'],['nu']],[['nu'],['L+']]] ~ [[['L+'],['nu']],[['nu'],['L-']]], [[['L-'],['nu']],[['nu'],['L+']]] ~ [[['L+'],['nu']],[['L-'],['nu']]], [[['L-'],['nu']],[['nu'],['L+']]] ~ [[['nu'],['L+']],[['nu'],['L-']]],[[['L-'],['nu']],[['nu'],['L+']]] > 2.7*[[['ta-'],['nu']],['nu'],['L+']]],[[['L-'],['nu']],[['nu'],['L+']]] > 2.7*[[['L-'],['nu']],[['nu'],['ta+']]], [[['L+'],['nu']],[['nu'],['L-']]] > 2.7*[[['ta+'],['nu']],[['nu'],['L-']]], [[['L+'],['nu']],[['nu'],['L-']]] > 2.7*[[['L+'],['nu']],[['nu'],['ta-']]], [[['L+'],['nu']],[['L-'],['nu']]] > 2.7*[[['ta+'],['nu']],[['L-'],['nu']]], [[['L+'],['nu']],[['L-'],['nu']]] > 2.7*[[['L+'],['nu']],[['ta-'],['nu']]], [[['nu'],['L+']],[['nu'],['L-']]] > 2.7*[[['nu'],['ta+']],[['nu'],[L-']]], [[['nu'],['L+']],[['nu'],['L-']]] > 2.7*[[['nu'],['L+']],[['nu'],[ta-']]], [[['L-'],['nu']],[['nu'],['L+']]] > 2.7*[[['e-'],['nu']],[['nu'],['L+']]], [[['L-'],['nu']],[['nu'],['L+']]] > 2.7*[[['L-'],['nu']],[['nu'],['e+']]], [[['L+'],['nu']],[['nu'],['L-']]] > 2.7*[[['e+'],['nu']],[['nu'],['L-']]],[[['L+'],['nu']],[['nu'],['L-']]] > 2.7*[[['L+'],['nu']],[['nu'],['e-']]],[[['L+'],['nu']],[['L-'],['nu']]] > 2.7*[[['e+'],['nu']],[['L-'],['nu']]],[[['L+'],['nu']],[['L-'],['nu']]] > 2.7*[[['L+'],['nu']],[['e-'],['nu']]],[[['nu'],['L+']],[['nu'],['L-']]] > 2.7*[[['nu'],['e+']],[['nu'],['L-']]], [[['nu'],['L+']],[['nu'],['L-']]] > 2.7*[[['nu'],['L+']],[['nu'],['e-']]]"
TChipChimSlepSnu.condition ="Csim([[['L-'],['nu']],[['nu'],['L+']]],[[['L+'],['nu']],[['nu'],['L-']]],[[['L+'],['nu']],[['L-'],['nu']]],[[['nu'],['L+']],[['nu'],['L-']]]); Cgtr([[['L-'],['nu']],[['nu'],['L+']]],3.*[[['ta-'],['nu']],[['nu'],['L+']]]); Cgtr([[['L-'],['nu']],[['nu'],['L+']]],3.*[[['L-'],['nu']],[['nu'],['ta+']]]); Cgtr([[['L+'],['nu']],[['nu'],['L-']]],3.*[[['ta+'],['nu']],[['nu'],['L-']]]); Cgtr([[['L+'],['nu']],[['nu'],['L-']]],3.* [[['L+'],['nu']],[['nu'],['ta-']]]); Cgtr([[['L+'],['nu']],[['L-'],['nu']]],3.*[[['ta+'],['nu']],[['L-'],['nu']]]); Cgtr([[['L+'],['nu']],[['L-'],['nu']]],3.*[[['L+'],['nu']],[['ta-'],['nu']]]); Cgtr([[['nu'],['L+']],[['nu'],[L-']]],3.*[[['nu'],['ta+']],[['nu'],[L-']]]); Cgtr([[['nu'],['L+']],[['nu'],[L-']]],3.*[[['nu'],['L+']],[['nu'],[ta-']]]); Cgtr([[['L-'],['nu']],[['nu'],['L+']]],3.*[[['e-'],['nu']],[['nu'],['L+']]]); Cgtr([[['L-'],['nu']],[['nu'],['L+']]],3.*[[['L-'],['nu']],[['nu'],['e+']]]); Cgtr([[['L+'],['nu']],[['nu'],['L-']]],3.*[[['e+'],['nu']],[['nu'],['L-']]]); Cgtr([[['L+'],['nu']],[['nu'],['L-']]],3.* [[['L+'],['nu']],[['nu'],['e-']]]); Cgtr([[['L+'],['nu']],[['L-'],['nu']]],3.*[[['e+'],['nu']],[['L-'],['nu']]]); Cgtr([[['L+'],['nu']],[['L-'],['nu']]],3.*[[['L+'],['nu']],[['e-'],['nu']]]); Cgtr([[['nu'],['L+']],[['nu'],[L-']]],3.*[[['nu'],['e+']],[['nu'],[L-']]]); Cgtr([[['nu'],['L+']],[['nu'],[L-']]],3.*[[['nu'],['L+']],[['nu'],[e-']]])"

TChipChimSlepSnu.source = "CMS"
#+++++++ next mass plane block ++++++++++++++
TChipChimSlepSnu_1 = TChipChimSlepSnu.addMassPlane(2*[[x, 0.5*x+0.5*y, y]])
TChipChimSlepSnu_1.figure = 'Fig. 8a'
TChipChimSlepSnu_1.figureUrl ='http://cms-results.web.cern.ch/cms-results/public-results/publications/SUS-17-010/CMS-SUS-17-010_Figure_008-a.png'
TChipChimSlepSnu_1.dataUrl = 'http://cms-results.web.cern.ch/cms-results/public-results/publications/SUS-17-010/'
TChipChimSlepSnu_1.histoDataUrl ='http://cms-results.web.cern.ch/cms-results/public-results/publications/SUS-17-010/CMS-SUS-17-010_Figure_008-a.root'

TChipChimSlepSnu_1.exclusionDataUrl ='http://cms-results.web.cern.ch/cms-results/public-results/publications/SUS-17-010/CMS-SUS-17-010_Figure_008-a.root'
TChipChimSlepSnu_1.setSources(dataLabels= ['expExclusion','expExclusionM1','expExclusionP1','obsExclusionP1','obsExclusion','obsExclusionM1','expectedUpperLimits','upperLimits'],
                 dataFiles= ['orig/CMS-SUS-17-010_Figure_008-a.root', 'orig/CMS-SUS-17-010_Figure_008-a.root', 'orig/CMS-SUS-17-010_Figure_008-a.root', 'orig/CMS-SUS-17-010_Figure_008-a.root', 'orig/CMS-SUS-17-010_Figure_008-a.root', 'orig/CMS-SUS-17-010_Figure_008-a.root', 'orig/CMS-SUS-17-010_Figure_008-a.root', 'orig/CMS-SUS-17-010_Figure_008-a.root'],
                 dataFormats= ['root', 'root','root', 'root','root', 'root','root', 'root'],objectNames= ['gr_Exp','gr_ExpM','gr_ExpP','gr_ObsP','gr_Obs','gr_ObsM','hXsec_exp','hXsec_obs'],indices= [1, 2, 3, 4, 5, 6, 7, 8],units= [None,None,None,None,None,None,'pb','pb'])



databaseCreator.create()
