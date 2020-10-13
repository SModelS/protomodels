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
info = MetaInfoInput('CMS-SUS-17-005')
info.url = 'https://cms-results.web.cern.ch/cms-results/public-results/publications/SUS-17-005/'
info.sqrts = 13
info.lumi = 35.9
info.prettyName = '1L + multijets + Etmiss, top tagging'
info.private = False
info.arxiv = 'https://arxiv.org/abs/1805.05784'
info.contact ='cms-phys-conveners-sus@cern.ch'
info.publication = 'JHEP 1809 (2018) 065'
info.comment = 'using 1l cut and count analysis'



#+++++++ dataset block ++++++++++++++
dataset = DataSetInput('data')
dataset.setInfo(dataType = 'upperLimit', dataId = None)

#+++++++ next txName block ++++++++++++++
T2bbffff = dataset.addTxName('T2bbffff')
T2bbffff.checked = ''
T2bbffff.constraint = "3.375*[[['b','l','nu']],[['b','jet','jet']]]"
#T2bbffff.constraint = "5.06*([[['b','l','nu']],[['b','q','q']]]+[[['b','l','nu']],[['b','l','nu']]])"
#T2bbffff.constraint = "1.69*([[['b','l','nu']],[['b','jet','jet']]]+[[['b','jet','jet']],[['b','jet','jet']]])"
T2bbffff.conditionDescription = None
T2bbffff.condition = None
T2bbffff.source = "CMS"


#+++++++ next mass plane block ++++++++++++++
T2bbffff_1 = T2bbffff.addMassPlane(2*[[x, x-y]])
T2bbffff_1.figure = 'Fig.7-a'
T2bbffff_1.figureUrl = 'http://cms-results.web.cern.ch/cms-results/public-results/publications/SUS-17-005/CMS-SUS-17-005_Figure_007-a.png'
T2bbffff_1.dataUrl = 'http://cms-results.web.cern.ch/cms-results/public-results/publications/SUS-17-005/CMS-SUS-17-005_Figure_007-a.root'
T2bbffff_1.histoDataUrl = 'http://cms-results.web.cern.ch/cms-results/public-results/publications/SUS-17-005/CMS-SUS-17-005_Figure_007-a.root'
T2bbffff_1.dataUrl = 'http://cms-results.web.cern.ch/cms-results/public-results/publications/SUS-17-005/CMS-SUS-17-005_Figure_007-a.root'
T2bbffff_1.exclusionDataUrl = 'http://cms-results.web.cern.ch/cms-results/public-results/publications/SUS-17-005/CMS-SUS-17-005_Figure_007-a.root'
T2bbffff_1.setSources(dataLabels= ['obsExclusionM1','expExclusionM1', 'expExclusionP1','obsExclusion', 'expExclusion','upperLimits','expectedUpperLimits','obsExclusionP1'],
                 dataFiles= ['orig/CMS-SUS-17-005_Figure_007-a.root', 'orig/CMS-SUS-17-005_Figure_007-a.root','orig/CMS-SUS-17-005_Figure_007-a.root',
                 'orig/CMS-SUS-17-005_Figure_007-a.root','orig/CMS-SUS-17-005_Figure_007-a.root', 'orig/CMS-SUS-17-005_Figure_007-a.root','orig/CMS-SUS-17-005_Figure_007-a.root',
                 'orig/CMS-SUS-17-005_Figure_007-a.root'],
                 dataFormats= ['root','root','root','root','root','root','root','root'],objectNames= ['contour_ObsMinus1','contour_ExpMinus1','contour_ExpPlus1','contour_Obs','contour_Exp','xsecUL_Obs','xsecUL_Exp','contour_ObsPlus1'])


"""
#+++++++ next mass plane block ++++++++++++++
T2bbffff_2 = T2bbffff.addMassPlane(2*[[x, x-y]])
T2bbffff_2.figure = 'Fig.9-a'
T2bbffff_2.figureUrl = 'http://cms-results.web.cern.ch/cms-results/public-results/publications/SUS-17-005/CMS-SUS-17-005_Figure_009-a.png'
T2bbffff_2.dataUrl = 'http://cms-results.web.cern.ch/cms-results/public-results/publications/SUS-17-005/CMS-SUS-17-005_Figure_009-a.root'
T2bbffff_2.histoDataUrl = 'http://cms-results.web.cern.ch/cms-results/public-results/publications/SUS-17-005/CMS-SUS-17-005_Figure_009-a.root'
T2bbffff_2.dataUrl = 'http://cms-results.web.cern.ch/cms-results/public-results/publications/SUS-17-005/CMS-SUS-17-005_Figure_009-a.root'
T2bbffff_2.exclusionDataUrl = 'http://cms-results.web.cern.ch/cms-results/public-results/publications/SUS-17-005/CMS-SUS-17-005_Figure_009-a.root'
T2bbffff_2.setSources(dataLabels= ['obsExclusionM1','expExclusionM1', 'expExclusionP1','obsExclusion', 'expExclusion','upperLimits','expectedUpperLimits','obsExclusionP1'],
                 dataFiles= ['orig/CMS-SUS-17-005_Figure_009-a.root', 'orig/CMS-SUS-17-005_Figure_009-a.root','orig/CMS-SUS-17-005_Figure_009-a.root',
                 'orig/CMS-SUS-17-005_Figure_009-a.root','orig/CMS-SUS-17-005_Figure_009-a.root',
                 'orig/CMS-SUS-17-005_Figure_009-a.root','orig/CMS-SUS-17-005_Figure_009-a.root','orig/CMS-SUS-17-005_Figure_009-a.root'],
                 dataFormats= ['root','root','root','root','root','root','root','root'],objectNames= ['contour_ObsMinus1','contour_ExpMinus1','contour_ExpPlus1','contour_Obs','contour_Exp','xsecUL_Obs','xsecUL_Exp','contour_ObsPlus1'])
"""

#+++++++ next txName block ++++++++++++++
T6bbWWoff = dataset.addTxName('T6bbWWoff')
T6bbWWoff.checked = ''
T6bbWWoff.constraint = "3.375*[[['b'],['l','nu']],[['b'],['jet','jet']]]"
#T6bbWWoff.constraint = "5.06*([[['b'],['l','nu']],[['b'],['q','q']]]+[[['b'],['l','nu']],[['b'],['l','nu']]])"
#T6bbWWoff.constraint = "1.69*([[['b'],['l','nu']],[['b'],['jet','jet']]]+[[['b'],['jet','jet']],[['b'],['jet','jet']]])"
T6bbWWoff.conditionDescription = None
T6bbWWoff.condition = None
T6bbWWoff.source = "CMS"


#+++++++ next mass plane block ++++++++++++++
T6bbWWoff_1 = T6bbWWoff.addMassPlane(2*[[x,x - 0.5*y, x-y]])
T6bbWWoff_1.figure = 'Fig.8'
T6bbWWoff_1.figureUrl = 'http://cms-results.web.cern.ch/cms-results/public-results/publications/SUS-17-005/CMS-SUS-17-005_Figure_008.png'
T6bbWWoff_1.dataUrl = 'http://cms-results.web.cern.ch/cms-results/public-results/publications/SUS-17-005/CMS-SUS-17-005_Figure_008.root'
T6bbWWoff_1.histoDataUrl = 'http://cms-results.web.cern.ch/cms-results/public-results/publications/SUS-17-005/CMS-SUS-17-005_Figure_008.root'
T6bbWWoff_1.dataUrl = 'http://cms-results.web.cern.ch/cms-results/public-results/publications/SUS-17-005/CMS-SUS-17-005_Figure_008.root'
T6bbWWoff_1.exclusionDataUrl = 'http://cms-results.web.cern.ch/cms-results/public-results/publications/SUS-17-005/CMS-SUS-17-005_Figure_008.root'
T6bbWWoff_1.setSources(dataLabels= ['obsExclusionM1','expExclusionM1', 'expExclusionP1','obsExclusion', 'expExclusion','upperLimits','expectedUpperLimits','obsExclusionP1'],
                 dataFiles= ['orig/CMS-SUS-17-005_Figure_008.root', 'orig/CMS-SUS-17-005_Figure_008.root','orig/CMS-SUS-17-005_Figure_008.root',
                 'orig/CMS-SUS-17-005_Figure_008.root','orig/CMS-SUS-17-005_Figure_008.root', 'orig/CMS-SUS-17-005_Figure_008.root','orig/CMS-SUS-17-005_Figure_008.root',
                 'orig/CMS-SUS-17-005_Figure_008.root'],
                 dataFormats= ['root','root','root','root','root','root','root','root'],objectNames= ['contour_ObsMinus1','contour_ExpMinus1','contour_ExpPlus1','contour_Obs','contour_Exp','xsecUL_Obs','xsecUL_Exp','contour_ObsPlus1'])

"""
#+++++++ next mass plane block ++++++++++++++
T6bbWWoff_2 = T6bbWWoff.addMassPlane(2*[[x,x-0.5*y, x-y]])
T6bbWWoff_2.figure = 'Fig.9-b'
T6bbWWoff_2.figureUrl = 'http://cms-results.web.cern.ch/cms-results/public-results/publications/SUS-17-005/CMS-SUS-17-005_Figure_009-b.png'
T6bbWWoff_2.dataUrl = 'http://cms-results.web.cern.ch/cms-results/public-results/publications/SUS-17-005/CMS-SUS-17-005_Figure_009-b.root'
T6bbWWoff_2.histoDataUrl = 'http://cms-results.web.cern.ch/cms-results/public-results/publications/SUS-17-005/CMS-SUS-17-005_Figure_009-b.root'
T6bbWWoff_2.dataUrl = 'http://cms-results.web.cern.ch/cms-results/public-results/publications/SUS-17-005/CMS-SUS-17-005_Figure_009-b.root'
T6bbWWoff_2.exclusionDataUrl = 'http://cms-results.web.cern.ch/cms-results/public-results/publications/SUS-17-005/CMS-SUS-17-005_Figure_009-b.root'
T6bbWWoff_2.setSources(dataLabels= ['obsExclusionM1','expExclusionM1', 'expExclusionP1','obsExclusion', 'expExclusion','upperLimits','expectedUpperLimits','obsExclusionP1'],
                 dataFiles= ['orig/CMS-SUS-17-005_Figure_009-b.root', 'orig/CMS-SUS-17-005_Figure_009-b.root','orig/CMS-SUS-17-005_Figure_009-b.root',
                 'orig/CMS-SUS-17-005_Figure_009-b.root','orig/CMS-SUS-17-005_Figure_009-b.root', 'orig/CMS-SUS-17-005_Figure_009-b.root','orig/CMS-SUS-17-005_Figure_009-b.root',
                 'orig/CMS-SUS-17-005_Figure_009-b.root'],
                 dataFormats= ['root','root','root','root','root','root','root','root'],objectNames= ['contour_ObsMinus1','contour_ExpMinus1','contour_ExpPlus1','contour_Obs','contour_Exp','xsecUL_Obs','xsecUL_Exp','contour_ObsPlus1'])
"""
databaseCreator.create()
