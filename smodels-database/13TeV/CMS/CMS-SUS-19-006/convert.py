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
info             = MetaInfoInput('CMS-SUS-19-006')
info.url         = 'http://cms-results.web.cern.ch/cms-results/public-results/publications/SUS-19-006/index.html'
info.sqrts       = 13
info.lumi        = 137
info.prettyName  = '0L + jets, MHT'
info.private     = False
info.arxiv       = 'https://arxiv.org/abs/1908.04722'
info.contact     = 'cms-phys-conveners-sus@cern.ch'
info.publication = 'JHEP 10 (2019) 244'
info.comment     = 'Observed and expected limits; dense grid.'


#+++++++ dataset block ++++++++++++++
dataset = DataSetInput('data')
dataset.setInfo(dataType = 'upperLimit', dataId = None)


#====================  T1tttt  ==============================  

#+++++++txName block++++++++++++++++++++

T1tttt=dataset.addTxName('T1tttt')
T1tttt.checked=''
T1tttt.constraint="[[['t','t']],[['t','t']]]"
T1tttt.condition=None
T1tttt.conditionDescription = None
T1tttt.source="CMS"
T1tttt.massConstraint=[['dm>=338.0'],['dm>=338.0']]

T1ttttoff=dataset.addTxName('T1ttttoff')
T1ttttoff.checked=''
T1ttttoff.constraint="[[['b','W','b','W']],[['b','W','b','W']]]"
T1ttttoff.condition=None
T1ttttoff.source="CMS"
T1ttttoff.massConstraint=[['160 < dm < 338.0'], ['160 < dm < 338.0']]

#++++++mass plane block+++++++++

T1tttt_1 = T1tttt.addMassPlane(2*[[x,y]])
T1tttt_1.figure='Fig. 13a'
T1tttt_1.figureUrl='http://cms-results.web.cern.ch/cms-results/public-results/publications/SUS-19-006/CMS-SUS-19-006_Figure_013-a.png'
T1tttt_1.dataUrl='https://doi.org/10.17182/hepdata.90835.v1/t18'
T1tttt_1.exclusionDataUrl='https://doi.org/10.17182/hepdata.90835.v1/t19'
T1tttt_1.setSources(dataLabels=['expExclusion','expExclusionM1','expExclusionP1','obsExclusion','obsExclusionM1','obsExclusionP1','upperLimits','expectedUpperLimits'],
                    dataFiles=['orig/T1ttttexpectedmasslimitcurve.csv','orig/T1ttttexpected-1s.d.masslimitcurve.csv','orig/T1ttttexpected+1s.d.masslimitcurve.csv','orig/T1ttttobservedmasslimitcurve.csv','orig/T1ttttobserved-1s.d.(theory)masslimitcurve.csv','orig/T1ttttobserved+1s.d.(theory)masslimitcurve.csv','orig/T1tttt_obs_upperlimits.csv','orig/T1tttt_exp_upperlimits.csv'],
                    dataFormats=['csv','csv','csv','csv','csv','csv','csv','csv'],units=[None,None,None,None,None,None,'pb','pb'])
T1ttttoff.addMassPlane(T1tttt_1)


#====================  T1bbbb  ==============================  

#+++++txName block +++++++++++++++++
T1bbbb=dataset.addTxName('T1bbbb')
T1bbbb.checked=''
T1bbbb.constraint="[[['b','b']],[['b','b']]]"
T1bbbb.condition=None
T1bbbb.conditionDescription = None
T1bbbb.source="CMS"

#++++++mass plane block+++++++++
T1bbbb_1 = T1bbbb.addMassPlane(2*[[x,y]])
T1bbbb_1.figure='Fig. 13b'
T1bbbb_1.figureUrl='http://cms-results.web.cern.ch/cms-results/public-results/publications/SUS-19-006/CMS-SUS-19-006_Figure_013-b.png'
T1bbbb_1.dataUrl='https://doi.org/10.17182/hepdata.90835.v1/t25'
T1bbbb_1.exclusionDataUrl='https://doi.org/10.17182/hepdata.90835.v1/t26'
T1bbbb_1.setSources(dataLabels=['expExclusion','expExclusionM1','expExclusionP1','obsExclusion','obsExclusionM1','obsExclusionP1','upperLimits','expectedUpperLimits'],
                    dataFiles=['orig/T1bbbbexpectedmasslimitcurve.csv','orig/T1bbbbexpected-1s.d.masslimitcurve.csv','orig/T1bbbbexpected+1s.d.masslimitcurve.csv','orig/T1bbbbobservedmasslimitcurve.csv','orig/T1bbbbobserved-1s.d.(theory)masslimitcurve.csv','orig/T1bbbbobserved+1s.d.(theory)masslimitcurve.csv','orig/T1bbbb_obs_upperlimits.csv','orig/T1bbbb_exp_upperlimits.csv'],
                    dataFormats=['csv','csv','csv','csv','csv','csv','csv','csv'],units=[None,None,None,None,None,None,'pb','pb'])


#====================   T1   ==============================  

#+++++txName block +++++++++++++++++
T1=dataset.addTxName('T1')
T1.checked=''
T1.constraint="[[['jet','jet']],[['jet','jet']]]"
T1.condition=None
T1.conditionDescription = None
T1.source="CMS"

#++++++next mass plane block+++++++++

T1_1 = T1.addMassPlane(2*[[x,y]])
T1_1.figure='Fig. 13c'
T1_1.figureUrl='http://cms-results.web.cern.ch/cms-results/public-results/publications/SUS-19-006/CMS-SUS-19-006_Figure_013-c.png'
T1_1.dataUrl='https://doi.org/10.17182/hepdata.90835.v1/t32'
T1_1.exclusionDataUrl='https://doi.org/10.17182/hepdata.90835.v1/t33'
T1_1.setSources(dataLabels=['expExclusion','expExclusionM1','expExclusionP1','obsExclusion','obsExclusionM1','obsExclusionP1','upperLimits','expectedUpperLimits'],
                    dataFiles=['orig/T1qqqqexpectedmasslimitcurve.csv','orig/T1qqqqexpected-1s.d.masslimitcurve.csv','orig/T1qqqqexpected+1s.d.masslimitcurve.csv','orig/T1qqqqobservedmasslimitcurve.csv','orig/T1qqqqobserved-1s.d.(theory)masslimitcurve.csv','orig/T1qqqqobserved+1s.d.(theory)masslimitcurve.csv','orig/T1qqqq_obs_upperlimits.csv','orig/T1qqqq_exp_upperlimits.csv'],
                    dataFormats=['csv','csv','csv','csv','csv','csv','csv','csv'],units=[None,None,None,None,None,None,'pb','pb'])


#====================   T2tt   ==============================  

#+++++txName block +++++++++++++++++

T2tt=dataset.addTxName('T2tt')
T2tt.checked=''
T2tt.constraint="[[['t']],[['t']]]"
T2tt.condition=None
T2tt.conditionDescription = None
T2tt.massConstraint=[['dm>169'],['dm>169']]
T2tt.source="CMS"

T2ttoff=dataset.addTxName('T2ttoff')
T2ttoff.checked=''
T2ttoff.constraint="[[['b','W']],[['b','W']]]"
T2ttoff.condition=None
T2ttoff.conditionDescription=None
T2ttoff.massConstraint=[['80<dm<169'],['80<dm<169']]
T2ttoff.source="CMS"

#++++++next mass plane block+++++++++

T2tt_1 = T2tt.addMassPlane(2*[[x,y]])
T2tt_1.figure='Fig. 14-a'
T2tt_1.figureUrl='http://cms-results.web.cern.ch/cms-results/public-results/publications/SUS-19-006/CMS-SUS-19-006_Figure_014-a.png'
T2tt_1.dataUrl='https://doi.org/10.17182/hepdata.90835.v1/t46'
T2tt_1.exclusionDataUrl='https://doi.org/10.17182/hepdata.90835.v1/t47'
T2tt_1.setSources(dataLabels=['expExclusion','expExclusionM1','expExclusionP1','obsExclusion','obsExclusionM1','obsExclusionP1','upperLimits','expectedUpperLimits'],
                    dataFiles=['orig/T2ttexpectedmasslimitcurve.csv','orig/T2ttexpected-1s.d.masslimitcurve.csv','orig/T2ttexpected+1s.d.masslimitcurve.csv','orig/T2ttobservedmasslimitcurve.csv','orig/T2ttobserved-1s.d.(theory)masslimitcurve.csv','orig/T2ttobserved+1s.d.(theory)masslimitcurve.csv','orig/T2tt_obs_upperlimits.csv','orig/T2tt_exp_upperlimits.csv'],
                    dataFormats=['csv','csv','csv','csv','csv','csv','csv','csv'],units=[None,None,None,None,None,None,'pb','pb'])
T2ttoff.addMassPlane(T2tt_1)


#====================   T2bb   ==============================  

#+++++txName block +++++++++++++++++

T2bb=dataset.addTxName('T2bb')
T2bb.checked=''
T2bb.constraint="[[['b']],[['b']]]"
T2bb.condition=None
T2bb.conditionDescription = None
T2bb.source="CMS"

#++++++next mass plane block+++++++++

T2bb_1 = T2bb.addMassPlane(2*[[x,y]])
T2bb_1.figure='Fig. 14b'
T2bb_1.figureUrl='http://cms-results.web.cern.ch/cms-results/public-results/publications/SUS-19-006/CMS-SUS-19-006_Figure_014-b.png'
T2bb_1.dataUrl='https://doi.org/10.17182/hepdata.90835.v1/t53'
T2bb_1.exclusionDataUrl='https://doi.org/10.17182/hepdata.90835.v1/t54'
T2bb_1.setSources(dataLabels=['expExclusion','expExclusionM1','expExclusionP1','obsExclusion','obsExclusionM1','obsExclusionP1','upperLimits','expectedUpperLimits'],
                    dataFiles=['orig/T2bbexpectedmasslimitcurve.csv','orig/T2bbexpected-1s.d.masslimitcurve.csv','orig/T2bbexpected+1s.d.masslimitcurve.csv','orig/T2bbobservedmasslimitcurve.csv','orig/T2bbobserved-1s.d.(theory)masslimitcurve.csv','orig/T2bbobserved+1s.d.(theory)masslimitcurve.csv','orig/T2bb_obs_upperlimits.csv','orig/T2bb_exp_upperlimits.csv'],
                    dataFormats=['csv','csv','csv','csv','csv','csv','csv','csv'],units=[None,None,None,None,None,None,'pb','pb'])


#====================   T2   ==============================  

#+++++txName block +++++++++++++++++

T2=dataset.addTxName('T2')
T2.checked=''
T2.constraint="[[['jet']],[['jet']]]"
T2.condition=None
T2.conditionDescription = None
T2.source="CMS"

#++++++next mass plane block+++++++++

T2_1 = T2.addMassPlane(2*[[x,y]])
T2_1.figure='Fig. 14c'
T2_1.figureUrl='http://cms-results.web.cern.ch/cms-results/public-results/publications/SUS-19-006/CMS-SUS-19-006_Figure_014-c.png'
T2_1.dataUrl='https://doi.org/10.17182/hepdata.90835.v1/t60'
T2_1.exclusionDataUrl='https://doi.org/10.17182/hepdata.90835.v1/t61'
T2_1.setSources(dataLabels=['expExclusion','expExclusionM1','expExclusionP1','obsExclusion','obsExclusionM1','obsExclusionP1','upperLimits','expectedUpperLimits'],
                    dataFiles=['orig/T2qqexpectedmasslimitcurve.csv','orig/T2qqexpected-1s.d.masslimitcurve.csv','orig/T2qqexpected+1s.d.masslimitcurve.csv','orig/T2qqobservedmasslimitcurve.csv','orig/T2qqobserved-1s.d.(theory)masslimitcurve.csv','orig/T2qqobserved+1s.d.(theory)masslimitcurve.csv','orig/T2qq_obs_upperlimits.csv','orig/T2qq_exp_upperlimits.csv'],
                    dataFormats=['csv','csv','csv','csv','csv','csv','csv','csv'],units=[None,None,None,None,None,None,'pb','pb'])

#==================== end ==============================  

databaseCreator.create()
