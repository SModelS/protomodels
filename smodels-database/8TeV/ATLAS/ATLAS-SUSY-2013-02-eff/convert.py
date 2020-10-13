#!/usr/bin/env python3

"""
.. module:: convert.py
   :synopsis: used to create info.txt and the <txname>.txt files.

"""
import sys
import os
import argparse

argparser = argparse.ArgumentParser(description =
    'create info.txt, txname.txt, twiki.txt and sms.py')
argparser.add_argument ('-utilsPath', '--utilsPath',
    help = 'path to the package smodels_utils',\
    type = str)
argparser.add_argument ('-smodelsPath', '--smodelsPath',
    help = 'path to the package smodels_utils',\
    type = str)
argparser.add_argument ('-no', '--noUpdate',
    help = 'do not update the lastUpdate field.',\
    action= "store_true" )
argparser.add_argument ('-t', '--ntoys',
    help = 'number of toys to throw [100000]',\
    type = int, default=200000  )
args = argparser.parse_args()

if args.noUpdate:
    os.environ["SMODELS_NOUPDATE"]="1"

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

DataSetInput.ntoys = args.ntoys

#+++++++ global info block ++++++++++++++
info = MetaInfoInput('ATLAS-SUSY-2013-02')
info.comment = 'The recast maps	have only 10 out of 15 SRs, ATLAS official have 15. So for the recast T2,T5,TGQ we have 5 less SR wrt the official T1.We recast the T2 to cover the compressed squark-LSP region'
info.sqrts = '8.0'
info.private = False
info.lumi = '20.3'
info.publication = 'http://link.springer.com/article/10.1007/JHEP09%282014%29176'
info.url = 'https://atlas.web.cern.ch/Atlas/GROUPS/PHYSICS/PAPERS/SUSY-2013-02/'
info.arxiv = 'http://arxiv.org/abs/1405.7875'
info.prettyName =  "jets and met"
info.supersedes = 'ATLAS-CONF-2013-047'




### T1: official maps from ATLAS (they provide 15 SRs , but in the recast (MA5) we have only 10 )

#+++++++ dataset block ++++++++++++++
dataset = DataSetInput("SR6jm")
dataset.setInfo(dataType = 'efficiencyMap', dataId = "SR6jm", observedN = 39, expectedBG = 33 , bgError = 6, upperLimit = '1.1173E+00*fb', expectedUpperLimit = '8.6116E-01*fb')
#+++++++ next txName block ++++++++++++++
T1 = dataset.addTxName('T1')
T1.constraint ="[[['jet','jet']],[['jet','jet']]]"
T1.conditionDescription ="None"
T1.condition ="None"
T1.source = 'ATLAS'
#+++++++ next mass plane block ++++++++++++++
T1 = T1.addMassPlane([[x,y]]*2)
T1.figure  = "figaux_43b"
T1.addSource('efficiencyMap',"orig/T1_SR6jm.dat","txt")
T1.addSource('obsExclusion',"orig/exclusion_T1.txt", "txt")
T1.figureUrl  = "https://atlas.web.cern.ch/Atlas/GROUPS/PHYSICS/PAPERS/SUSY-2013-02/figaux_43b.png"
T1.dataUrl  = "http://hepdata.cedar.ac.uk/view/ins1298722/d188"

"""
#+++++++ dataset block ++++++++++++++
dataset = DataSetInput("SR6jtp")
dataset.setInfo(dataType = 'efficiencyMap', dataId = "SR6jtp", observedN = 6, expectedBG = 4.9 , bgError = 1.6, upperLimit = '3.9922E-01*fb', expectedUpperLimit = '3.0218E-01*fb')
#+++++++ next txName block ++++++++++++++
T1 = dataset.addTxName('T1')
T1.constraint ="[[['jet','jet']],[['jet','jet']]]"
T1.conditionDescription ="None"
T1.condition ="None"
T1.source = 'ATLAS'
#+++++++ next mass plane block ++++++++++++++
T1 = T1.addMassPlane([[x,y]]*2)
T1.figure  = "figaux_45b"
T1.addSource('efficiencyMap',"orig/T1_SR6jt+.dat","txt")
T1.addSource('obsExclusion',"orig/exclusion_T1.txt", "txt")
T1.figureUrl  = "https://atlas.web.cern.ch/Atlas/GROUPS/PHYSICS/PAPERS/SUSY-2013-02/figaux_45b.png"
T1.dataUrl  = "http://hepdata.cedar.ac.uk/view/ins1298722/d194"


#+++++++ dataset block ++++++++++++++
dataset = DataSetInput("SR2jt")
dataset.setInfo(dataType = 'efficiencyMap', dataId = "SR2jt", observedN = 133, expectedBG = 125 , bgError = 10, upperLimit = '1.8181E+00*fb', expectedUpperLimit = '1.5124E+00*fb')
#+++++++ next txName block ++++++++++++++
T1 = dataset.addTxName('T1')
T1.constraint ="[[['jet','jet']],[['jet','jet']]]"
T1.conditionDescription ="None"
T1.condition ="None"
T1.source = 'ATLAS'
#+++++++ next mass plane block ++++++++++++++
T1 = T1.addMassPlane([[x,y]]*2)
T1.figure  = "figaux_33b"
T1.addSource('efficiencyMap',"orig/T1_SR2jt.dat","txt")
T1.addSource('obsExclusion',"orig/exclusion_T1.txt", "txt")
T1.figureUrl  = "https://atlas.web.cern.ch/Atlas/GROUPS/PHYSICS/PAPERS/SUSY-2013-02/figaux_33b.png"
T1.dataUrl  = "http://hepdata.cedar.ac.uk/view/ins1298722/d158"

#+++++++ dataset block ++++++++++++++
dataset = DataSetInput("SR5j")
dataset.setInfo(dataType = 'efficiencyMap', dataId = "SR5j", observedN = 121, expectedBG = 126 , bgError = 13, upperLimit = '1.5429E+00*fb', expectedUpperLimit = '1.7138E+00*fb')
#+++++++ next txName block ++++++++++++++
T1 = dataset.addTxName('T1')
T1.constraint ="[[['jet','jet']],[['jet','jet']]]"
T1.conditionDescription ="None"
T1.condition ="None"
T1.source = 'ATLAS'
#+++++++ next mass plane block ++++++++++++++
T1 = T1.addMassPlane([[x,y]]*2)
T1.figure  = "figaux_41b"
T1.addSource('efficiencyMap',"orig/T1_SR5j.dat","txt")
T1.addSource('obsExclusion',"orig/exclusion_T1.txt", "txt")
T1.figureUrl  = "https://atlas.web.cern.ch/Atlas/GROUPS/PHYSICS/PAPERS/SUSY-2013-02/figaux_41b.png"
T1.dataUrl  = "http://hepdata.cedar.ac.uk/view/ins1298722/d182"
"""

#+++++++ dataset block ++++++++++++++
dataset = DataSetInput("SR2jW")
dataset.setInfo(dataType = 'efficiencyMap', dataId = "SR2jW", observedN = 0, expectedBG = 2.3 , bgError = 1.4, upperLimit = '1.4709E-01*fb', expectedUpperLimit = '2.5070E-01*fb')
#+++++++ next txName block ++++++++++++++
T1 = dataset.addTxName('T1')
T1.constraint ="[[['jet','jet']],[['jet','jet']]]"
T1.conditionDescription ="None"
T1.condition ="None"
T1.source = 'ATLAS'
#+++++++ next mass plane block ++++++++++++++
T1 = T1.addMassPlane([[x,y]]*2)
T1.figure  = "figaux_34b"
T1.addSource('efficiencyMap',"orig/T1_SR2jW.dat","txt")
T1.addSource('obsExclusion',"orig/exclusion_T1.txt", "txt")
T1.figureUrl  = "https://atlas.web.cern.ch/Atlas/GROUPS/PHYSICS/PAPERS/SUSY-2013-02/figaux_34b.png"
T1.dataUrl  = "http://hepdata.cedar.ac.uk/view/ins1298722/d161"


#+++++++ dataset block ++++++++++++++
dataset = DataSetInput("SR4jW")
dataset.setInfo(dataType = 'efficiencyMap', dataId = "SR4jW", observedN = 16, expectedBG = 14 , bgError = 4, upperLimit = '6.7961E-01*fb', expectedUpperLimit = '5.9339E-01*fb')
#+++++++ next txName block ++++++++++++++
T1 = dataset.addTxName('T1')
T1.constraint ="[[['jet','jet']],[['jet','jet']]]"
T1.conditionDescription ="None"
T1.condition ="None"
T1.source = 'ATLAS'
#+++++++ next mass plane block ++++++++++++++
T1 = T1.addMassPlane([[x,y]]*2)
T1.figure  = "figaux_36b"
T1.addSource('efficiencyMap',"orig/T1_SR4jW.dat","txt")
T1.addSource('obsExclusion',"orig/exclusion_T1.txt", "txt")
T1.figureUrl  = "https://atlas.web.cern.ch/Atlas/GROUPS/PHYSICS/PAPERS/SUSY-2013-02/figaux_36b.png"
T1.dataUrl  = "http://hepdata.cedar.ac.uk/view/ins1298722/d167"

"""
#+++++++ dataset block ++++++++++++++
dataset = DataSetInput("SR4jt")
dataset.setInfo(dataType = 'efficiencyMap', dataId = "SR4jt", observedN = 0, expectedBG = 2.5 , bgError = 1.0, upperLimit = '1.4949E-01*fb', expectedUpperLimit = '2.4033E-01*fb')
#+++++++ next txName block ++++++++++++++
T1 = dataset.addTxName('T1')
T1.constraint ="[[['jet','jet']],[['jet','jet']]]"
T1.conditionDescription ="None"
T1.condition ="None"
T1.source = 'ATLAS'
#+++++++ next mass plane block ++++++++++++++
T1 = T1.addMassPlane([[x,y]]*2)
T1.figure  = "figaux_40b"
T1.addSource('efficiencyMap',"orig/T1_SR4jt.dat","txt")
T1.addSource('obsExclusion',"orig/exclusion_T1.txt", "txt")
T1.figureUrl  = "https://atlas.web.cern.ch/Atlas/GROUPS/PHYSICS/PAPERS/SUSY-2013-02/figaux_40b.png"
T1.dataUrl  = "http://hepdata.cedar.ac.uk/view/ins1298722/d179"
"""


#+++++++ dataset block ++++++++++++++
dataset = DataSetInput("SR2jl")
dataset.setInfo(dataType = 'efficiencyMap', dataId = "SR2jl", observedN = 12315, expectedBG = 13000 , bgError = 1000, upperLimit = '7.7800E+01*fb', expectedUpperLimit = '9.7112E+01*fb')
#+++++++ next txName block ++++++++++++++
T1 = dataset.addTxName('T1')
T1.constraint ="[[['jet','jet']],[['jet','jet']]]"
T1.conditionDescription ="None"
T1.condition ="None"
T1.source = 'ATLAS'
#+++++++ next mass plane block ++++++++++++++
T1 = T1.addMassPlane([[x,y]]*2)
T1.figure  = "figaux_31b"
T1.addSource('efficiencyMap',"orig/T1_SR2jl.dat","txt")
T1.addSource('obsExclusion',"orig/exclusion_T1.txt", "txt")
T1.figureUrl  = "https://atlas.web.cern.ch/Atlas/GROUPS/PHYSICS/PAPERS/SUSY-2013-02/figaux_31b.png"
T1.dataUrl  = "http://hepdata.cedar.ac.uk/view/ins1298722/d152"


"""
#+++++++ dataset block ++++++++++++++
dataset = DataSetInput("SR2jm")
dataset.setInfo(dataType = 'efficiencyMap', dataId = "SR2jm", observedN = 715, expectedBG = 760 , bgError = 50, upperLimit = '4.2419E+00*fb', expectedUpperLimit = '5.5524E+00*fb')
#+++++++ next txName block ++++++++++++++
T1 = dataset.addTxName('T1')
T1.constraint ="[[['jet','jet']],[['jet','jet']]]"
T1.conditionDescription ="None"
T1.condition ="None"
T1.source = 'ATLAS'
#+++++++ next mass plane block ++++++++++++++
T1 = T1.addMassPlane([[x,y]]*2)
T1.figure  = "figaux_32b"
T1.addSource('efficiencyMap',"orig/T1_SR2jm.dat","txt")
T1.addSource('obsExclusion',"orig/exclusion_T1.txt", "txt")
T1.figureUrl  = "https://atlas.web.cern.ch/Atlas/GROUPS/PHYSICS/PAPERS/SUSY-2013-02/figaux_32b.png"
T1.dataUrl  = "http://hepdata.cedar.ac.uk/view/ins1298722/d155"

#+++++++ dataset block ++++++++++++++
dataset = DataSetInput("SR4jlm")
dataset.setInfo(dataType = 'efficiencyMap', dataId = "SR4jlm", observedN = 2169, expectedBG = 2120 , bgError = 110, upperLimit = '1.3292E+01*fb', expectedUpperLimit = '1.1561E+01*fb')
#+++++++ next txName block ++++++++++++++
T1 = dataset.addTxName('T1')
T1.constraint ="[[['jet','jet']],[['jet','jet']]]"
T1.conditionDescription ="None"
T1.condition ="None"
T1.source = 'ATLAS'
#+++++++ next mass plane block ++++++++++++++
T1 = T1.addMassPlane([[x,y]]*2)
T1.figure  = "figaux_37b"
T1.addSource('efficiencyMap',"orig/T1_SR4jl-.dat","txt")
T1.addSource('obsExclusion',"orig/exclusion_T1.txt", "txt")
T1.figureUrl  = "https://atlas.web.cern.ch/Atlas/GROUPS/PHYSICS/PAPERS/SUSY-2013-02/figaux_37b.png"
T1.dataUrl  = "http://hepdata.cedar.ac.uk/view/ins1298722/d170"

#+++++++ dataset block ++++++++
dataset = DataSetInput("SR4jl")
dataset.setInfo(dataType = 'efficiencyMap', dataId = "SR4jl", observedN = 608, expectedBG = 630 , bgError = 50, upperLimit = '4.7487E+00*fb', expectedUpperLimit = '5.4345E+00*fb')
#+++++++ next txName block ++++++++++++++
T1 = dataset.addTxName('T1')
T1.constraint ="[[['jet','jet']],[['jet','jet']]]"
T1.conditionDescription ="None"
T1.condition ="None"
T1.source = 'ATLAS'
#+++++++ next mass plane block ++++++++++++++
T1 = T1.addMassPlane([[x,y]]*2)
T1.figure  = "figaux_38b"
T1.addSource('efficiencyMap',"orig/T1_SR4jl.dat","txt")
T1.addSource('obsExclusion',"orig/exclusion_T1.txt", "txt")
T1.figureUrl  = "https://atlas.web.cern.ch/Atlas/GROUPS/PHYSICS/PAPERS/SUSY-2013-02/figaux_38b.png"
T1.dataUrl  = "http://hepdata.cedar.ac.uk/view/ins1298722/d173"

"""

#+++++++ dataset block ++++++++++++++
dataset = DataSetInput("SR4jm")
dataset.setInfo(dataType = 'efficiencyMap', dataId = "SR4jm", observedN = 24, expectedBG = 37 , bgError = 6, upperLimit = '5.0301E-01*fb', expectedUpperLimit = '8.8617E-01*fb')
#+++++++ next txName block ++++++++++++++
T1 = dataset.addTxName('T1')
T1.constraint ="[[['jet','jet']],[['jet','jet']]]"
T1.conditionDescription ="None"
T1.condition ="None"
T1.source = 'ATLAS'
#+++++++ next mass plane block ++++++++++++++
T1 = T1.addMassPlane([[x,y]]*2)
T1.figure  = "figaux_39b"
T1.addSource('efficiencyMap',"orig/T1_SR4jm.dat","txt")
T1.addSource('obsExclusion',"orig/exclusion_T1.txt", "txt")
T1.figureUrl  = "https://atlas.web.cern.ch/Atlas/GROUPS/PHYSICS/PAPERS/SUSY-2013-02/figaux_39b.png"
T1.dataUrl  = "http://hepdata.cedar.ac.uk/view/ins1298722/d176"

"""
#+++++++ dataset block ++++++++++++++
dataset = DataSetInput("SR3j")
dataset.setInfo(dataType = 'efficiencyMap', dataId = "SR3j", observedN = 7, expectedBG = 5 , bgError = 1.2, upperLimit = '4.3344E-01*fb', expectedUpperLimit = '3.3172E-01*fb')
#+++++++ next txName block ++++++++++++++
T1 = dataset.addTxName('T1')
T1.constraint ="[[['jet','jet']],[['jet','jet']]]"
T1.conditionDescription ="None"
T1.condition ="None"
T1.source = 'ATLAS'
#+++++++ next mass plane block ++++++++++++++
T1 = T1.addMassPlane([[x,y]]*2)
T1.figure  = "figaux_35b"
T1.addSource('efficiencyMap',"orig/T1_SR3j.dat","txt")
T1.addSource('obsExclusion',"orig/exclusion_T1.txt", "txt")
T1.figureUrl  = "https://atlas.web.cern.ch/Atlas/GROUPS/PHYSICS/PAPERS/SUSY-2013-02/figaux_35b.png"
T1.dataUrl  = "http://hepdata.cedar.ac.uk/view/ins1298722/d164"

#+++++++ dataset block ++++++++++++++
dataset = DataSetInput("SR6jt")
dataset.setInfo(dataType = 'efficiencyMap', dataId = "SR6jt", observedN = 5, expectedBG = 5.2 , bgError = 1.4, upperLimit = '3.3159E-01*fb', expectedUpperLimit = '3.3330E-01*fb')
#+++++++ next txName block ++++++++++++++
T1 = dataset.addTxName('T1')
T1.constraint ="[[['jet','jet']],[['jet','jet']]]"
T1.conditionDescription ="None"
T1.condition ="None"
T1.source = 'ATLAS'
#+++++++ next mass plane block ++++++++++++++
T1 = T1.addMassPlane([[x,y]]*2)
T1.figure  = "figaux_44b"
T1.addSource('efficiencyMap',"orig/T1_SR6jt.dat","txt")
T1.addSource('obsExclusion',"orig/exclusion_T1.txt", "txt")
T1.figureUrl  = "https://atlas.web.cern.ch/Atlas/GROUPS/PHYSICS/PAPERS/SUSY-2013-02/figaux_44b.png"
T1.dataUrl  = "http://hepdata.cedar.ac.uk/view/ins1298722/d191"

"""


def add ( dataset ):

    dataset_n = dataset._name.replace('SR','').replace('6jt+','6jtp').replace('4jl-','4jlm')

    #### T2
    T2 = dataset.addTxName('T2')
    T2.checked =''
    T2.constraint ="[[['jet']],[['jet']]]"
    T2.conditionDescription ="None"
    T2.condition ="None"
    T2.source = 'SModelS'
    #+++++++ next mass plane block ++++++++++++++
    t2FigureUrls = { "6jt": "28b", "6jl": "26b", "3j": "19b", "4jl": "22b", "4jt": "24b", "5j": "25b", "4jlm": "23b", "2jm": "16b", "2jt": "17b" }
    T2qq = T2.addMassPlane([[x,y]]*2)
    T2qq.figure    = None
    if dataset_n in t2FigureUrls.keys():
        T2qq.figureUrl = "https://atlas.web.cern.ch/Atlas/GROUPS/PHYSICS/PAPERS/SUSY-2013-02/figaux_%s.png" % t2FigureUrls[dataset_n]
    T2qq.dataUrl   = None
    T2qq.addSource('obsExclusion',"orig/exclusion_T2.txt", "txt")
    T2qq.addSource('efficiencyMap',"orig/T2/MA5_EM_T2_Results_%s.dat"% dataset_n, "txt" , objectName ="None", index = None)

    ### T5
    T5 = dataset.addTxName('T5')
    T5.checked = ''
    T5.constraint ="[[['jet'],['jet']],[['jet'],['jet']]]"
    T5.conditionDescription ="None"
    T5.condition ="None"
    T5.massConstraint = None
    T5.source = 'SModelS'
    T5.dataUrl = None
    T5.validated = "N/A"
    T5_x005 = T5.addMassPlane( [[x,0.05*x + 0.95*y,y]]*2 )
    T5_x005.addSource ( "efficiencyMap", "orig/MA5_T5_TGQ/T5_X005_atlas_2013_02/MA5_EM_T5_MAPS_%s.dat" % dataset_n, "txt", objectName ="None", index = None )
    T5_x05 = T5.addMassPlane( [[x,0.5*x + 0.5*y,y]]*2 )
    T5_x05.addSource ( "efficiencyMap", "orig/MA5_T5_TGQ/T5_X05_atlas_2013_02/MA5_EM_T5_MAPS_%s.dat" % dataset_n, "txt", objectName ="None", index = None )
    T5_x095 = T5.addMassPlane( [[x,0.95*x + 0.05*y,y]]*2 )
    T5_x095.addSource ( "efficiencyMap", "orig/MA5_T5_TGQ/T5_X095_atlas_2013_02/MA5_EM_T5_MAPS_%s.dat" % dataset_n, "txt", objectName ="None", index = None )
    T5_DeltaGluSq_5 = T5.addMassPlane( [[x, x-5 ,y]]*2 )
    T5_DeltaGluSq_5.addSource ( "efficiencyMap", "orig/MA5_T5_TGQ/T5_DeltaGluSq_5_atlas_2013_02/MA5_EM_T5_MAPS_%s.dat" % dataset_n, "txt", objectName ="None", index = None )
    T5_DeltaSqNeu_5 = T5.addMassPlane( [[x, y+5 ,y]]*2 )
    T5_DeltaSqNeu_5.addSource ( "efficiencyMap", "orig/MA5_T5_TGQ/T5_DeltaSqNeu_5_atlas_2013_02/MA5_EM_T5_MAPS_%s.dat" % dataset_n, "txt", objectName ="None", index = None )

    ### TGQ
    TGQ = dataset.addTxName('T3GQ')
    TGQ.checked =''
    TGQ.constraint ="[[['jet']],[['jet'],['jet']]]"
    TGQ.conditionDescription ="None"
    TGQ.condition ="None"
    TGQ.source = 'SModelS'
    TGQ.validated = "N/A"
    TGQ.dataUrl = None
    for num in  [name.split('_')[1] for name in os.listdir('orig/MA5_T5_TGQ') if 'TGQ' in name]: # listing the gluino masses from the name in the directories
        a = TGQ.addMassPlane( [ [x,y], [eval(num),x,y] ] )
        a.addSource('efficiencyMap',"orig/MA5_T5_TGQ/TGQ_"+num+"_atlas_2013_02/MA5_EM_TGQon_MAPS_%s.dat" % dataset_n, "txt", objectName ="None", index = None )




# commented SRs do not have the recast implementation so they are skipped

datasets = { "SR6jl" : ( 121,111 ,11 , '1.9230E+00*fb', '1.5312E+00*fb'),
             #"SR6jm": (39,36,3      , upperLimit = '1.1173E+00*fb', expectedUpperLimit = '8.6116E-01*fb'),
             #"SR2jW": ( 0 , 2.3 , 1.4 , upperLimit = '1.4709E-01*fb', expectedUpperLimit = '2.5070E-01*fb'),
             #"SR4jW": (14 , 14 , 4 , upperLimit = '6.7961E-01*fb', expectedUpperLimit = '5.9339E-01*fb'),
             #"SR2jl" : (12315 , 13000 , 1000 , upperLimit = '7.7800E+01*fb', expectedUpperLimit = '9.7112E+01*fb'),
             "SR6jtp": (6 , 4.9 , 1.6 , '3.9922E-01*fb', '3.0218E-01*fb'),
             "SR2jt" : (133, 125 , 10, '1.8181E+00*fb',   '1.5124E+00*fb'),
             "SR5j"  : (121, 126, 13 , '1.5429E+00*fb',   '1.7138E+00*fb'),
             "SR4jt" : (0, 2.5, 1.0 , '1.4949E-01*fb',    '2.4033E-01*fb'),
             "SR2jm" : (715,760,50, '4.2419E+00*fb',      '5.5524E+00*fb'),
             "SR4jlm": (2169, 2120, 110, '1.3292E+01*fb', '1.1561E+01*fb'),
             "SR4jl" : (608 , 630, 50, '4.7487E+00*fb',   '5.4345E+00*fb'),
             #"SR4jm" : (24, 37, 6 , upperLimit = '5.0301E-01*fb', expectedUpperLimit = '8.8617E-01*fb'),
             "SR3j": (7 , 5 , 1.2 ,  '4.3344E-01*fb', '3.3172E-01*fb'),
             "SR6jt": (5 , 5.2 , 1.4 , '3.3159E-01*fb', '3.3330E-01*fb')

             }

dses = {}

for name, numbers in datasets.items():
    #+++++++ dataset block ++++++++++++++
    dataset = DataSetInput( name )
    dses[name] = dataset
    name = name.replace('SR_','')
    dataset.setInfo(dataType = 'efficiencyMap', dataId = name,
                    observedN = numbers[0], expectedBG = numbers[1], bgError = numbers[2],
                    upperLimit = numbers[3] , expectedUpperLimit = numbers[4] )
    add ( dataset )

T1 = dses["SR6jtp"].addTxName('T1')
T1.constraint ="[[['jet','jet']],[['jet','jet']]]"
T1.conditionDescription ="None"
T1.condition ="None"
T1.source = 'ATLAS'
#+++++++ next mass plane block ++++++++++++++
T1 = T1.addMassPlane([[x,y]]*2)
T1.figure  = "figaux_45b"
T1.addSource('efficiencyMap',"orig/T1_SR6jt+.dat","txt")
T1.addSource('obsExclusion',"orig/exclusion_T1.txt", "txt")
T1.figureUrl  = "https://atlas.web.cern.ch/Atlas/GROUPS/PHYSICS/PAPERS/SUSY-2013-02/figaux_45b.png"
T1.dataUrl  = "http://hepdata.cedar.ac.uk/view/ins1298722/d194"

T1 = dses["SR6jl"].addTxName('T1')
T1.constraint ="[[['jet','jet']],[['jet','jet']]]"
T1.conditionDescription ="None"
T1.condition ="None"
T1.source = 'ATLAS'
#+++++++ next mass plane block ++++++++++++++
T1 = T1.addMassPlane([[x,y]]*2)
T1.figure  = "figaux_42b"
T1.addSource('efficiencyMap',"orig/T1_SR6jl.dat","txt")
T1.addSource('obsExclusion',"orig/exclusion_T1.txt", "txt")
T1.figureUrl  = "https://atlas.web.cern.ch/Atlas/GROUPS/PHYSICS/PAPERS/SUSY-2013-02/figaux_42b.png"
T1.dataUrl  = "http://hepdata.cedar.ac.uk/view/ins1298722/d185"

#+++++++ next txName block ++++++++++++++
T1 = dses["SR4jlm"].addTxName('T1')
T1.constraint ="[[['jet','jet']],[['jet','jet']]]"
T1.conditionDescription ="None"
T1.condition ="None"
T1.source = 'ATLAS'
#+++++++ next mass plane block ++++++++++++++
T1 = T1.addMassPlane([[x,y]]*2)
T1.figure  = "figaux_37b"
T1.addSource('efficiencyMap',"orig/T1_SR4jl-.dat","txt")
T1.addSource('obsExclusion',"orig/exclusion_T1.txt", "txt")
T1.figureUrl  = "https://atlas.web.cern.ch/Atlas/GROUPS/PHYSICS/PAPERS/SUSY-2013-02/figaux_37b.png"
T1.dataUrl  = "http://hepdata.cedar.ac.uk/view/ins1298722/d170"

databaseCreator.create()
