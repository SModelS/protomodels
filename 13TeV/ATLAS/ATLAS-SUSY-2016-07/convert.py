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
help = 'path to the package smodels_utils',\
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



#+++++++ global info block ++++++++++++++
info = MetaInfoInput('ATLAS-SUSY-2016-07')
info.url = "https://atlas.web.cern.ch/Atlas/GROUPS/PHYSICS/PAPERS/SUSY-2016-07/"
info.sqrts = 13
info.lumi = 36.1
info.prettyName = "0L + jets + Etmiss"
info.private = False
info.arxiv = 'https://arxiv.org/abs/1712.02332'
info.contact = 'atlas-phys-susy-conveners@cern.ch'
info.publication = 'Phys. Rev. D 97, 112001 (2018)'

#+++++++ dataset block ++++++++++++++
dataset = DataSetInput('data')
dataset.setInfo(dataType = 'upperLimit', dataId = None)

topos = [ "T1", "T2", "T6WW", "T5WW", "T5ZZ", "T5WZh", "T6WZh", "T5WWoff", "T6WWoff" ]
# topos = [ "T5WZh", "T6WZh" ]

constraints = { "T2": "[[[jet]],[[jet]]]", "T1": "[[[jet,jet]],[[jet,jet]]]", "T6WW": "[[[jet],[W]],[[jet],[W]]]", "T6WWoff": "2.23 * [[[jet],[jet,jet]],[[jet],[jet,jet]]]", "T5WW": "[[[jet,jet],[W]],[[jet,jet],[W]]]", "T5ZZ": "[[[jet,jet],[Z]],[[jet,jet],[Z]]]", "T5WZh": "[[[jet,jet],[Z]],[[jet,jet],[W]]]+[[[jet,jet],[higgs]],[[jet,jet],[W]]]+[[[jet,jet],[W]],[[jet,jet],[W]]]+[[[jet,jet],[Z]],[[jet,jet],[Z]]]+[[[jet,jet],[higgs]],[[jet,jet],[higgs]]]+[[[jet,jet],[higgs]],[[jet,jet],[Z]]]", "T6WZh": "[[[jet],[Z]],[[jet],[W]]]+[[[jet],[higgs]],[[jet],[W]]]+[[[jet],[W]],[[jet],[W]]]+[[[jet],[Z]],[[jet],[Z]]]+[[[jet],[higgs]],[[jet],[higgs]]]+[[[jet],[Z]],[[jet],[higgs]]]", "T5WWoff": "2.23*[[[jet,jet],[jet,jet]],[[jet,jet],[jet,jet]]]" }

conditions = { "T5WZh": "Csim([[[jet,jet],[Z]],[[jet,jet],[W]]]+[[[jet,jet],[higgs]],[[jet,jet],[W]]],2.*([[[jet,jet],[W]],[[jet,jet],[W]]]))", "T6WZh": "Csim([[[jet],[Z]],[[jet],[W]]]+[[[jet],[higgs]],[[jet],[W]]],2.*([[[jet],[W]],[[jet],[W]]]))" }

figure = { "T2": "Fig. 76a", "T1": "Fig. 76b", "T6WW": ('Fig. 77a', 'Fig. 77b' ), "T5WW": ( 'Fig. 77c', 'Fig. 77d' ), "T5ZZ": ('Fig.79', None) , "T5WZh": ('Fig. 78b', None ), "T6WZh":  ('Fig. 78a', None ), "T5WWoff": ( 'Fig. 77c', 'Fig. 77d' ), "T6WWoff": ('Fig. 77a', 'Fig. 77b' ) }

baseUrl = 'https://atlas.web.cern.ch/Atlas/GROUPS/PHYSICS/PAPERS/SUSY-2016-07/'
figureUrl = { "T2": 'figaux_076a.png', "T1": 'figaux_076b.png',
              "T6WW": ( 'figaux_077a.png', 'figaux_077b.png' ),
              "T5WW": ( 'figaux_077c.png', 'figaux_077d.png' ),
              "T5WWoff": ( 'figaux_077c.png', 'figaux_077d.png' ),
              "T5WZh": ('figaux_078b.png', None),
              "T6WZh": ('figaux_078a.png', None ),
              "T5ZZ": ('figaux_079.png', None),
              "T6WWoff": ( 'figaux_077a.png', 'figaux_077b.png' ),
              "TGQ": ( 'figaux_080a.png', 'figaux_80b.png', 'figaux_80c.png' )
}
dataBaseUrl = "'https://www.hepdata.net/record/ins1641270?version=5&table=X-section U.L. & best SR @@NUM@@  : @@MEFFRJR@@'"
exclusionIndex = { "T2": 1, "T1": 2, "T6WW": ( 3, 4 ), "T5WW": (5, 6 ),
                   "T5ZZ": (7, None), "T6WWoff": ( 3, 4 ),
                   "T5WZh": (9,None), "T6WZh": (8, None ), "T5WWoff": (5, 6 ),
                   "TGQ": (10,11,12) }
exclusionLineIndex = { "T2": 1, "T1": 2, "T6WW": ( 3, 4 ), "T5WW": (5, 6 ),
                   "T5ZZ": (7, None), "T6WWoff": ( 13, 14 ),
                   "T5WZh": (9,None), "T6WZh": (8, None ), "T5WWoff": (15, 16 ),
                   "TGQ": (10,11,12) }

massPlanes = { "T2": 2*[[x, y]], "T1": 2*[[x, y]],
        "T6WW": ( 2*[[x, 0.5*(x+y), y]], 2*[[x, y , 60.]] ),
        "T6WWoff": ( 2*[[x, 0.5*(x+y), y]], 2*[[x, y , 60.]] ),
        "T5WZh": (2*[[x, y , 60.]],None),
        "T6WZh": (2*[[x, y , 60.]], None ), 
        "T5ZZ": (2*[[x, y, 1.]], None),
        "T5WW": ( 2*[[x, 0.5*(x+y), y]], 2*[[x, y , 60.]] ),
        "T5WWoff": ( 2*[[x, 0.5*(x+y), y]], 2*[[x, y, 60.]] )
#        "T5WW": ( 2*[[x, 0.5*(x+y), y]], 2*[[x, 60. + y * (x - 60.), 60.]] ),
#        "T5WWoff": ( 2*[[x, 0.5*(x+y), y]], 2*[[x, 60. + y * (x - 60.), 60.]] )
#        "T6WW": ( 2*[[x, 0.5*(x+y), y]], 2*[[x, 60. + y * (x - 60.), 60.]] ),
#        "T6WWoff": ( 2*[[x, 0.5*(x+y), y]], 2*[[x, 60. + y * (x - 60.), 60.]] ),
#        "T5WZh": (2*[[x, 60. + y * (x - 60.), 60.]],None),
#        "T6WZh": (2*[[x, 60. + y * (x - 60.), 60.]], None ), 
}

rjrs = [ 1, 2, 3, 5 ]

for topo in topos:
    #+++++++ next txName block ++++++++++++++
    Tx = dataset.addTxName( topo )
    Tx.checked = ''
    Tx.constraint = constraints[topo]
    Tx.conditionDescription = None
    Tx.condition = None
    if topo in [ "T5WWoff", "T6WWoff" ]:
            Tx.massConstraint = [['dm >= 0.0','dm <= 76.']]*2
    if topo in [ "T6WW", "T6WWoff" ]:
        Tx.validationTarball = topo+"left.tar.gz"
    if topo in conditions.keys():
        Tx.condition = conditions[topo]
    Tx.source = "ATLAS"
    mPlanes = massPlanes[topo]
    if type (mPlanes) == list:
        Tx_1 = Tx.addMassPlane(massPlanes[topo])
        Tx_1.figure = figure[topo]
        Tx_1.figureUrl = baseUrl + figureUrl[topo]
        eIdx = exclusionIndex[topo]
        eLIdx = exclusionLineIndex[topo]
        meffrjr = "Meff"
        meffrjrurl = "Meff"
        if eIdx in rjrs:
            meffrjr = "Meff+RJR"
            meffrjrurl = "Meff%2BRJR"
        Tx_1.dataUrl = dataBaseUrl.replace("@@NUM@@",str(eIdx) ).replace("@@MEFFRJR@@",meffrjrurl )

        ulfile = 'orig/X-sectionU.L.&bestSR%d:%s.csv' % ( eIdx, meffrjr )
        expCont = 'orig/Exclusioncontour(exp.)%d:%s.csv' % (eLIdx,meffrjr)
        obsCont = 'orig/Exclusioncontour(obs.)%d:%s.csv' % (eLIdx,meffrjr)
        units = [ None, None, 'fb' ]
        if "60." in str(massPlanes[topo]):
            units = [ ("GeV","X:60"), ("GeV","X:60" ), None ]
        Tx_1.setSources(dataLabels= ['expExclusion', 'obsExclusion', 'upperLimits'],
            objectNames = [ None, None, None ],
            dataFiles= [ expCont, obsCont, ulfile ],
            units = units,
            dataFormats= ['csv', 'csv', 'csv'])
    if type (mPlanes) == tuple:
        eIdx = exclusionIndex[topo][0]
        eLIdx = exclusionLineIndex[topo][0]
        for ctr,mp in enumerate(mPlanes):
            if mp == None:
                continue
            Tx_m = Tx.addMassPlane(mp)
            Tx_m.figure = figure[topo][0]
            if ctr == 0:
                    Tx_m.figureUrl = baseUrl + figureUrl[topo][0]
            else:
                pass
            eIdx = exclusionIndex[topo][ctr]
            eLIdx = exclusionLineIndex[topo][ctr]
            Tx_m.dataUrl = dataBaseUrl.replace("@@NUM@@",str(eIdx) ).replace("@@MEFFRJR@@",meffrjrurl )

            meffrjr = "Meff"
            if eIdx in rjrs:
                meffrjr = "Meff+RJR"
            ulfile = 'orig/X-sectionU.L.&bestSR%d:%s.csv' % ( eIdx, meffrjr )
            expCont = 'orig/Exclusioncontour(exp.)%d:%s.csv' % (eLIdx,meffrjr)
            obsCont = 'orig/Exclusioncontour(obs.)%d:%s.csv' % (eLIdx,meffrjr)
            units = [ None, None, 'fb' ]
            if "60." in str(mp):
                units = [ ("GeV","X:60"), ("GeV","X:60" ), None ]
            if eIdx in [ 4, 6, 8, 9 ]:
                ulfile = 'orig/UL_SR%d.csv' % eIdx
                expCont = 'orig/expSR%d.csv' % eIdx
                obsCont = 'orig/obsSR%d.csv' % eIdx
                units = [ None, None, 'fb' ]
            Tx_m.setSources(dataLabels= ['expExclusion', 'obsExclusion', 'upperLimits' ],
                objectNames = [ None, None, None ],
                dataFiles= [ expCont, obsCont, ulfile ],
                units = units,
                dataFormats= ['csv', 'csv', 'csv' ])

databaseCreator.create()
