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
from smodels_utils.dataPreparation.inputObjects import MetaInfoInput,DataSetInput, getSignalRegionsEMBaked, getStatsEMBaked
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
info.contact = 'smodels-users@lists.oeaw.ac.at'
info.publication = 'Phys. Rev. D 97, 112001 (2018)'
info.comment = 'Recast with MadAnalysis5, PAD code http://doi.org/10.7484/INSPIREHEP.DATA.56DC.PPE2'

stats = getStatsEMBaked ()

topos = []
topos = [ "T1", "T2", "T6WW", "T5WW", "T5ZZ", "T5WZh", "T6WZh", "T5WWoff", "T6WWoff", "TGQ" ]
topos.append ( "T3GQ" )
topos.append ( "T5GQ" )
# topos = [ "T6WWoff" ]
constraints = { "T2": "[[[jet]],[[jet]]]", "T1": "[[[jet,jet]],[[jet,jet]]]", "T6WW": "[[[jet],[W]],[[jet],[W]]]", "T6WWoff": "2.23 * [[[jet],[jet,jet]],[[jet],[jet,jet]]]", "T5WW": "[[[jet,jet],[W]],[[jet,jet],[W]]]", "T5ZZ": "[[[jet,jet],[Z]],[[jet,jet],[Z]]]", "T5WZh": "[[[jet,jet],[Z]],[[jet,jet],[W]]]+[[[jet,jet],[higgs]],[[jet,jet],[W]]]+[[[jet,jet],[higgs]],[[jet,jet],[higgs]]]+[[[jet,jet],[higgs]],[[jet,jet],[Z]]]", "T6WZh": "[[[jet],[Z]],[[jet],[W]]]+[[[jet],[higgs]],[[jet],[W]]]+[[[jet],[Z]],[[jet],[Z]]]+[[[jet],[higgs]],[[jet],[higgs]]]+[[[jet],[Z]],[[jet],[higgs]]]", "T5WWoff": "2.23*[[[jet,jet],[jet,jet]],[[jet,jet],[jet,jet]]]", "TGQ": "[[[jet]],[[jet,jet]]]", "T3GQ": "[[[jet]],[[jet],[jet]]]", "T5GQ": "[[[jet],[jet,jet]],[[jet,jet]]]" }
conditions = { "T5WZh": "Csim([[[jet,jet],[Z]],[[jet,jet],[W]]]+[[[jet,jet],[higgs]],[[jet,jet],[W]]],2.*([[[jet,jet],[W]],[[jet,jet],[W]]]))", "T6WZh": "Csim([[[jet],[Z]],[[jet],[W]]]+[[[jet],[higgs]],[[jet],[W]]],2.*([[[jet],[W]],[[jet],[W]]]))" }
figure = { "T2": "Fig. 76a", "T1": "Fig. 76b", "T6WW": ('Fig. 77a', 'Fig. 77b' ), "T5WW": ( 'Fig. 77c', 'Fig. 77d' ), "T5ZZ": ('Fig.79', None) , "T5WZh": ('Fig. 78b', None ), "T6WZh":  ('Fig. 78a', None ), "T5WWoff": ( 'Fig. 77c', 'Fig. 77d' ), "T6WWoff": ('Fig. 77a', 'Fig. 77b' ), "TGQ": ('Fig. 80a', 'Fig. 80b', 'Fig. 80c' ), "T3GQ": None, "T5GQ": None }

baseUrl = 'https://atlas.web.cern.ch/Atlas/GROUPS/PHYSICS/PAPERS/SUSY-2016-07/'
figureUrl = { "T2": 'figaux_076a.png', "T1": 'figaux_076b.png',
              "T6WW": ( 'figaux_077a.png', 'figaux_077b.png' ),
              "T5WW": ( 'figaux_077c.png', 'figaux_077d.png' ),
              "T5WWoff": ( 'figaux_077c.png', 'figaux_077d.png' ),
              "T5WZh": ('figaux_078b.png', None),
              "T6WZh": ('figaux_078a.png', None ),
              "T5ZZ": ('figaux_079.png', None),
              "T6WWoff": ( 'figaux_077a.png', 'figaux_077b.png' ),
              "T3GQ": None,
              "T5GQ": None,
#              "T3GQ": ( None, None, None ),
#              "T5GQ": ( None, None, None ),
              "TGQ": ( 'figaux_080a.png', 'figaux_80b.png', 'figaux_80c.png' )
}
dataBaseUrl = "'https://www.hepdata.net/record/ins1641270?version=5&table=X-section U.L. & best SR @@NUM@@  : @@MEFFRJR@@'"
exclusionIndex = { "T2": 1, "T1": 2, "T6WW": ( 3, 4 ), "T5WW": (5, 6 ),
                   "T5ZZ": (7, None), "T6WWoff": ( 13, 14 ),
                   "T5WZh": (9,None), "T6WZh": (8, None ), "T5WWoff": ( 15, 16 ),
                   "TGQ": (10,11,12), "T3GQ": None, "T5GQ": None }
massPlanes = { "T2": 2*[[x, y]], "T1": 2*[[x, y]],
    "T6WW": ( 2*[[x, 0.5*(x+y), y]], 2*[[x, y, 60.]]),
    "T6WWoff": ( 2*[[x, 0.5*(x+y), y]], 2*[[x, y, 60.]]),
    "T5WZh": (2*[[x, y, 60.]],None),
    "T6WZh": (2*[[x, y, 60.]], None ), "T5ZZ": (2*[[x, y, 1.]], None),
    "T5WW": ( 2*[[x, 0.5*(x+y), y]], 2*[[x, y, 60.]]),
    "T5WWoff": ( 2*[[x, 0.5*(x+y), y]], 2*[[x, y, 60.]]),
    "T3GQ": [[y,z],[x,y,z]],
    "T5GQ": [[x,y,z],[y,z]],
#    "T3GQ": ( [[y,0.],[x,y,0.]], [[y,695.],[x,y,695.]], [[y,995.],[x,y,995.]] ),
#    "T5GQ": ( [[x,y,0.],[y, 0.]], [[x, y, 695.],[y, 695.]], [[x, y, 995.],[y, 995.]] ),
    "TGQ": ( [[x,0.],[y, 0.]], [[x, 695.],[y, 695.]], [[x, 995.],[y, 995.]] )
}

datasets={}

dsnames=getSignalRegionsEMBaked("orig/%s.embaked" % "T2" )

for dsname in dsnames:
    #+++++++ dataset block ++++++++++++++
    if not dsname in stats:
        print ( "cannot find stats for %s" % dsname )
        continue
    dataset = DataSetInput( dsname )
    dst = stats[dsname]
    dataset.setInfo(dataType = 'efficiencyMap', dataId = dsname,
         observedN = dst["nobs"], expectedBG = dst["nb"] , bgError = dst["deltanb"] )
    datasets[dsname]=dataset


rjrs = [ 1, 2, 3, 5 ]

for topo in topos:
    for dsname in dsnames:
        #+++++++ next txName block ++++++++++++++
        Tx = datasets[dsname].addTxName( topo )
        Tx.checked = ''
        Tx.constraint = constraints[topo]
        Tx.conditionDescription = None
        Tx.condition = None
        if topo in [ "T5WWoff", "T6WWoff" ]:
            Tx.massConstraint = [['dm >= 0.0','dm <= 76.']]*2
        if topo in [ "T6WZh", "T6WW", "T6WWoff" ]:
            Tx.validationTarball = topo+"left.tar.gz"
        if topo in conditions.keys():
            Tx.condition = conditions[topo]
        Tx.source = "SModelS"
        #+++++++ next mass plane block ++++++++++++++
        mPlanes = massPlanes[topo]
        if type (mPlanes) == list:
            Tx_1 = Tx.addMassPlane(massPlanes[topo])
            Tx_1.figure = figure[topo]
            if figureUrl[topo]!=None:
                Tx_1.figureUrl = baseUrl + figureUrl[topo]
            if topo in [ "TGQ" ]:
                Tx_1.axes = "[[x, 0.0], [y, 0.0]]; [[x, 695.0], [y, 695.0]]; [[x, 995.0], [y, 995.0]]"
            if topo in [ "T3GQ" ]:
                Tx_1.axes = "[[y, 0.0], [x, y, 0.0]]; [[y, 695.0], [x, y, 695.0]]; [[y, 995.0], [x, y, 995.0]]"
            if topo in [ "T5GQ" ]:
                Tx_1.axes = "[[x, y, 0.0], [y, 0.0]]; [[x, y, 695.0], [y, 695.0]]; [[x, y, 995.0], [y, 995.0]]"
            eIdx = exclusionIndex[topo]
            meffrjr = "Meff"
            meffrjrurl = "Meff"
            if eIdx in rjrs:
                meffrjr = "Meff+RJR"
                meffrjrurl = "Meff%2BRJR"
            Tx_1.dataUrl = None # dataBaseUrl.replace("@@NUM@@",str(eIdx) ).replace("@@MEFFRJR@@",meffrjrurl )
            units = [ None, None, None ]
            if "60." in str(massPlanes[topo]):
                units = [ ("GeV","X:60"), ("GeV","X:60" ), None ]
            if eIdx == None:
                Tx_1.setSources(dataLabels= [ 'efficiencyMap'],
                    objectNames = [ dsname ],
                    dataFiles= [ 'orig/%s.embaked' % topo ],
                    units = [ None ],
                    dataFormats= ['embaked'])

            else:
                expCont = 'orig/Exclusioncontour(exp.)%d:%s.csv' % (eIdx,meffrjr)
                obsCont = 'orig/Exclusioncontour(obs.)%d:%s.csv' % (eIdx,meffrjr)
                if eIdx in [ 4, 6, 8, 9 ]:
                    ulfile = 'orig/UL_SR%d.csv' % eIdx
                    expCont = 'orig/expSR%d.csv' % eIdx
                    obsCont = 'orig/obsSR%d.csv' % eIdx
                Tx_1.setSources(dataLabels= ['expExclusion', 'obsExclusion', 'efficiencyMap'],
                    objectNames = [ None, None, dsname ],
                    dataFiles= [ expCont, obsCont, 'orig/%s.embaked' % topo ],
                    units = units,
                    dataFormats= ['csv', 'csv', 'embaked'])
        if type (mPlanes) == tuple:
            mArray=2*[[x,y,z]]
            if topo == "TGQ":
                mArray=[[x,z],[y,z]]
            Tx_1 = Tx.addMassPlane(mArray)
            eIdx = exclusionIndex[topo][0]
            Tx_1.setSources(dataLabels= ['efficiencyMap'],
                objectNames = [ dsname ],
                dataFiles= [ 'orig/%s.embaked' % topo ],
                units = [ None],
                dataFormats= [ 'embaked'] )
            for ctr,mp in enumerate(mPlanes):
                if mp == None:
                    continue
                Tx_m = Tx.addMassPlane(mp)
                Tx_m.figure = figure[topo][0]
                if ctr == 0:
                    Tx_m.figureUrl = baseUrl + figureUrl[topo][0]
                eIdx = exclusionIndex[topo][ctr]
                meffrjr = "Meff"
                meffrjrurl = "Meff"
                if eIdx in rjrs:
                    meffrjr = "Meff+RJR"
                    meffrjrurl = "Meff%2BRJR"

                Tx_m.dataUrl = None # dataBaseUrl.replace("@@NUM@@",str(eIdx) ).replace("@@MEFFRJR@@",meffrjrurl )
                expCont = 'orig/Exclusioncontour(exp.)%d:%s.csv' % (eIdx,meffrjr)
                obsCont = 'orig/Exclusioncontour(obs.)%d:%s.csv' % (eIdx,meffrjr)
                units = [ None, None ]
                if "60." in str(mp):
                    units = [ ("GeV","X:60"), ("GeV","X:60" ) ]
                Tx_m.setSources(dataLabels= ['expExclusion', 'obsExclusion' ],
                    objectNames = [ None, None ],
                    dataFiles= [ expCont, obsCont ],
                    units = units,
                    dataFormats= ['csv', 'csv' ])

databaseCreator.create()
