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
info                 = MetaInfoInput('ATLAS-SUSY-2018-31')
info.url             = 'https://atlas.web.cern.ch/Atlas/GROUPS/PHYSICS/PAPERS/SUSY-2018-31/'
info.sqrts             = 13
info.lumi             = 139.
info.prettyName     = '2b + 2H(bb) + Etmiss'
info.private         = False
info.arxiv             = 'arXiv:1908.03122'
info.contact         = 'atlas-phys-susy-conveners@cern.ch'
info.publication     = ''

# +++++++ dataset block ++++++++++++++
dataset = DataSetInput('data')
dataset.setInfo( dataType = 'upperLimit', dataId = None )

# +++++++ next txName block ++++++++++++++
T6bbHH = dataset.addTxName("T6bbHH" )
T6bbHH.checked              = 'NO'
T6bbHH.constraint           = "[[['b'],['higgs']],[['b'],['higgs']]]"
T6bbHH.conditionDescription = None
T6bbHH.condition            = None
T6bbHH.source               = 'ATLAS'

# +++++++ next mass plane block ++++++++++++++
T6bbHH60           = T6bbHH.addMassPlane( 2*[[ x, y, 60. ]] )
T6bbHH60.figure    = "figaux_01a"
T6bbHH60.figureUrl = "https://atlas.web.cern.ch/Atlas/GROUPS/PHYSICS/PAPERS/SUSY-2018-31/figaux_01a.png"
T6bbHH60.dataUrl   = "https://www.hepdata.net/record/ins1748602?version=1&table=M60_bestCombined_UpperLimit"
T6bbHH60.setSources(dataLabels = [ 'obsExclusion', 'upperLimits' ],
    dataFiles = [ "orig/HEPData-ins1748602-v1-M60_Obs.csv", \
                  "orig/HEPData-ins1748602-v1-M60_bestCombined_UpperLimit.csv" ],
    units = [ None, 'pb' ], dataFormats = [ 'csv', 'csv'] )

T6bbHH130           = T6bbHH.addMassPlane( 2*[[ x, y, y-130. ]] )
T6bbHH130.figure    = "figaux_01b"
T6bbHH130.figureUrl = "https://atlas.web.cern.ch/Atlas/GROUPS/PHYSICS/PAPERS/SUSY-2018-31/figaux_01b.png"
T6bbHH130.dataUrl   = "https://www.hepdata.net/record/ins1748602?version=1&table=DM130_bestCombined_UpperLimit"
T6bbHH130.setSources(dataLabels = [ 'obsExclusion', 'upperLimits' ],
    dataFiles = [ "orig/HEPData-ins1748602-v1-DM130_Obs.csv", \
                  "orig/HEPData-ins1748602-v1-DM130_bestCombined_UpperLimit.csv" ],
    units = [ None, 'pb' ], dataFormats = [ 'csv', 'csv'] )

databaseCreator.create()
