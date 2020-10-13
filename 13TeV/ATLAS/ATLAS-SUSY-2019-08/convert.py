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
info                 = MetaInfoInput('ATLAS-SUSY-2019-08')
info.url             = 'https://atlas.web.cern.ch/Atlas/GROUPS/PHYSICS/PAPERS/SUSY-2019-08/'
info.sqrts             = 13
info.lumi             = 139.
info.prettyName     = '1L + higgs + Etmiss (EWino)'
info.private         = False
info.arxiv             = 'arXiv:1909.09226'
info.contact         = 'atlas-phys-susy-conveners@cern.ch'
info.publication     = ''

# +++++++ dataset block ++++++++++++++
dataset = DataSetInput('data')
dataset.setInfo( dataType = 'upperLimit', dataId = None )

# +++++++ next txName block ++++++++++++++
TChiWH = dataset.addTxName("TChiWH" )
TChiWH.checked              = 'no'
TChiWH.constraint           = "[[['W']],[['higgs']]]"
TChiWH.conditionDescription = None
TChiWH.condition            = None
TChiWH.source               = 'ATLAS'

# +++++++ next mass plane block ++++++++++++++
TChiWH1           = TChiWH.addMassPlane( 2*[[ x, y ]] )
TChiWH1.figure    = "fig_06"
TChiWH1.figureUrl = "https://atlas.web.cern.ch/Atlas/GROUPS/PHYSICS/PAPERS/SUSY-2019-08/figaux_03.png"
TChiWH1.dataUrl   = "https://doi.org/10.17182/hepdata.90607.v1/t17"

TChiWH1.setSources(dataLabels = [ 'obsExclusion', 'obsExclusionP1', 'obsExclusionM1', 'expExclusion', 'upperLimits' ],
    dataFiles = [ "orig/HEPData-ins1755298-v1-Observed_limit_1lbb.csv",
                  "orig/HEPData-ins1755298-v1-Observed_limit_1lbb_(Up).csv",
                  "orig/HEPData-ins1755298-v1-Observed_limit_1lbb_(Down).csv",
                  "orig/HEPData-ins1755298-v1-Expected_limit_1lbb.csv",
                  "orig/HEPData-ins1755298-v1-Upper_limits_1Lbb.csv" ],
    units = [ None, None, None, None, 'pb' ], dataFormats = [ 'csv' ]*5 )

databaseCreator.create()
