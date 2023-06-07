#!/usr/bin/env python3

""" restore hiscore.hi from the dictionary file """

# from ptools import hiscoreTools
from walker.hiscore import Hiscore
from builder.manipulator import Manipulator
from builder.protomodel import ProtoModel
from tester.predictor import Predictor
from smodels.tools.smodelsLogging import logger
logger.setLevel("ERROR")
import subprocess

## maybe add copying of real*.dict to hiscores.dict
## maybe add check for database pickle file
subprocess.getoutput ( "rm H*hi" )
subprocess.getoutput ( "rm Kold.conf" )

pmodel = ProtoModel ( 0 )
pr = Predictor ( 0, do_combine = False )
ma = Manipulator( pmodel )
ma.initFromDictFile ( "hiscores.dict", initTestStats=True )
print ( "The previous K value was", ma.M.K )
pr.predict ( ma.M )
print ( "We end up with K=", ma.M.K )
hi = Hiscore ( 0, True, picklefile = "H1.hi", hiscores = [ ma.M ] )
hi.writeListToPickle()

#a = subprocess.getoutput ( "./upHi.py" ) # no idea what that was
#print ( a )
