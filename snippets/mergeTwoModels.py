#!/usr/bin/env python

import time
from ptools.hiscoreTools import mergeTwoModels
model = mergeTwoModels ( "pmodel_ewkino.dict", "pmodel_hadrons.dict" )
print ( "model", model )
f=open("merged.py","wt")
f.write ( "# this model is due to a merge of the ewkino and the hadrons models\n" )
f.write ( f"# created {time.asctime()}\n" )
f.write ( str(model)+"\n" )
f.close()
