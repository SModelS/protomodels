#!/usr/bin/env python3

import IPython

from builder.manipulator import Manipulator
from tester.predictor import Predictor

m = Manipulator( "hiscores.dict" )
predictor = Predictor(0, do_combine=True )
predictor.predict ( m.M, keep_predictions=True )

IPython.embed ( colors = "neutral" )
