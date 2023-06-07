#!/usr/bin/env python3

"""
.. module:: testPredictor
   :synopsis: Test computation of model predictions

.. moduleauthor:: Andre Lessar lessa.a.p@gmail.com>

"""

import sys,os,glob
sys.path.insert(0,"../")
try:
    import smodels
except:
    from ptools import setPath
import unittest
import pickle
from tester.predictor import Predictor
import numpy as np


class PredictionsTest(unittest.TestCase):

    def testPredictions(self):

        with open('randomModels_default.pcl','rb') as f:
            protomodel = pickle.load(f)[0]

        protomodel.dbversion = None
        protomodel.codeversion = None
        pNew = protomodel.copy()
        pNew.templateSLHA = os.path.abspath('../builder/templates/template1g.slha')
        pNew.Z = None
        pNew.llhd = None
        pNew.bestCombo = None
        pNew.muhat = None
        pNew.mumax = None
        pNew.tpList = []


        #Test against old version:
        predictor =  Predictor( 0, dbpath='./database.pcl',
                              expected=False, select='all', do_combine=False )
        predictor.rthreshold = 1.7
        predictor.predict(pNew)
        #OBS: Since the original protomodel already has all of its cross-sections rescaled, we do not
        #need to rescale pNew again (since muhat should be 1)

        #Hard coded combarison against (old) version (minor xsec differences and the change in prior):
        oldPrior = predictor.combiner.computePrior(pNew,name="expo2",nll=True)
        newPrior = predictor.combiner.computePrior(pNew,name="expo1",nll=True)
        newK = pNew.K - 2*oldPrior + 2*newPrior
        self.assertTrue(abs(3.9-newK)/3.9 < 0.1)
        self.assertTrue(abs(2.3-pNew.Z)/2.3 < 0.1)
        self.assertEqual(15,len(pNew.rvalues))
        self.assertEqual('BDGK',pNew.letters)
        #From previous (pre-refac) version:
        rvalues = [1.725, 1.216, 1.133,0.824,0.536, 0.285, 0.27 , 0.198, 0.147,0.112, 0.074, 0.073, 0.052, 0.049, 0.00242]
        for i,r in enumerate(rvalues):
            self.assertTrue(abs(pNew.rvalues[i]-r)/r < 0.05) #Allow no more than 5% differences (due to xsec)


        #Compare against new default:
        predictor.rthreshold = 1.3
        predictor.predict(pNew)
        self.assertAlmostEqual(protomodel.Z,pNew.Z,3)
        self.assertAlmostEqual(protomodel.K,pNew.K,3)
        self.assertAlmostEqual(protomodel.muhat,pNew.muhat,3)
        self.assertAlmostEqual(protomodel.mumax,pNew.mumax,3)
        self.assertAlmostEqual(protomodel.llhd,pNew.llhd,3)


        np.testing.assert_almost_equal(protomodel.rvalues,pNew.rvalues,3)
        self.assertEqual(len(protomodel.bestCombo),len(pNew.bestCombo))
        for i,pred in enumerate(protomodel.bestCombo):
            self.assertEqual(str(pred.expResult),str(pNew.bestCombo[i].expResult))
            self.assertEqual(str(pred.dataset),str(pNew.bestCombo[i].dataset))
            self.assertAlmostEqual(pred.xsection.value.asNumber(),pNew.bestCombo[i].xsection.value.asNumber(),3)
            self.assertAlmostEqual(pred.upperLimit.asNumber(),pNew.bestCombo[i].upperLimit.asNumber(),3)
        for i,pred in enumerate(protomodel.tpList):
            self.assertAlmostEqual(pred[0],pNew.tpList[i][0])
            self.assertEqual(str(pred[2].expResult),str(pNew.tpList[i][2].expResult))
            self.assertAlmostEqual(pred[2].xsection.value.asNumber(),pNew.tpList[i][2].xsection.value.asNumber(),3)


        #Remove files generated during run
        for f in glob.glob('.cur*slha'):
            os.remove(f)



if __name__ == "__main__":
    unittest.main()
