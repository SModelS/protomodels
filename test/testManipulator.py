#!/usr/bin/env python3

"""
.. module:: testManipulator
   :synopsis: Test random changes in model

.. moduleauthor:: Andre Lessar lessa.a.p@gmail.com>

"""

import sys,os
sys.path.insert(0,"../")
try:
    import smodels
except:
    from tools import setPath
import unittest
import pickle,copy
from tools import helpers
from builder.manipulator import Manipulator

class PredictionsTest(unittest.TestCase):

    def testRandomModel(self):

        helpers.seedRandomNumbers ( 456 )
        with open('randomModels_default.pcl','rb') as f:
            pList= pickle.load(f) #List with original models and modified ones

        pNew = copy.deepcopy(pList[0]) #Obs: Can not use ProtoModel.copy(), since it calls random
        pNew.templateSLHA = os.path.abspath('../builder/templates/template1g.slha')
        m = Manipulator(pNew)

        for p in pList:
            self.assertEqual(pNew.masses,p.masses)
            self.assertEqual(pNew.decays,p.decays)
            self.assertEqual(pNew.ssmultipliers,p.ssmultipliers)
            self.assertEqual(pNew._stored_xsecs,p._stored_xsecs)
            m.randomlyChangeModel()


    def testRandomUnfreeze(self):
        helpers.seedRandomNumbers ( 456 )
        with open('randomModels_default.pcl','rb') as f:
            model= pickle.load(f)[0]

        m = Manipulator(model)
        m.randomlyUnfreezeParticle(sigma=0.5)
        self.assertEqual(sorted(model.unFrozenParticles()), [1000001, 1000021, 1000022])
        masses = {1000022: 215.71, 1000001: 941.23, 1000021: 1068.91}
        for pid,mass in masses.items():
            self.assertAlmostEqual(model.masses[pid],mass,places=1)

        decays = {1000022: {},
                1000001: {1000022: 1.0},
                1000021: {(1000022, 6): 0.0, (1000022, 5): 0.5035,
                (1000022, 21): 0.248, (1000022, 1): 0.202,
                (1000001, 1): 0.0, (1000022,4): 0.0456}
                }
        for pid,decay in decays.items():
            for dpid,br in decay.items():
                self.assertAlmostEqual(model.decays[pid][dpid],br,places=1)

        m.randomlyUnfreezeParticle(force= True)
        self.assertEqual(sorted(model.unFrozenParticles()), [1000001, 1000021, 1000022, 1000024])

    def testUnfreeze(self):
        helpers.seedRandomNumbers ( 123 )
        with open('randomModels_default.pcl','rb') as f:
            modelList = pickle.load(f)

        m = Manipulator(modelList[4])
        r = m.unFreezeParticle(pid = 2000006) #Should not work since 1000006 is frozen!
        self.assertEqual(r,0)
        self.assertEqual(sorted(m.M.unFrozenParticles()), [1000022, 1000024])
        r = m.unFreezeParticle(pid = 1000006) #Should not work since 1000006 is frozen!
        r = m.unFreezeParticle(pid = 2000006) #Now it should work since 1000006 is unfrozen!
        self.assertEqual(r,1)
        self.assertEqual(sorted(m.M.unFrozenParticles()), [1000006, 1000022, 1000024, 2000006])
        self.assertTrue(m.M.masses[2000006] > m.M.masses[1000006]) #stop2 mass should be higher than stop1
        minMass = min([mass for mass in m.M.masses.values()])
        self.assertTrue(minMass == m.M.masses[m.M.LSP]) #all masses should smaller than the LSP mass

        #Check if the only open decays are created:
        m = Manipulator(modelList[1])
        r = m.unFreezeParticle(pid = 1000021, force=True)
        self.assertEqual(r,1)
        self.assertEqual(sorted(m.M.unFrozenParticles()), [1000001, 1000021, 1000022])
        self.assertEqual(sorted(m.M.decays.keys()), [1000001, 1000021, 1000022])
        self.assertEqual(sorted(m.M.decays[1000021].keys()),[(1000001, 1), (1000022, 1), (1000022, 4), (1000022, 5),
                                                  (1000022, 6), (1000022, 21)])
        for pid in m.M.decays[1000021]:
            self.assertTrue(m.M.masses[1000021] > m.M.masses[max(pid)])



if __name__ == "__main__":
    unittest.main()
