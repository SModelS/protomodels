#!/usr/bin/env python3

"""
.. module:: testVertically
   :synopsis: Testing "vertically", meaning we run a walker with a defined database and
              seed.

.. moduleauthor:: Wolfgang Waltenberger <wolfgang.waltenberger@gmail.com>

"""

import sys,os
sys.path.insert(0,"../")
import unittest
from walker.randomWalker import RandomWalker
from ptools import helpers

class VerticalTest(unittest.TestCase):

    def testRun(self):

        if os.path.isfile('H0.hi'):
            os.remove('H0.hi')
        helpers.seedRandomNumbers ( 1 )
        walker = RandomWalker ( nsteps=11, dbpath="./database.pcl", nevents = 10000 )
        walker.predictor.rthreshold = 1.3 #Make sure to use the correct threshold
        walker.walk()

        self.assertAlmostEqual ( walker.protomodel.K, 2.396, 2 )
        self.assertEqual ( len(walker.hiscoreList.hiscores), 1 )
        self.assertAlmostEqual ( walker.hiscoreList.hiscores[0].Z, 2.22952,2 )
        self.assertEqual ( walker.hiscoreList.hiscores[0].step, 4 )
        self.assertAlmostEqual ( walker.protomodel.masses[1000001], 791.1794, 2 )
        self.assertAlmostEqual ( walker.protomodel.masses[1000022], 240.3092, 2 )
        self.assertAlmostEqual ( walker.hiscoreList.hiscores[0].masses[1000001], 791.1794, 2 )
        self.assertAlmostEqual ( walker.hiscoreList.hiscores[0].masses[1000022], 240.3092, 2 )


if __name__ == "__main__":
    unittest.main()
