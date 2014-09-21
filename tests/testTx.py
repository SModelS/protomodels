#!/usr/bin/env python

"""
.. module:: testTx
   :synopsis: Tests with Tx slha input files.
    
.. moduleauthor:: Wolfgang Waltenberger <wolfgang.waltenberger@gmail.com>
    
"""

import setPath
from smodels.theory import slhaDecomposer
from smodels.tools import xsecComputer
from smodels.tools.xsecComputer import NLL
from smodels.tools.physicsUnits import GeV, fb, TeV
import unittest
import logging

class TxTest(unittest.TestCase):
    logger = logging.getLogger(__name__)

    def testT1(self):
        self.logger.info ( "T1" )
        """ test with the T1 slha input file """
        slhafile="../inputFiles/slha/T1xsecs.slha"
        topos = slhaDecomposer.decompose ( slhafile, .1*fb, False, False, 5.*GeV )
        for topo in topos:
            for element in topo.elementList:
                masses=element.getMasses()
                # print "e=",element,"masses=",masses
                mgluino=masses[0][0]
                mLSP=masses[0][1]
                self.assertEqual ( str(element), "[[[jet,jet]],[[jet,jet]]]" )
                self.assertEqual ( int ( mgluino / GeV ), 675 )
                self.assertEqual ( int ( mLSP / GeV ), 600 )

if __name__ == "__main__":
    unittest.main()
