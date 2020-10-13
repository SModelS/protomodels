#!/usr/bin/env python3

"""
.. module:: testServer
   :synopsis: Tests the database server

.. moduleauthor:: Wolfgang Waltenberger <wolfgang.waltenberger@gmail.com>

"""

import sys,os,glob
sys.path.insert(0,"../")
import unittest
import glob
from smodels.tools import crashReport
from smodels.tools.timeOut import NoTime
from smodels.experiment.databaseObj import Database
from unitTestHelpers import equalObjs, runMain, importModule
import time, random
import subprocess

from smodels.tools.smodelsLogging import logger, setLogLevel
from smodels.tools.databaseClient import DatabaseClient
from smodels.tools.databaseServer import DatabaseServer
import warnings

setLogLevel ( "warn" )

class ServerTest(unittest.TestCase):
    def mestRun(self):
        # port = random.choice ( range(31700, 42000 ) )
        port = 31771
        # pclfile = glob.glob("database/db*.pcl" )[0]
        pclfile = "database/db30.pcl"
        #server.run ( nonblocking = True )
        startserver = f"../smodels/tools/smodelsTools.py proxydb -p {port} -i {pclfile} -o ./proxy.pcl -r"
        cmd = startserver.split(" ")
        print ( "starting server %s" % startserver )
        popen = subprocess.Popen ( cmd, stdout = subprocess.PIPE )
        print ( "started server" )
        time.sleep(3)
        filename = "./testFiles/slha/simplyGluino.slha"
        db = Database("./proxy.pcl" )
        outputfile = runMain(filename,suppressStdout = False, overridedatabase = db )
        # for stdout_line in iter(popen.stdout.readline, ""):
        #    print ( stdout_line )
        popen.stdout.close()
        client = DatabaseClient ( port = port ) 
        client.send_shutdown()

    def testCompare(self):
        warnings.simplefilter("ignore", ResourceWarning)
        from simplyGluino_default import smodelsOutputDefault
        filename = "./testFiles/slha/simplyGluino.slha"
        port = random.choice ( range(31700, 42000 ) )
        # port = 31744

        startserver = f"../smodels/tools/smodelsTools.py proxydb -p {port} -i database/db30.pcl -o ./proxy.pcl -r -v error"
        cmd = startserver.split(" ")
        # print ( "starting server %s" % startserver )
        subprocess.Popen ( cmd )

        time.sleep ( 2 )

        db = Database("./proxy.pcl" )
        outputfile = runMain(filename,suppressStdout = True, overridedatabase = db )
        smodelsOutput = importModule ( outputfile )
        
        client = DatabaseClient ( port = port, verbose = "warn" ) 
        client.send_shutdown()

        ignoreFields = ['input file','smodels version', 'ncpus', 'Element', 
                        'database version', 'Total missed xsec',
                        'Missed xsec long-lived', 'Missed xsec displaced', 
                        'Missed xsec MET', 'Total outside grid xsec',
                        'Total xsec for missing topologies (fb)',
                        'Total xsec for missing topologies with displaced decays (fb)',
                        'Total xsec for missing topologies with prompt decays (fb)',
                        'Total xsec for topologies outside the grid (fb)' ]
        smodelsOutputDefault['ExptRes'] = sorted(smodelsOutputDefault['ExptRes'],
                    key=lambda res: res['r'], reverse=True)
        equals = equalObjs(smodelsOutput,smodelsOutputDefault,allowedDiff=0.08,
                           ignore=ignoreFields, fname = outputfile )
        if not equals:
            e =  "simplyGluino.py != simplyGluino_default.py"
            logger.error( e )
            # raise AssertionError( e )

        self.assertTrue(equals)
        self.removeOutputs( outputfile )

    def removeOutputs( self, f ):
        """ remove cruft outputfiles """
        for i in [ f, f.replace(".py",".pyc") ]:
            if os.path.exists( i ): os.remove( i )

if __name__ == "__main__":
    unittest.main()
