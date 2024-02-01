#!/usr/bin/env python3

""" a class meant for recording the history of a walk.
When executed, it runs a typical walker with the history recorder. """

import time
from colorama import Fore

class History:
    def __init__ ( self, outfile = "history.list" ):
        """
        :param outfile: the file to story history in
        """
        self.outfile = outfile
        self.handle = open ( outfile, "wt" )
        self.handle.write ( f"# history recording {time.asctime()}\n" )
        self.handle.write ( "[" )
        self.pprint ( f"starting history recorder at {self.outfile}." )

    def add ( self, ma ):
        """ add protomodel from step
        :param ma: manipulator object
        """
        protomodel = ma.M
        D = protomodel.dict()
        D["K"]=protomodel.K
        D["Z"]=protomodel.Z
        D["step"]=protomodel.step
        bc=[]
        if protomodel.bestCombo != None:
            for tp in protomodel.bestCombo:
                dI = tp.dataId()
                if dI == None:
                    dI="ul"
                bc.append ( tp.analysisId()+":"+dI )
        D["bestCombo"]=bc
        if len(ma.recording)>0:
            D["actions"]=ma.recording
            ma.recording=[]
        self.pprint ( f"adding protomodel to {self.outfile}." )
        self.handle.write ( str(D)+",\n" )
        self.handle.flush()

    def pprint ( self, *args ):
        """ pprint """
        print ( f"{Fore.LIGHTBLUE_EX}[history] {' '.join(map(str,args))}{Fore.RESET}" )

    def save ( self ):
        self.pprint ( f"writing history to {self.outfile}." )
        self.handle.write ( "]\n" )
        self.handle.close()

if __name__ == "__main__":
    import argparse
    from argparse import RawTextHelpFormatter
    argparser = argparse.ArgumentParser( description='history recorder' )
    argparser.add_argument ( '-d', '--database',
            help='database to use [./default.pcl]',
            type=str, default="./default.pcl" )
    argparser.add_argument ( '-o', '--outfile',
            help='file to write out recoreded history. ["history.list"]',
            type=str, default="history.list" )
    args = argparser.parse_args()
    from protomodels.csetup import setup
    setup()
    from walker.randomWalker import RandomWalker
    walker = RandomWalker ( 0, 2000, "aggressive", False, 0, args.database, record_history=True )
    walker.walk()
    hi.writeListToPickle()
