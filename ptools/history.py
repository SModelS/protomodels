#!/usr/bin/env python3

<<<<<<< HEAD
""" a class meant for recording the history of a walk.
=======
""" a class meant for recording the history of a walk. 
>>>>>>> 96cd173644c9e15dfbe797f626e22f1adf1ddb87
When executed, it runs a typical walker with the history recorder. """

import colorama, time

class History:
    def __init__ ( self, outfile = "history.list" ):
        """
        :param outfile: the file to story history in
        """
        self.outfile = outfile
        self.handle = open ( outfile, "wt" )
        self.handle.write ( f"# history recording {time.asctime()}\n" )
        self.handle.write ( "[" )
<<<<<<< HEAD
        self.pprint ( f"starting history recorder at {self.outfile}." )

    def add ( self, ma ):
        """ add protomodel from step
        :param ma: manipulator object
        """
        protomodel = ma.M
=======
        self.pprint ( f"staring history recorder at {self.outfile}." )

    def add ( self, protomodel ):
        """ add protomodel from step """
>>>>>>> 96cd173644c9e15dfbe797f626e22f1adf1ddb87
        D = protomodel.dict()
        D["K"]=protomodel.K
        D["Z"]=protomodel.Z
        D["step"]=protomodel.step
<<<<<<< HEAD
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
=======
>>>>>>> 96cd173644c9e15dfbe797f626e22f1adf1ddb87
        self.pprint ( f"adding protomodel to {self.outfile}." )
        self.handle.write ( str(D)+",\n" )
        self.handle.flush()

    def pprint ( self, *args ):
        """ pprint """
        print ( "%s[history] %s%s" % \
                ( colorama.Fore.LIGHTBLUE_EX, " ".join(map(str,args)), colorama.Fore.RESET ) )

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
    import csetup
    csetup.setup()
    from walker.randomWalker import RandomWalker
<<<<<<< HEAD
    walker = RandomWalker ( 0, 2000, "aggressive", False, 0, args.database, record_history=True )
=======
    hi = History()
    walker = RandomWalker ( 0, 2000, "aggressive", False, 0, args.database )
    walker.recorder = hi
>>>>>>> 96cd173644c9e15dfbe797f626e22f1adf1ddb87
    walker.walk()
    hi.save()
