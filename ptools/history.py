#!/usr/bin/env python3

""" a class meant for recording the history of a walk. 
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
        self.pprint ( f"staring history recorder at {self.outfile}." )

    def add ( self, protomodel ):
        """ add protomodel from step """
        D = protomodel.dict()
        D["K"]=protomodel.K
        D["Z"]=protomodel.Z
        D["step"]=protomodel.step
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
    hi = History()
    walker = RandomWalker ( 0, 2000, "aggressive", False, 0, args.database )
    walker.recorder = hi
    walker.walk()
    hi.save()
