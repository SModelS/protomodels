#!/usr/bin/env python3

""" create tex table for the particle contents of protomodel runs """

import glob, os, sys
from smodels_utils.helper.sparticleNames import SParticleNames
# from builder import protomodel

class Table:
    def __init__ ( self, pattern, realpattern ):
        """
        :param pattern: pattern of dict files to consider
        """
        self.pattern = pattern
        self.realpattern = realpattern
        self.namer = SParticleNames ( susy = False )

    def write ( self, outfile ):
        outh = open ( outfile, "wt" )
        files = glob.glob ( self.pattern )
        files += glob.glob ( self.realpattern )
        files.sort()
        for f in files:
            label = "signal"
            labels = [ "fake", "signal", "real" ]
            nr = os.path.basename ( f )
            nr = nr.replace(".dict","")
            for l in labels:
                nr = nr.replace(l,"")
                if l in f:
                    label = l
            nr = int(nr)
            h = open ( f, "rt" )
            txt = eval(h.read())
            h.close()
            K = txt[0]["K"]
            pids = list ( txt[0]["masses"].keys() )
            pids.remove ( 1000022 )
            def sorter ( x ):  ## we sort such that stop is first and the lbp is last
                ## which is not the real order but a close bet
                if x == 1000006:
                    x-=100000
                if x == 1000022:
                    x+= 1000000
                return x
            pids.sort( key= sorter )
            if True:
                label = label.replace("signal", "fake" )
            prtcles = self.namer.texName ( pids, addDollars=True )
            line = "%s \\#%d & %.3f & %s \\\\ \n" % \
                   ( label, nr, K, prtcles )  
            outh.write ( line )
        outh.close()


def main():
    import argparse
    argparser = argparse.ArgumentParser(description="tex table creator")
    argparser.add_argument ( '-d', '--dictfiles', nargs='*',
            help='input dictionary files [./signal*.dict]',
            type=str, default='./signal*.dict' )
    argparser.add_argument ( '-r', '--realfiles', nargs='*',
            help='real dictionary files [./real*.dict]',
            type=str, default='./real*.dict' )
    argparser.add_argument ( '-o', '--outfile', nargs='?',
            help='output file [./table.tex]',
            type=str, default='./table.tex' )
    args=argparser.parse_args()
    table = Table ( args.dictfiles, args.realfiles )
    table.write( args.outfile )

if __name__ == "__main__":
    main()
