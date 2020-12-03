#!/usr/bin/env python3

""" analyse the meta statistics of database.dict """

import matplotlib
matplotlib.use('agg')
from matplotlib import pyplot as plt
import numpy as np
import os, glob, pickle, sys
import scipy.stats
import matplotlib.mlab as mlab

class Analyzer:
    def __init__ ( self, pathname ):
        """
        :param pathname: filename of dictionary
        :param reset: if true, then dont recycle pickle files
        """
        self.filenames = []
        for pname in pathname:
            if os.path.isdir ( pname ):
                pname = pname + "/db*dict"
            self.filenames += glob.glob ( pname )

    def analyze ( self ):
        for filename in self.filenames:
            self.analyzeFile ( filename )

    def read ( self, fname ):
        """ read in content of filename """
        with open( fname,"rt") as f:
            tmp=f.readlines()
        lines = []
        for line in tmp:
            if line.startswith("#"):
                continue
            lines.append ( line )
        basename = os.path.basename ( fname ).replace(".dict","")
        meta = eval(lines[0])
        nan=float("nan")
        data = eval("\n".join(lines[1:]))
        newdata = {}
        for i,v in data.items():
            if "expectedBG" in v and v["expectedBG"]>=0.:
                newdata[i]=v
            else:
                if i.endswith ( ":ul" ):
                    pass
                    #print ( f"[plotDBDict] removing {basename}:{i} (is an UL)" )
                else:
                    eBG=None
                    if "expectedBG" in v:
                        eBG = v["expectedBG"]
                    #print ( f"[plotDBDict] removing {basename}:{i} (eBG is {eBG})" )
        # print ( f"[plotDBDict] keeping {len(newdata)}/{len(data)} for {basename}" )
        return meta,newdata

    def getTopos ( self, values ):
        topos=[]
        for vk,vv in values.items():
            if vk.startswith("sigNT"):
                p = vk.find("T")
                topo = vk[p:]
                topos.append ( topo )
        return topos

    def analyzeFile ( self, filename ):
        print ( f"analyzing {filename}" )
        meta, data = self.read ( filename )
        byS = {}
        for anaid, values in data.items():
            if "S" in values:
                byS[ values["S"] ] = ( anaid, values )
        keys = list ( byS.keys() )
        keys.sort( reverse = True )
        for ctr,k in enumerate(keys[:5]):
            values = byS[k][1]
            topos = self.getTopos ( values )
            print( "S=%.2f: %s; %s" % ( k, byS[k][0], ",".join(topos) ) )
        print ()
        for ctr,k in enumerate(keys[-3:]):
            values = byS[k][1]
            topos = self.getTopos ( values )
            print( "S=%.2f: %s; %s" % ( k, byS[k][0], ",".join(topos) ) )
        


def main():
    import argparse
    argparser = argparse.ArgumentParser(description="meta statistics analyzer")
    argparser.add_argument ( '-d', '--dictfile', nargs='*',
            help='input dictionary file(s) [../data/database/]',
            type=str, default='.,/data/database/' )
    args=argparser.parse_args()
    analyzer = Analyzer ( args.dictfile )
    analyzer.analyze ( )

if __name__ == "__main__":
    main()
