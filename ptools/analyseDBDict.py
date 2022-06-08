#!/usr/bin/env python3

""" analyse the meta statistics of database.dict """

import matplotlib
matplotlib.use('agg')
from matplotlib import pyplot as plt
import numpy as np
import os, glob, pickle, sys
import scipy.stats
import matplotlib.mlab as mlab
from smodels_utils.helper.various import getSqrts

class Analyzer:
    def __init__ ( self, pathname, topos ):
        """
        :param pathname: filename of dictionary
        :param topos: topologies to filter for
        """
        self.filenames = []
        self.topos = set()
        if topos not in [ None, "", [] ]:
            for to in topos:
                for k in to.split(","):
                    self.topos.add ( k.strip() )
        for pname in pathname:
            if os.path.isdir ( pname ):
                pname = pname + "/db*dict"
            self.filenames += glob.glob ( pname )

    def analyze ( self, nsmallest, nlargest ):
        for filename in self.filenames:
            self.analyzeFile ( filename, nsmallest, nlargest )

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

    def getTopos ( self, values, ana ):
        # we filter with self.topos
        if "txns" in values:
            ret = values["txns"]
            tret = ret.split(",")
            isIn = False
            if len(self.topos)==0:
                isIn = True
            else:
                for t in tret:
                    if t in self.topos:
                        isIn = True
            if not isIn:
                return None
            if len(ret)>15:
                for i in range(15,5,-1):
                    if ret[i]==",":
                        break
                ret=ret[:i+1]+" ..."
            if ret == "":
                print ( f"empty txns in {ana}? >>{values['txns']}<<" )
            return ret
        topos=[]
        for vk,vv in values.items():
            if vk.startswith("sigNT"):
                p = vk.find("T")
                topo = vk[p:]
                topos.append ( topo )
        return ",".join(topos)

    def analyzeFile ( self, filename, nsmallest, nlargest ):
        print ( f"[analyzeDBDict] {filename}" )
        meta, data = self.read ( filename )
        byS, byp = {}, {}
        for anaid, values in data.items():
            topos = self.getTopos ( values, anaid )
            if topos == None:
                continue
            if "origS" in values:
                byS[ values["origS"] ] = ( anaid, values )
            if "orig_p" in values:
                byp[ values["orig_p"] ] = ( anaid, values )
        keys = list ( byp.keys() )
        keys.sort( reverse = False )
        #keys.sort( reverse = True )
        pavg, pavg13 = [], []
        for ctr,k in enumerate(keys):
            ana = byp[k][0]
            sqrts = getSqrts ( ana )
            values = byp[k][1]
            p = values["orig_p"]
            pavg.append ( p )
            if sqrts > 10:
                pavg13.append ( p )
        for ctr,k in enumerate(keys[:nsmallest]):
            values = byp[k][1]
            ana = byp[k][0]
            topos = self.getTopos ( values, ana )
            # p = 1. - scipy.stats.norm.cdf(k)
            p = values["orig_p"]
            obsN = values["origN"]
            expBG = values["expectedBG"]
            bgErr = values["bgError"]
            print( "p=%.2f: %s %s (obsN=%d, bg=%.2f+-%.2f)" % ( k, ana, topos, obsN, expBG, bgErr ) )
        print ()
        for ctr,k in enumerate(keys[-nlargest:]):
            values = byp[k][1]
            ana = byp[k][0]
            topos = self.getTopos ( values, ana )
            obsN = values["origN"]
            expBG = values["expectedBG"]
            bgErr = values["bgError"]
            print( "p=%.2f: %s %s (obsN=%d, bg=%.2f+-%.2f)" % ( k, ana, topos, obsN, expBG, bgErr ) )
        
        print ( f"[analyzeDBDict] pavg={np.mean(pavg):.2f}" )
        print ( f"[analyzeDBDict] pavg(13tev)={np.mean(pavg13):.2f}" )


def main():
    import argparse
    argparser = argparse.ArgumentParser(description="meta statistics analyzer")
    argparser.add_argument ( '-d', '--dictfile', nargs='*',
            help='input dictionary file(s) [../data/database/]',
            type=str, default='.,/data/database/' )
    argparser.add_argument ( '-t', '--topos', nargs='*',
            help='filter for topologies, comma separated list or multiple arguments [None]',
            type=str, default=None )
    argparser.add_argument ( '-n', '--nsmallest',
            help='number of result to list with small p values [10]',
            type=int, default=10 )
    argparser.add_argument ( '-N', '--nlargest',
            help='number of result to list with large p values [3]',
            type=int, default=3 )
    args=argparser.parse_args()
    analyzer = Analyzer ( args.dictfile, args.topos )
    analyzer.analyze ( args.nsmallest, args.nlargest )

if __name__ == "__main__":
    main()
