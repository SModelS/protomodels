#!/usr/bin/env python3

""" analyse the meta statistics of database.dict """

import matplotlib
matplotlib.use('agg')
from matplotlib import pyplot as plt
import numpy as np
import os, glob, pickle, sys
import scipy.stats
import matplotlib.mlab as mlab
from smodels_utils.helper.various import getSqrts, findCollaboration
from ptools.moreHelpers import namesForSetsOfTopologies
from typing import Union, Text, List

class Analyzer:
    def __init__ ( self, pathname : str, topos : Union[Text,None,List], 
                   nocolors : bool = False ):
        """
        :param pathname: filename of dictionary
        :param topos: topologies to filter for
        """
        self.setColors ( nocolors )
        self.reportZvalues = True
        self.filenames = []
        topos, _ = namesForSetsOfTopologies ( topos )
        
        self.topos = set()
        if topos != None:
            for k in topos.split(","):
                self.topos.add ( k.strip() )
        for pname in pathname:
            if os.path.isdir ( pname ):
                pname = pname + "/db*dict"
            self.filenames += glob.glob ( pname )
        self.pvalues = { 8:[], 13:[] }
        self.Zvalues = { 8:[], 13:[] }

    def setColors ( self, nocolors ):
        self.green,self.reset="",""
        if not nocolors:
            try:
                import colorama
                self.green = colorama.Fore.GREEN
                self.reset = colorama.Fore.RESET
            except ImportError as e:
                print ( "no colorama" )

    def summarize ( self ):
        if self.reportZvalues:
            Ztot = self.Zvalues[8]+self.Zvalues[13]
            self.pprint ( f"Zavg(total) =  {np.mean(Ztot):.2f}+-{np.std(Ztot)/np.sqrt(len(Ztot)):.3f}" )
            self.pprint ( f"Zavg( 8tev) = {np.mean(self.Zvalues[8]):.2f}+-{np.std(self.Zvalues[8])/np.sqrt(len(self.Zvalues[8])):.3f}" )
            self.pprint ( f"Zavg(13tev) =  {np.mean(self.Zvalues[13]):.2f}+-{np.std(self.Zvalues[13])/np.sqrt(len(self.Zvalues[13])):.3f}" )
        else:
            self.pprint ( f"pavg={np.mean(pavg):.2f}" )
            self.pprint ( f"pavg(13tev)={np.mean(pavg13):.2f}" )

    def analyze ( self, nlargest : int, nsmallest : int, enum : bool ):
        for filename in self.filenames:
            self.analyzeFile ( filename, nlargest, nsmallest, enum )

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

    def pprint ( self, *args ):
        print ( f"[analyseDBDict] {self.green}{' '.join(args)}{self.reset}" )

    def analyzeFile ( self, filename : str, nlargest : int, nsmallest : int,
           enum : bool ):
        self.pprint ( f"reading {filename}" )
        meta, data = self.read ( filename )
        byS, byp = {}, {}
        anas = set()
        for anaid, values in data.items():
            topos = self.getTopos ( values, anaid )
            if topos == None:
                continue
            anas.add ( anaid[:anaid.find(":")] )

            if "origS" in values:
                byS[ values["origS"] ] = ( anaid, values )
            if "orig_p" in values:
                byp[ values["orig_p"] ] = ( anaid, values )
        colls = [ findCollaboration(x) for x in anas ]
        self.pprint ( f"{colls.count('CMS')} CMS and {colls.count('ATLAS')} ATLAS results" )
        keys = list ( byp.keys() )
        keys.sort( reverse = False )
        #keys.sort( reverse = True )
        for ctr,k in enumerate(keys):
            ana = byp[k][0]
            sqrts = getSqrts ( ana )
            values = byp[k][1]
            p = values["orig_p"]
            Z = - scipy.stats.norm.ppf ( p )
            self.pvalues[sqrts].append(p)
            self.Zvalues[sqrts].append(Z)
        for ctr,k in enumerate(keys[:nlargest]):
            values = byp[k][1]
            ana = byp[k][0]
            topos = self.getTopos ( values, ana )
            # p = 1. - scipy.stats.norm.cdf(k)
            p = values["orig_p"]
            obsN = values["origN"]
            expBG = values["expectedBG"]
            bgErr = values["bgError"]
            Z = - scipy.stats.norm.ppf ( p )
            line = ""
            if enum:
                line += f"#{ctr+1:2d}: "
            if self.reportZvalues:
                line += f"Z={Z:.2f}:"
            else:
                line += f"p={k:.2f}:"
            line += f" {ana} {topos} (obsN={obsN:.0f}, bg={expBG:.2f}+-{bgErr:.2f})"
            print ( line )
        print ()
        for ctr,k in enumerate(keys[-nsmallest:]):
            values = byp[k][1]
            ana = byp[k][0]
            topos = self.getTopos ( values, ana )
            p = values["orig_p"]
            obsN = values["origN"]
            expBG = values["expectedBG"]
            bgErr = values["bgError"]
            Z = - scipy.stats.norm.ppf ( p )
            line = ""
            if enum:
                line += f"#{ctr+1:2d}: "
            if self.reportZvalues:
                line += f"Z={Z:.2f}:"
            else:
                line += f"p={k:.2f}:"
            line += f" {ana} {topos} (obsN={obsN:.0f}, bg={expBG:.2f}+-{bgErr:.2f})"
            print ( line )
        
        self.summarize()



def main():
    import argparse
    argparser = argparse.ArgumentParser(description="meta statistics analyzer")
    argparser.add_argument ( '-d', '--dictfile', nargs='*',
            help='input dictionary file(s) [../data/database/]',
            type=str, default='.,/data/database/' )
    argparser.add_argument ( '-t', '--topos', nargs='*',
            help='filter for topologies, comma separated list or multiple arguments [None]',
            type=str, default=None )
    argparser.add_argument ( '-n', '--nlargest',
            help='number of result to list with largest Z values [10]',
            type=int, default=10 )
    argparser.add_argument ( '-N', '--nsmallest',
            help='number of result to list with smallest Z values [3]',
            type=int, default=3 )
    argparser.add_argument ( '-e', '--enumerate',
            help='enumerate the list', action="store_true" )
    argparser.add_argument ( '--nocolors',
            help='dont use colors in output', action="store_true" )
    args=argparser.parse_args()
    analyzer = Analyzer ( args.dictfile, args.topos, args.nocolors )
    analyzer.analyze ( args.nlargest, args.nsmallest, args.enumerate )

if __name__ == "__main__":
    main()
