#!/usr/bin/env python3

""" analyse the meta statistics of database.dict """

import matplotlib
matplotlib.use('agg')
from matplotlib import pyplot as plt
import numpy as np
import os, glob, pickle, sys
sys.path.insert(0,"../")
import scipy.stats
import matplotlib.mlab as mlab
from smodels_utils.helper.various import getSqrts, findCollaboration
from ptools.moreHelpers import namesForSetsOfTopologies
from typing import Union, Text, List, Dict
from protomodels.builder.loggerbase import LoggerBase
import subprocess, shutil

class Analyzer ( LoggerBase ):
    def __init__ ( self, pathname : str, topos : Union[Text,None,List], 
                   nocolors : bool = False, latex : bool = False,
                   enum : bool = False ):
        """
        :param pathname: filename of dictionary
        :param topos: topologies to filter for
        """
        super ( Analyzer, self ).__init__ ( 0 )
        self.setColors ( nocolors )
        self.latex = latex
        self.enum = enum
        self.latexHeader()
        self.reportZvalues = True
        self.filenames = []
        topos, _ = namesForSetsOfTopologies ( topos )
        if topos == "all":
            topos = None
        
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

    def latexHeader ( self ):
        if not self.latex:
            return
        self.latexfile = open ( "regions.tex", "wt" )
        self.latexfile.write ( f"% list of most significant excesses\n" )
        import time
        self.latexfile.write ( f"% file created {time.asctime()}\n" )
        self.latexfile.write ( f"\n" )
        lformat = "{|l|l|l|l|r|r|}"
        if self.enum:
            lformat = lformat[:2]+"l|"+lformat[3:]
        # self.latexfile.write ( r"\resizebox{\textwidth}{!}{" )
        self.latexfile.write ( r"\begin{tabular}"+lformat )
        self.latexfile.write ( "\n" )
        self.latexfile.write ( r"\hline" )
        if self.enum:
            self.latexfile.write ( r"{\bf nr} & " )
        self.latexfile.write ( r"{\bf Z} & {\bf analysis} & {\bf SR} & {\bf topo} & {\bf obsN} & {\bf expected} \\" )
        self.latexfile.write ( "\n" )
        self.latexfile.write ( r"\hline" )
        self.latexfile.write ( "\n" )

    def latexFooter( self ):
        if not self.latex:
            return
        self.latexfile.write ( r"\hline" )
        self.latexfile.write ( r"\end{tabular}" )
        self.latexfile.write ( "\n" )
        self.latexfile.close()
        if not os.path.exists ( "region_list.tex" ):
            cmd = "ln -s share/ExcessesListTemplate.tex ./region_list.tex"
            subprocess.getoutput ( cmd )
        cmd = "pdflatex region_list.tex"
        subprocess.getoutput ( cmd )
        if shutil.which ( "timg" ) != None:
            cmd  = "timg region_list.pdf"
            o = subprocess.getoutput ( cmd )
            print ( o )

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
            if len ( self.Zvalues[8] ) > 0:
                self.pprint ( f"Zavg( 8tev) = {np.mean(self.Zvalues[8]):.2f}+-{np.std(self.Zvalues[8])/np.sqrt(len(self.Zvalues[8])):.3f}" )
            if len ( self.Zvalues[13] ) > 0:
                self.pprint ( f"Zavg(13tev) =  {np.mean(self.Zvalues[13]):.2f}+-{np.std(self.Zvalues[13])/np.sqrt(len(self.Zvalues[13])):.3f}" )
        else:
            self.pprint ( f"pavg={np.mean(pavg):.2f}" )
            self.pprint ( f"pavg(13tev)={np.mean(pavg13):.2f}" )

    def analyze ( self, nlargest : int, nsmallest : int ):
        for filename in self.filenames:
            self.analyzeFile ( filename, nlargest, nsmallest )

    def read ( self, fname ):
        from ptools.expResModifier import readDictFile 
        d = readDictFile ( fname )
        return ( d["meta"], d["data"] )

    def topoIsIn ( self, topo : str ) -> Union[None,bool]:
        """ check if the given topo is compatible with self.topos.

        :param topo: e.g. TRV1
        :returns: true if topo is in selection
        """
        if len(self.topos)==0:
            return True
        hasOnlyNegativeTopos = True
        ntopo = f"^{topo}"
        for t in self.topos:
            if not t.startswith ( "^" ):
                hasOnlyNegativeTopos = False
            if topo == t:
                return True
            if ntopo == t: # we explictly veto these
                return False 
        if hasOnlyNegativeTopos:
            return True
        return None

    def getTopos ( self, values : Dict, ana : str ) -> str:
        """ get the topologies. """
        # we filter with self.topos
        if "txns" in values:
            ret = values["txns"]
            tret = ret.split(",")
            isIn = False
            for t in tret:
                if self.topoIsIn ( t ) == True:
                    isIn = True
                if self.topoIsIn ( t ) == False: # there was a veto!
                    isIn = False
                    break
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

    def writeLatex ( self, line ):
        if not self.latex:
            return
        self.latexfile.write ( line )

    def analyzeFile ( self, filename : str, nlargest : int, nsmallest : int ):
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
        if len(colls)==0:
            return
        keys = list ( byp.keys() )
        keys.sort( reverse = False )
        #keys.sort( reverse = True )
        for ctr,k in enumerate(keys):
            ana = byp[k][0]
            sqrts = getSqrts ( ana )
            values = byp[k][1]
            p = values["orig_p"]
            # Z = - scipy.stats.norm.ppf ( p )
            Z = values["orig_Z"]
            # print ( "@@2 Z", Z, values["orig_Z"] )
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
            if self.enum:
                line += f"#{ctr+1:2d}: "
                self.writeLatex ( f"{ctr+1:2d} & " )
            if self.reportZvalues:
                line += f"Z={Z:.2f}:"
                self.writeLatex ( f"{Z:.2f} & " )
            else:
                line += f"p={k:.2f}:"
                self.writeLatex ( f"{p:.2f} & " )
            line += f" {ana} {topos} (obsN={obsN:.0f}, bg={expBG:.2f}+-{bgErr:.2f})"
            anaonly = ana[:ana.find(":")]
            sr = ana[ana.find(":")+1:].replace("_",r"\_")
            self.writeLatex ( f"{anaonly} & {sr} & {topos} & {obsN:.0f} & {expBG:.2f}$\\pm${bgErr:.2f} \\\\\n" )
            print ( line )
        print ()
        if nsmallest>0:
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
                if self.enum:
                    line += f"#{ctr+1:2d}: "
                    self.writeLatex ( f"{ctr+1:2d} & " )
                if self.reportZvalues:
                    line += f"Z={Z:.2f}:"
                    self.writeLatex ( f"{Z:.2f} & " )
                else:
                    line += f"p={k:.2f}:"
                line += f" {ana} {topos} (obsN={obsN:.0f}, bg={expBG:.2f}+-{bgErr:.2f})"
                self.writeLatex ( f"{anaonly} & {sr} & {topos} & {obsN:.0f} & {expBG:.2f}+-{bgErr:.2f} \\\\\n" )
                print ( line )
        
        self.summarize()
        self.latexFooter()



def main():
    import argparse
    argparser = argparse.ArgumentParser(description="meta statistics analyzer")
    argparser.add_argument ( '-d', '--dictfile', nargs='*',
            help='input dictionary file(s) [../data/database/]',
            type=str, default='.,/data/database/' )
    argparser.add_argument ( '-t', '--topos', nargs='*',
            help='filter for topologies, comma separated list or multiple arguments. prefix with ^ is negation. [None]',
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
    argparser.add_argument ( '-l', '--latex', help='create a latex version',
            action="store_true" )
    args=argparser.parse_args()
    analyzer = Analyzer ( args.dictfile, args.topos, args.nocolors, args.latex,
                          args.enumerate )
    analyzer.analyze ( args.nlargest, args.nsmallest )

if __name__ == "__main__":
    main()
