#!/usr/bin/env python3

""" plot the meta statistics of database.dict """

from matplotlib import pyplot as plt
import numpy as np
import os, glob, pickle
import scipy.stats
import matplotlib.mlab as mlab

class Plotter:
    def __init__ ( self, pathname, filtervalue: float, comment, lognormal = True ):
        """
        :param filename: filename of dictionary
        :param filtervalue: filter out signal regions with expectedBG < filtervalue
        :param comment: an optional comment, to write in the plot
        :param lognormal: if False, use Gauss, else lognormal
        """
        self.lognormal = lognormal ## False: gauss, True: lognormal
        self.filenames = []
        if comment in [ "None", "", "none" ]:
            comment = None
        self.comment = comment
        for pname in pathname:
            self.filenames += glob.glob ( pname )
        self.filter = filtervalue
        self.meta = {}
        self.data = {}
        self.read()

    def read ( self ):
        """ read in content of filename """
        for fname in self.filenames:
            with open( fname,"rt") as f:
                tmp=f.readlines()
            lines = []
            for line in tmp:
                if line.startswith("#"):
                    continue
                lines.append ( line )
            basename = os.path.basename ( fname ).replace(".dict","")
            self.meta.update (  eval(lines[0]) )
            nan=float("nan")
            data = eval("\n".join(lines[1:]))
            newdata = {}
            for i,v in data.items():
                if "expectedBG" in v and v["expectedBG"]>=self.filter:
                    newdata[i]=v
                else:
                    if i.endswith ( ":ul" ):
                        print ( f"[plotDBDict] removing {basename}:{i} (is an UL)" )
                    else:
                        eBG=None
                        if "expectedBG" in v:
                            eBG = v["expectedBG"]
                        print ( f"[plotDBDict] removing {basename}:{i} (eBG is {eBG})" )
            print ( f"[plotDBDict] keeping {len(newdata)}/{len(data)} for {basename}" )
            self.data[basename] = newdata

    def computeP ( self, obs, bg, bgerr ):
        """ compute p value, for now we assume Gaussanity """
        simple = False ## approximation as Gaussian
        if simple:
            x = (obs - bg ) / np.sqrt ( bgerr**2 + bg )
            p = scipy.stats.norm.cdf ( x )
        else:
            return self.computePWithToys ( obs, bg, bgerr )
        return p

    def computePWithToys ( self, obs, bg, bgerr ):
        """ compute p value, for now we assume Gaussanity """
        fakes = []
        bigger = 0
        n= 10000
        if self.lognormal: ## lognormal
            loc = bg**2 / np.sqrt ( bg**2 + bgerr**2 )
            stderr = np.sqrt ( np.log ( 1 + bgerr**2 / bg**2 ) )
            lmbda = scipy.stats.lognorm.rvs ( s=[stderr]*n, scale=[loc]*n )
        else: ## Gauss
            lmbda = scipy.stats.norm.rvs ( loc=[bg]*n, scale=[bgerr]*n )
            lmbda = lmbda[lmbda>0.]
        fakeobs = scipy.stats.poisson.rvs ( lmbda )
        return sum(fakeobs>obs) / len(fakeobs)

    def compute ( self, variable, fakeVariable, store ):
        """ compute the p-values """
        S,Sfake,P,Pfake=[],[],[],[]
        for filename in self.filenames:
            selfbase = os.path.basename ( filename )
            fname = selfbase.replace(".dict",".pcl")
            print ( f"[plotDBDict] looking for {fname}" )
            hasPickle = False
            if os.path.exists ( fname ):
                print ( f"[plotDBDict] found {fname}. Using data therein." )
                with open ( fname, "rb" ) as f:
                    pname = os.path.basename ( pickle.load ( f ) )
                    filtervalue = pickle.load ( f )
                    if fname != pname or abs (self.filter - filtervalue ) > 1e-6:
                        print ( f"[plotDBDict] we want {fname} pickle has {pname}. Wont use." )
                    else:
                        S += pickle.load ( f )
                        Sfake += pickle.load ( f )
                        P += pickle.load ( f )
                        Pfake += pickle.load ( f )
                        f.close()
                        hasPickle = True
            if not hasPickle:
                print ( f"[plotDBDict] not found {fname}. Creating." )
                S_,Sfake_,P_,Pfake_=[],[],[],[]
                data = self.data [ fname.replace(".pcl","") ]
                for k,v in data.items():
                    if not ":ul" in k:
                        s = v[variable]
                        sfake = v[fakeVariable]
                        S.append( s )
                        S_.append ( s )
                        Sfake.append( sfake )
                        Sfake_.append ( sfake )
                        obs = v["origN"]
                        if not "orig" in variable:
                            obs = v["newObs"]
                        fakeobs = v["newObs"]
                        vexp = v["expectedBG"]
                        if vexp < self.filter:
                            continue
                        p = self.computeP ( obs, vexp, v["bgError"] )
                        P.append( p )
                        P_.append ( p )
                        pfake = self.computeP ( fakeobs, vexp, v["bgError"] )
                        Pfake.append( pfake )
                        Pfake_.append ( pfake )
                        P.append( scipy.stats.norm.cdf ( s ) )
                        cfake = scipy.stats.norm.cdf ( sfake )
                        Pfake.append( cfake )
                if store:
                    print ( f"[plotDBDict] dumping to {fname}" )
                    with open ( fname, "wb" ) as f:
                        pickle.dump ( os.path.basename ( fname ), f )
                        pickle.dump ( self.filter, f )
                        pickle.dump ( S_, f )
                        pickle.dump ( Sfake_, f)
                        pickle.dump ( P_, f )
                        pickle.dump ( Pfake_, f )
                        f.close()
        # print ( "P", P )
        return S,Sfake,P,Pfake

    def plot( self, variable, fakeVariable, outfile ):
        """ plot the p-values """
        S,Sfake,P,Pfake=self.compute ( variable, fakeVariable, True )
        mean,std = np.mean ( S), np.std ( S )
        #minX, maxX = min(S), max(S)
        #x = np.linspace( minX, maxX,100 )
        # plt.legend()
        dbname = os.path.basename ( self.meta["database"] )
        title = f"$p$-values, SModelS database v{dbname}"
        fudge = 1.
        if "fudge" in self.meta:
            fudge = self.meta["fudge"]
        if abs ( fudge - 1. ) > 1e-3:
            title += ", fudge=%.2f" % fudge
        nbins = 10 ## change the number of bins
        fig, ax = plt.subplots()
        plt.hist ( P, weights = [ 1. / len(self.filenames) ]*len(P), bins=nbins,
                   label="real", facecolor="tab:blue" )
        plt.hist ( Pfake, weights = [ 1. / len(self.filenames) ]*len(P), bins=nbins,
                   label="fake", edgecolor="red", linewidth=3, histtype="step" )
        print ( "real Ps %d entries at %.3f +/- %.2f" % 
                ( len(P), np.mean(P), np.std(P)  ) )
        print ( "fake Ps %d entries at %.3f +/- %.2f" % 
                ( len(Pfake), np.mean(Pfake), np.std(Pfake) ) )
        plt.legend()
        plt.title  ( title )
        plt.xlabel ( "$p$-values" )
        plt.ylabel ( "# Signal Regions" )
        print ( f"[plotDBDict.py] plotting {outfile}"  )
        if self.comment != None:
            plt.text ( .65, -.11, self.comment, transform=ax.transAxes, 
                       style="italic" )
        plt.savefig ( outfile )
        plt.clf()

def main():
    import argparse
    argparser = argparse.ArgumentParser(description="meta statistics plotter, i.e. the thing that plots pDatabase.png")
    argparser.add_argument ( '-d', '--dictfile', nargs='*',
            help='input dictionary file [./database.dict]',
            type=str, default='./database.dict' )
    argparser.add_argument ( '-o', '--outfile', nargs='?',
            help='output file [./pDatabase.png]',
            type=str, default='./pDatabase.png' )
    argparser.add_argument ( '-c', '--comment', nargs='?',
            help='an optional comment, to put in the plot [None]',
            type=str, default="(lognormal)" )
    argparser.add_argument ( '-f', '--filter', nargs='?',
            help='filter out signal regions with expectedBG<x [x=-1.]',
            type=float, default=-1. )
    args=argparser.parse_args()
    plotter = Plotter ( args.dictfile, args.filter, args.comment, True )
    plotter.plot( "origS", "S", args.outfile )

if __name__ == "__main__":
    main()
