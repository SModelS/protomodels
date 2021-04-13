#!/usr/bin/env python3

""" plot the meta statistics of database.dict """

import matplotlib
matplotlib.use('agg')
from matplotlib import pyplot as plt
import numpy as np
import os, glob, pickle, sys
import scipy.stats
import matplotlib.mlab as mlab

class Plotter:
    def __init__ ( self, pathname, filtervalue: float, comment, likelihood: str, reset,
                   topologies, unscale, signalmodel, ntoys: int, filtersigma: float,
                   collaboration : str ):
        """
        :param filename: filename of dictionary
        :param filtervalue: filter out signal regions with expectedBG < filtervalue
        :param comment: an optional comment, to write in the plot
        :param likelihood: form of likelihood: "gauss", "gauss+poisson", or 
                           "lognormal+poisson"
                           "gauss" or "g" means only a Gaussian for everything
                           "gauss+poisson" or "gp" means Gauss * Poisson
                           "lognormal+poisson" or "lp" means Lognormal * Poisson
        :param reset: if true, then dont recycle pickle files
        :param topologies: if not Not, then filter for these topologies (e.g. T2tt)
        :param unscale: unscale, i.e. use the fudged bgError also for computing likelihoods
        :param signalmodel: use the signal+bg model for computing likelihoods
        :param ntoys: number of MC toys
        :param filtersigma: filter out signal regions with expectedBG/bgErr < filtersigma
        :param collaboration: select a specific collaboration
        """
        collaboration = collaboration.upper()
        if collaboration in [ "", "*" ]:
            collaboration = "ALL"
        if not collaboration in [ "CMS", "ATLAS", "ALL" ]:
            print ( "error: collaboration must be either CMS, ATLAS, or ALL." )
            sys.exit(-1)
        self.collaboration = collaboration
        abbreviations = { "g": "gauss", "gp": "gauss+poisson", "lp": "lognormal+poisson" }
        if likelihood in abbreviations:
            likelihood = abbreviations[likelihood]
        if likelihood not in [ "gauss", "gauss+poisson", "lognormal+poisson" ]:
            print ( "error, likelihood is to be one of: gauss, gauss+poisson, lognormal+poisson" )
            sys.exit()
        self.likelihood = likelihood ## False: gauss, True: lognormal
        self.ntoys = ntoys
        self.reset = reset
        self.unscale = unscale
        self.signalmodel = signalmodel
        self.topologies = []
        self.filtersigma = filtersigma
        if topologies not in [ None, "" ]:
           self.topologies = topologies.split(",")
        self.filenames = []
        if comment in [ "None", "", "none" ]:
            comment = None
        self.comment = comment
        for pname in pathname:
            if os.path.isdir ( pname ):
                pname = pname + "/db*dict"
            self.filenames += glob.glob ( pname )
        self.filter = filtervalue
        self.meta = {}
        self.data = {}
        self.read()

    def selectedCollaboration( self, anaid ):
        """ does anaid pass the collaboration selection? """
        if self.collaboration == "ALL":
            return True
        if self.collaboration in anaid:
            return True
        return False

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
                if not self.selectedCollaboration ( i ):
                    continue
                if "expectedBG" in v and v["expectedBG"]>=self.filter and \
                        v["expectedBG"]/v["bgError"]>=self.filtersigma:
                    newdata[i]=v
                else:
                    if i.endswith ( ":ul" ):
                        print ( f"[plotDBDict] removing {basename}:{i} (is an UL)" )
                    else:
                        eBG,bgerr=None,None
                        if "expectedBG" in v:
                            eBG = v["expectedBG"]
                            bgerr = v["bgError"]
                        print ( f"[plotDBDict] removing {basename}:{i} (eBG is {eBG}+-{bgerr})" )
            print ( f"[plotDBDict] keeping {len(newdata)}/{len(data)} for {basename}" )
            self.data[basename] = newdata

    def computeP ( self, obs, bg, bgerr, sigN ):
        """ compute p value, gaussian or lognormal nuisance model """
        #simple = False ## approximation as Gaussian
        simple = ( self.likelihood == "gauss" )
        loc = bg
        if self.signalmodel and sigN != None:
            loc = bg + sigN
        if simple:
            dn = obs - loc
            x = dn / np.sqrt ( bgerr**2 + loc )
            p = scipy.stats.norm.cdf ( x )
        else:
            return self.computePWithToys ( obs, bg, bgerr, sigN )
        return p

    def computePWithToys ( self, obs, bg, bgerr, sigN ):
        """ compute p value, for now we assume Gaussanity """
        fakes = []
        bigger = 0
        n = self.ntoys
        central = bg
        if self.signalmodel and sigN != None:
            central = bg + sigN
        if "lognormal" in self.likelihood and central > ( bgerr / 4. ):
            loc = central**2 / np.sqrt ( central**2 + bgerr**2 )
            stderr = np.sqrt ( np.log ( 1 + bgerr**2 / central**2 ) )
            lmbda = scipy.stats.lognorm.rvs ( s=[stderr]*n, scale=[loc]*n )
        else: ## Gauss
            lmbda = scipy.stats.norm.rvs ( loc=[central]*n, scale=[bgerr]*n )
            lmbda = lmbda[lmbda>0.]
            # lmbda = np.array( [x if x>0. else 0. for x in lmbda ] )
        fakeobs = scipy.stats.poisson.rvs ( lmbda )
        return ( sum(fakeobs>obs) + .5*sum(fakeobs==obs) ) / len(fakeobs)

    def getSqrts ( self, anaid ):
        """ get the sqrts of anaid """
        ret = 13
        t = anaid.replace("CMS-","").replace("ATLAS-","").replace("SUSY-","").\
                  replace("SUS-","").replace("PAS-","").replace("EXO-","").replace("CONF-","")
        t = t[:t.find("-")]
        t = int(t) % 2000
        if t < 15:
            ret = 8
        return ret

    def compute ( self, variable, fakeVariable, store ):
        """ compute the p-values """
        S,Sfake,P,Pfake={8:[], 13:[] },{8:[], 13:[] },{8:[], 13:[] },{8:[], 13:[] }
        for filename in self.filenames:
            selfbase = os.path.basename ( filename )
            fname = selfbase.replace(".dict",".pcl")
            hasPickle = False
            #if not self.reset:
            #    print ( f"[plotDBDict] looking for {fname}" )
            if False: # os.path.exists ( fname ) and not self.reset:
                print ( f"[plotDBDict] found {fname}. Using data therein." )
                with open ( fname, "rb" ) as f:
                    pname = os.path.basename ( pickle.load ( f ) )
                    topologies = pickle.load ( f )
                    unscale = pickle.load ( f )
                    signalmodel = pickle.load ( f )
                    if fname != pname or self.topologies != topologies or unscale != self.unscale or signalmodel != self.signalmodel:
                        print ( f"[plotDBDict] we want {fname} pickle has {pname}. Wont use." )
                    else:
                        tS = pickle.load ( f )
                        for k,v in tS.items():
                            S[k] += v
                        tSfake = pickle.load ( f )
                        for k,v in tSfake.items():
                            Sfake[k] += v
                        tP = pickle.load ( f )
                        for k,v in tP.items():
                            P[k] += v
                        tPfake = pickle.load ( f )
                        for k,v in tPfake.items():
                            Pfake[k] += v
                        f.close()
                        hasPickle = True
            if not hasPickle:
                #if self.reset:
                #    print ( f"[plotDBDict] or ordered reset: creating {fname}." )
                #else:
                #    print ( f"[plotDBDict] not found {fname}. Creating." )
                S_,Sfake_,P_,Pfake_={8:[], 13:[] },{8:[], 13:[] },{8:[], 13:[] },{8:[], 13:[] }
                data = self.data [ fname.replace(".pcl","") ]
                for k,v in data.items():
                    txns = []
                    if "txns" in v:
                        txns = v["txns"].split(",")
                    passesTx=False
                    if len(self.topologies)==0:
                        passesTx=True
                    for tx in self.topologies:
                        if tx in txns:
                            passesTx=True
                            break
                    if not passesTx:
                        print ( f"[plotDBDict] skipping {k}: does not pass Tx filter" )
                        continue

                    if not ":ul" in k:
                        s = v[variable]
                        sfake = v[fakeVariable]
                        sqrts = self.getSqrts ( k )
                        S[sqrts].append( s )
                        S_[sqrts].append ( s )
                        Sfake_[sqrts].append ( sfake )
                        obs = v["origN"]
                        if not "orig" in variable:
                            obs = v["newObs"]
                        fakeobs = v["newObs"]
                        vexp = v["expectedBG"]
                        bgErr = v["bgError"]/v["fudge"]
                        if self.unscale:
                            bgErr = v["bgError"]
                        if vexp < self.filter:
                            continue
                        if vexp / bgErr < self.filtersigma:
                            continue
                        sigN = None
                        if "sigN" in v:
                            sigN = v["sigN"]
                        # bgErr = v["bgError"]# /v["fudge"]
                        if "orig_p" in v:
                            p = v["orig_p"]
                        else:
                            p = self.computeP ( obs, vexp, bgErr, sigN )
                        P[sqrts].append( p )
                        P_[sqrts].append ( p )
                        if "new_p" in v:
                            pfake = v["new_p"]
                        else:
                            pfake = self.computeP ( fakeobs, vexp, bgErr, sigN )
                        Pfake[sqrts].append( pfake )
                        Pfake_[sqrts].append ( pfake )
                if False: # store:
                    print ( f"[plotDBDict] dumping to {fname}" )
                    with open ( fname, "wb" ) as f:
                        pickle.dump ( os.path.basename ( fname ), f )
                        pickle.dump ( self.topologies, f )
                        pickle.dump ( self.unscale, f )
                        pickle.dump ( self.signalmodel, f )
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
        mean,std = np.mean ( S[13]+S[8] ), np.std ( S[13]+S[8] )
        #minX, maxX = min(S), max(S)
        #x = np.linspace( minX, maxX,100 )
        # plt.legend()
        dbname = os.path.basename ( self.meta["database"] )
        title = f"SModelS database v{dbname}"
        # title = f"$p$-values, SModelS database v{dbname}"
        fudge = 1.
        if "fudge" in self.meta:
            fudge = self.meta["fudge"]
        if abs ( fudge - 1. ) > 1e-3:
            title += ", fudge=%.2f" % fudge
        if len (self.topologies )>0:
            title += f", selecting {','.join(self.topologies)}"
        if self.unscale:
            title += f" (unscaling)"
        if self.signalmodel:
            title += f" (signalmodel)"
        nbins = 10 ## change the number of bins
        fig, ax = plt.subplots()
        x = [ P[8], P[13] ]
        avgp8,stdp8=np.mean( P[8] ),np.std( P[8] )
        bin8=int(avgp8*nbins)
        avgp13,stdp13=np.mean( P[13] ),np.std( P[13] )
        bin13=int(avgp13*nbins)
        nm1 = 1. / len(self.filenames)
        weights = [ [ nm1 ]*len(P[8]), [ nm1 ]*len(P[13]) ]
        bins = np.arange ( 0., 1+1e-7, 1/nbins )
        H1 = plt.hist ( x, weights = weights, bins=bins, histtype="bar",
                   label=[ "real, 8 TeV", "real, 13 TeV" ], color=[ "tab:green", "tab:blue" ], stacked=True )
        h8 = H1[0][0][bin8]
        h13 = H1[0][1][bin13]
        l13 = H1[0][0][bin13]
        l8 = plt.plot ( [ avgp8, avgp8 ], [0, h8 ], color = "darkgreen", zorder=1 )
        l13 = plt.plot ( [ avgp13, avgp13 ], [ l13, h13 ], color = "darkblue", zorder=1 )
        fweights = [ nm1 ]*len(Pfake[8]+Pfake[13])
        # fweights = [ [ nm1 ]*len(Pfake[8]), [ nm1 ]*len(Pfake[13]) ]
        H2 = plt.hist ( [ Pfake[8] + Pfake[13] ], weights = fweights,
                        bins=bins, stacked=True, zorder=9,
                        label="fake", color=["red" ], linewidth=3, histtype="step" )
        #plt.hist ( Pfake[8]+Pfake[13], weights = fweights, bins=nbins, 
        #           label="fake", edgecolor="red", linewidth=3, histtype="step" )
        print ( "real Ps %d entries at %.3f +/- %.2f" % 
                ( len(P[13]), np.mean(P[13]), np.std(P[13])  ) )
        print ( "fake Ps %d entries at %.3f +/- %.2f" % 
                ( len(Pfake[13]), np.mean(Pfake[13]), np.std(Pfake[13]) ) )
        plt.legend()
        if self.likelihood == "lognormal+poisson":
            title += " (lognormal)"
        if self.likelihood == "gauss":
            title += " (simple)"
        if self.collaboration != "ALL":
            title += f" {self.collaboration} only"
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
            help='input dictionary file(s) [../data/database/]',
            type=str, default='.,/data/database/' )
    argparser.add_argument ( '-o', '--outfile', nargs='?',
            help='output file [./pDatabase.png]',
            type=str, default='./pDatabase.png' )
    argparser.add_argument ( '-c', '--comment', nargs='?',
            help='an optional comment, to put in the plot [None]',
            type=str, default=None )
    argparser.add_argument ( '-r', '--reset', 
            help='reset, dont recycle pickle files', action='store_true' )
    argparser.add_argument ( '-u', '--unscale', 
            help='unscale, i.e. use the fudged bgError also for computing likelihoods', action='store_true' )
    argparser.add_argument ( '-S', '--signalmodel', 
            help='use the signal+bg model for computing likelihoods', action='store_true' )
    argparser.add_argument ( '-l', '--likelihood', nargs='?',
            help='likelihood: gauss (g), gauss+poisson (gp), or lognormal+poisson (lp) [gauss+poisson]',
            type=str, default="gauss+poisson" )
    argparser.add_argument ( '-t', '--topologies', nargs='?',
            help='filter for certain topologies, e.g. T1, T2tt. Comma separated. The signal region must have a map for any one of the given topologies. [None]',
            type=str, default=None )
    argparser.add_argument ( '-f', '--filter', nargs='?',
            help='filter out signal regions with expectedBG<x [x=3.]',
            type=float, default=3.5 )
    argparser.add_argument ( '-s', '--filtersigma', nargs='?',
            help='filter out signal regions with expectedBG/bgErr<x [x=2.]',
            type=float, default=3. )
    argparser.add_argument ( '-C', '--select_collaboration', nargs='?',
            help='select a specific collaboration CMS, ATLAS, all [all]',
            type=str, default="all" )
    argparser.add_argument ( '-n', '--ntoys', nargs='?',
            help='number of MC toys to throw [50000]',
            type=float, default=50000 )
    args=argparser.parse_args()
    plotter = Plotter ( args.dictfile, args.filter, args.comment, args.likelihood, 
                        args.reset, args.topologies, args.unscale, args.signalmodel,
                        args.ntoys, args.filtersigma, args.select_collaboration )
    plotter.plot( "origS", "S", args.outfile )

if __name__ == "__main__":
    main()
