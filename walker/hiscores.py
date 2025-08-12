#!/usr/bin/env python3

""" A class that centralizes access to the hiscore list over multiple threads.
"""

__all__ = [ "Hiscores" ]

import random, copy, pickle, os, fcntl, time, subprocess, colorama
import numpy as np
from scipy import stats
from csetup import setup
setup()
from builder.manipulator import Manipulator
from tester.combiner import Combiner
from ptools import helpers
from ptools import sparticleNames
from ptools.helpers import py_dumps
from typing import Union
from os import PathLike
from base.loggerbase import LoggerBase

class Hiscores ( LoggerBase ):
    """ encapsulates the hiscore list. """
    def __init__ ( self, walkerid: int = 0, save_hiscores: bool = False,
                   picklefile: PathLike="hiscores.cache", backup : bool = True, keep_separate_hiscores = False,
                   hiscores = None, predictor = None ):
        """ the constructor
        :param save_hiscores: if true, then assume you want to save, not just read.
        :param picklefile: path of pickle file name to connect hiscore list with
        :param backup: if True, make a backup pickle file old_<name>.pcl
        :param hiscores: if None, try to get them from file, if a list,
                         then these are the hiscore protomodels.
        """
        super ( Hiscores, self ).__init__ ( walkerid )
        self.walkerid = walkerid
        self.save_hiscores = save_hiscores
        self.backup = backup ## backup hiscore lists?
        self.nkeep = 3 ## how many do we keep.
        self.hiscores = [ None ]*self.nkeep
        self.predictor = predictor
        self.fileAttempts = 0 ## unsucessful attempts at reading or writing
        self.pickleFile = picklefile
        self.mtime = 0 ## last modification time of current list
        self.namer = sparticleNames.SParticleNames ( susy = False )
        self.keep_separate_hiscores = keep_separate_hiscores
        if hiscores == None:
            self.updateListFromPickle ( )
        else:
            self.hiscores = hiscores
            self.mtime = time.time()

    @classmethod
    def writeDictionariesToFile ( cls, filename : os.PathLike,
           objs : list ) -> bool:
        """ class method, write the dictionaries objs in a 
        formatted manner to file filename

        :returns: True if worked
        """
        ds = []
        for obj in objs:
            d = py_dumps ( obj, level = 1 )
            d = " "*4 + d
            ds.append ( d )
        from ptools.locker import lock, unlock
        lock ( filename )
        with open ( filename, "wt" ) as f:
            f.write("[\n")
            first = True
            for d in ds:
                if not first:
                    f.write ( ",\n" )
                f.write ( f"{d}" )
                first = False
            f.write("\n]\n")
            f.close()
        unlock ( filename )
        return True
        
    def currentMinTL ( self ):
        """ the current minimum TL to make it into the list. """
        if self.hiscores[-1] == None:
            return 0.
        return self.hiscores[-1].TL

    def currentMinK ( self, zeroIsMin=False ):
        """ the current minimum K to make it into the list.
        :param zeroIsMin:  if false, min k can become negative
        """
        if self.hiscores[-1] == None or len(self.hiscores)<10:
            if zeroIsMin:
                return 0.
            return -30.
        mk = -10.
        if hasattr ( self.hiscores[-1], "K" ):
            mk = self.hiscores[-1].K
            if mk == None:
                mk = -10.
        if zeroIsMin:
            return max ( mk, 0. )
        return mk

    def globalMaxTL ( self ):
        """ globally (across all walkers), the highest TL """
        ret = 0.
        if self.hiscores[0] != None:
            if self.hiscores[0].TL > ret:
                ret = self.hiscores[0].TL
        TLoldfile = "TLold.conf"
        if os.path.exists ( TLoldfile ):
            with open ( TLoldfile, "rt" ) as f:
                lines = f.readlines()
                if len(lines)>0:
                    ret = float(lines[0])
                f.close()
        return ret

    def globalMaxK ( self ):
        """ globally (across all walkers), the highest K """
        ret = -3. ## set to negative if no hiscore exists
        if self.hiscores[0] != None:
            if self.hiscores[0].K > ret:
                ret = self.hiscores[0].K
        Koldfile = "Kold.conf"
        if os.path.exists ( Koldfile ):
            with open ( Koldfile, "rt" ) as f:
                lines = f.readlines()
                if len(lines)>0:
                    ret = float(lines[0])
                f.close()
        return ret

    def globalMinK ( self ):
        """ the minimum K needed to make it into hiscore list """
        Kminfile = "Kmin.conf"
        ret = -99
        if not os.path.exists ( Kminfile ):
            return ret
        with open ( Kminfile, "rt" ) as f:
            lines = f.readlines()
            if len(lines)>0:
                ret = float(lines[0])
            f.close()
        return ret

    def similarDicts ( self, a : dict, b : dict ) -> bool:
        """ are models a and b similar? 
        
        :returns: True if similar
        """
        dK = 0.
        if a["K"] is None and b["K"] is not None:
            return False
        if b["K"] is None and a["K"] is not None:
            return False
        if a["K"] is not None and b["K"] is not None:
            dK = abs ( a["K"] - b["K"] )
        if dK > 1e-5:
            return False
        dTL = abs ( a["TL"] - b["TL"] )
        if dTL > 1e-5:
            return False
        if a["masses"].keys() != b["masses"].keys():
            return False
        massDiff = [ abs(m-b["masses"][pid])/m for pid,m in a["masses"].items() if m ]
        if max(massDiff) > 1e-5:
            return False
        return True

    def insertHiscore ( self, L : list, hi : dict ):
        """ insert hiscore <hi> into list <L> at the appropriate place """
        K = hi["K"]
        ret = []
        for oldhi in L: ## as long as the Ks are above the new K, append
            if self.similarDicts ( oldhi, hi ):
                ## already exists in list? skip insertion
                self.log("Protomodel already exists in hiscore file. Skip")
                return L, False
            if K is None: # if K is None, append them all
                ret.append ( oldhi )
            elif oldhi["K"] is not None and oldhi["K"] >= K:
                ret.append ( oldhi ) # K is not None, oldK is greater
            else: # oldK is None or oldK < K
                break
        ## now the new hiscore
        ret.append ( hi )
        ## now fill up
        for oldhi in L[len(ret)-1:]:
            ret.append ( oldhi )
        ret = ret[:10] ## cut off, max ten
        return ret, True

    @classmethod
    def fromDictionaryFile ( cls, path : PathLike,
           firstn : Union[None,int] = 0, 
           dbpath : PathLike = "official",
           walkerid : Union[str,int] = 0 ):
        """ initialise from a dictionary file

        :param path: filename of .dict file
        :param firstn: initialise only first n entries
        :param dbpath: path to database
        :param walkerid: log everything as walker #walkerid
        :returns: Hiscores object
        """
        assert firstn == 0, "firstn != 0 not yet working"
        from tester.predictor import Predictor
        predictor = Predictor(walkerid, do_srcombine=True, dbpath = dbpath )
        hiscores = []
        c = 0
        while True:
            m = Manipulator( path, nth = c, walkerid = walkerid )
            m.M.walkerid = walkerid
            predictor.predict ( m, keep_predictions=True )
            hiscores.append ( m.M )
            c+=1
            if type(firstn) == int and c > firstn:
                break
        return cls ( hiscores= hiscores, predictor = predictor,
                     walkerid = walkerid )

        # assert False, "implement me"
    
    def updateGlobalHiscoreFile ( self, m : Manipulator,
           hiscorefile : PathLike = "hiscores_global.dict" ) -> bool:
        """
        we have a new hiscore, add to hiscores.dict
        :param m: manipulator
        :param hiscorefile: the hiscore dict file to update
        :returns: true, if successful
        """
        oldhiscores=[]
        if os.path.exists ( hiscorefile ):
            tryRead=0
            success=False
            ## stop loop at success or when tryRead is at least 5
            while (not success) and tryRead<5:
                tryRead+=1
                try:
                    with open ( hiscorefile, "rt" ) as h:
                        txt = h.read()
                        txt = txt.replace("inf","float('inf')" )
                        oldhiscores = eval( txt )
                        h.close()
                        success=True
                except SyntaxError as e:
                    time.sleep( .1+3*tryRead )
        D=m.writeDictFile ( None, cleanOut = False, ndecimals = 6 )
        newlist, added = self.insertHiscore ( oldhiscores, D )
        self.writeListToDictFile ( hiscorefile, newlist )
        """ ## old version
        self.log (f"Write model to {hiscorefile}" )
        with open ( hiscorefile, "wt" ) as f:
            f.write ( "[\n" )
            for ctr,l in enumerate(newlist):
                f.write ( f"{l}" )
                #if ctr < len(newlist)-1:
                #    f.write ( ",\n" % ( l ) )
            f.write ( "\n]\n" )
            f.close()
        """
        with open ( "Kold.conf", "wt" ) as f:
            f.write ( f"{m.M.K}\n" )
            f.close()
        with open ( "Kmin.conf", "wt" ) as f:
            f.write ( f"{newlist[-1]['K']}\n" )
            f.close()
        return True
    
    def updateTopHiscoreFile ( self, m : Manipulator,
           hiscorefile : PathLike = "hiscores_top.dict" ) -> bool:
        """
        Update the top hiscore model from each walk
        :param m: manipulator
        :param hiscorefile: the hiscore dict file to update
        :returns: true, if successful
        """
        hiscore_top = []
        import glob
        high_files = glob.glob("all_hiscores/hiscore*.dict")
        for file in high_files:
            if os.path.exists ( file ):
                tryRead=0
                success=False
                ## stop loop at success or when tryRead is at least 5
                while (not success) and tryRead<5:
                    tryRead+=1
                    try:
                        with open ( file, "rt" ) as h:
                            txt = h.read()
                            hiscores = eval( txt )
                            hiscore_top.append(hiscores[0])
                            h.close()
                            success=True
                    except SyntaxError as e:
                        if tryRead > 10:
                            raise e
                        time.sleep( .1+3*tryRead )
        
        newlist = sorted(hiscore_top, key=lambda val:val['K'], reverse=True)
        self.writeListToDictFile ( hiscorefile, newlist )
        """
        self.log(f"Updating  {hiscorefile}" )
        # self.writeListToDictFile ( hiscorefile, newlist ) # FIXME use this
        with open ( hiscorefile, "wt" ) as f:
            f.write ( "[" )
            for ctr,l in enumerate(newlist):
                f.write ( f"{l}" )
                if ctr < len(newlist)-1:
                    f.write ( ",\n" % ( l ) )
            f.write ( "]\n" )
            f.close()
        """
        return True
        
    def updateHiscoreFile ( self, m : Manipulator,
           hiscorefile : PathLike = "hiscores.dict" ) -> bool:
        """
        we have a new hiscore, add to hiscores.dict
        :param m: manipulator
        :param hiscorefile: the hiscore dict file to update
        :returns: true, if successful
        """
        oldhiscores=[]
        helpers.mkdir ( "all_hiscores" )
        if os.path.exists ( f"all_hiscores/{hiscorefile}" ):
            tryRead=0
            success=False
            ## stop loop at success or when tryRead is at least 5
            while (not success) and tryRead<5:
                tryRead+=1
                try:
                    with open ( f"all_hiscores/{hiscorefile}" , "rt" ) as h:
                        txt = h.read()
                        oldhiscores = eval( txt )
                        h.close()
                        success=True
                except SyntaxError as e:
                    time.sleep( .1+3*tryRead )
        D=m.writeDictFile ( None, cleanOut = False, ndecimals = 6 )
        newlist, added = self.insertHiscore ( oldhiscores, D )
        if added:
            self.log(f"Write model to {hiscorefile}" )
            with open ( f"all_hiscores/{hiscorefile}", "wt" ) as f:
                f.write ( "[" )
                for ctr,l in enumerate(newlist):
                    f.write ( f"{l}" )
                    if ctr < len(newlist)-1:
                        f.write ( ",\n" % ( l ) )
                f.write ( "]\n" )
                f.close()
            return True
        return False

    def addResult ( self, ma ):
        """ add a result to the list
        :param ma: the manipulator object
        :returns: true, if result was added
        """
        
        #if ma.M.K < self.currentMinK():        #SN: removed zeroIsMin for now
        #    self.log(f"K {ma.M.K} less than Min K {self.currentMinK()}. Not adding to hiscore list.")
        #    return False ## doesnt pass minimum requirement
        #if ma.M.K == 0.:
        #    return False ## just to be sure, should be taken care of above, though

        # Kold = self.globalMaxK()
        #Kmin = self.globalMinK()
        # self.pprint ( f"adding results Kold is {Kold} Knew is {ma.M.K}" )
        ## FIXME we should only write into this file in the first maxstep/3 steps
        #if ma.M.K > Kmin:              #FIXME!
        # self.pprint ( "WARNING we shouldnt write into hiscore file afte maxstep/3 steps!!" )
        added = self.updateHiscoreFile( ma, hiscorefile = f"hiscores{self.walkerid}.dict")
        if added:
            self.updateGlobalHiscoreFile( ma)
            self.updateTopHiscoreFile(ma)
        ## we have a new hiscore?
        ## compute the particle contributions
        #if not hasattr ( ma.M, "particleContributions" ):
        #    self.pprint ( "particleContributions missing, compute them!" )
        #    self.computeParticleContributions(m)
        ## compute the analysis contributions
        #if not hasattr ( ma.M, "analysisContributions" ):
        #    self.pprint ( "analysisContributions missing, compute them!" )
        #    self.computeAnalysisContributions(m)
        protomodel = ma.M
        protomodel.getXsecs() #Make sure cross-sections have been computed

        for i,mi in enumerate(self.hiscores):
            if mi!=None and mi.almostSameAs ( ma.M ):
                ### this ma.M is essentially the ma.M in hiscorelist.
                ### Skip!
                self.pprint ( "the protomodel seems to be already in highscore list. skip" )
                return False

            if mi==None or ma.M.K > mi.K: ## ok, <i>th best result!
                self.demote ( i )
                self.hiscores[i] = copy.deepcopy ( ma.M )
                self.hiscores[i].cleanBestCombo( )
                break
        return True

    def computeParticleContributions ( self, manipulator ):
        """ this function sequentially removes all particles to compute
            their contributions to K """
        if False:
            from smodels.base import runtime
            runtime._experimental = True

        #Make sure the model is backep up
        manipulator.backupModel()

        unfrozen = manipulator.M.unFrozenParticles( withLSP=False )
        oldTL = manipulator.M.TL
        oldK = manipulator.M.K
        particleContributions = {} ## save the scores for the non-discarded particles.
        #particleContributionsTL = {} ## save the scores for the non-discarded particles, TLs

        #Make sure predictor is accesible
        if not self.predictor:
            self.pprint( "asked to compute particle contributions to score, but predictor has not been set")
            return

        pidsnmasses = [ (x,manipulator.M.masses[x]) for x in unfrozen ]
        pidsnmasses.sort ( key=lambda x: x[1], reverse=True )
        for cpid,(pid,mass) in enumerate(pidsnmasses):
            self.pprint ( f"computing contribution of {self.namer.asciiName(pid)} ({manipulator.M.masses[pid]:.1f}): [{int(cpid + 1)}/{len(unfrozen)}]" )

            #Remove particle and recompute SLHA file:
            manipulator.freezeParticle(pid, force=True )
            #Recompute cross-secions:
            manipulator.M.getXsecs()
            manipulator.M.K = 0.0
            manipulator.M.TL = 0.0
            self.predictor.predict( manipulator.M )
            if manipulator.M.K is None:
                self.pprint ( f"when removing {self.namer.asciiName(pid)}, K could not longer be computed. Setting to zero")
                manipulator.M.K = 0.0
                manipulator.M.TL = 0.0
            if oldK <= 0:
                percK = 0.
            else:
                percK = ( manipulator.M.K - oldK ) / oldK
                self.pprint ( "when removing %s, K changed: %.3f -> %.3f (%.1f%s), TL: %.3f -> %.3f (%d evts)" % \
                    ( self.namer.asciiName(pid), oldK, manipulator.M.K, 100.*percK, "%", oldTL,manipulator.M.TL, manipulator.M.nevents ) )

            #Store the new TL and K values in the original model:
            particleContributions[pid]=manipulator.M.K
            #particleContributionsTL[pid]=manipulator.M.TL
            #Make sure to restore the model to its initial (full particle content) state
            manipulator.restoreModel()
            #Store contributions in the protomodel:
            manipulator.M.particleContributions = particleContributions
            #manipulator.M.particleContributionsTL = particleContributionsTL

        self.pprint ( f"stored {len(manipulator.M.particleContributions)} particle contributions" )

    def computeAnalysisContributions( self, manipulator ):
        """ compute the contributions to TL of the individual analyses
        :returns: the model with the analysic constributions attached as
                  .analysisContributions
        """

        try:
            self.pprint ( "Now computing analysis contributions" )
            self.pprint ( f"Recompute the score. Old one at K={manipulator.M.K:.3f}"
                          f", TL={manipulator.M.TL:.2f}" )
            contributionsTL = {}
            contributionsK = {}
            combiner = Combiner()
            dTLtot, dKtot = 0., 0.
            bestCombo = copy.deepcopy ( manipulator.M.bestCombo )
            #self.pprint ( "we have %d entries in best combo" % len(bestCombo) )
            prior = combiner.computePrior ( manipulator.M )
            #self.pprint ( "the prior is %s" % prior )
            for ctr,pred in enumerate(bestCombo):
                #self.pprint ( "Now starting to compute for %d" % ctr )
                combo = bestCombo[:ctr]+bestCombo[ctr+1:]
                # combo = copy.deepcopy ( bestCombo )[:ctr]+copy.deepcopy ( bestCombo)[ctr+1:]
                #self.pprint ( "deep copy still worked: %d" % (len(combo)) )
                TL, muhat_ = combiner.getSignificance ( combo )
                #self.pprint ( "TL for %d is %s" % ( ctr, TL ) )
                K = combiner.computeK ( TL, prior )
                #self.pprint ( "K for %d is %s" % ( ctr, K ) )
                contributionsK [ ctr ] = K
            self.pprint ( "finished computing contributions" )

            contrsWithNames = {}
            for k,v in contributionsK.items():
                # self.pprint ( "contributionsK of %s reads %s" % ( k, v ) )
                contrsWithNames [ manipulator.M.bestCombo[k].analysisId() ] = v
            manipulator.M.analysisContributions = contrsWithNames
            self.pprint ( f"stored {len(manipulator.M.analysisContributions)} analyses contributions" )
        except Exception as e:
            self.pprint ( f"in computeAnalysisContributions caught {str(e)}" )

    def demote ( self, i ):
        """ demote everything from i+1 on,
            i.e (i+1)->(i+2), (i+2)->(i+3) and so on """
        for j in range(self.nkeep-1,i,-1):
            m = copy.deepcopy ( self.hiscores[j-1] )
            self.hiscores[j]= m
        if len(self.hiscores)>self.nkeep:
            self.hiscores = self.hiscores[:self.nkeep]

    def updateListFromPickle ( self ):
        """ fetch the list from the pickle file """
        if not os.path.exists ( self.pickleFile ) or \
            os.stat ( self.pickleFile ).st_size < 100:
            return
        mtime = os.stat ( self.pickleFile ).st_mtime
        if mtime > 0 and mtime == self.mtime:
            ## no modification. return
            return

        try:
            with open( self.pickleFile,"rb") as f:
                try:
                    #fcntl.flock ( f, fcntl.LOCK_EX | fcntl.LOCK_NB )
                    self.hiscores = pickle.load ( f )
                    self.timestamp = "?"
                    try:
                        self.timestamp = pickle.load ( f )
                    except EOFError:
                        pass
                    #fcntl.flock ( f, fcntl.LOCK_UN )
                    f.close()
                except (BlockingIOError,OSError) as e:
                    ## make sure we dont block!
                    #fcntl.flock( f, fcntl.LOCK_UN )
                    raise e
            self.mtime = mtime
            nhs = 0
            for i in self.hiscores:
                if i != None:
                    nhs += 1
            pfname = helpers.simplifyUnixPath ( self.pickleFile )
            self.pprint ( f"loaded {nhs} hiscores from {pfname}" )
            self.fileAttempts=0
        except Exception as e:
        # except OSError or BlockingIOError or EOFError or pickle.UnpicklingError or TypeError as e:
            self.fileAttempts+=1
            if self.fileAttempts<20: # try again
                self.pprint ( f"Exception[X] {e!s}: type({type(e)}), Waiting for {self.pickleFile} file, {int(self.fileAttempts)}" )
                time.sleep ( (.2 + np.random.uniform(0.,1.))*self.fileAttempts )
                self.updateListFromPickle()
                self.pprint ( "Loading hiscores worked this time" )
            else:
                self.pprint ( "Timed out when try to get hiscores!" )

    def clean ( self ):
        """ clean hiscore list, i.e. remove cruft from protomodels.
            leave first one as it is """
        for ctr,h in enumerate(self.hiscores[1:]):
            if h != None:
                m=Manipulator ( h )
                m.rescaleSignalBy(m.M.muhat)
                m.delBackup ( )
                m.M.cleanBestCombo ()
                self.hiscores[ctr+1]=m.M

    def writeListToDictFile ( self, dictFile : Union[None,str] = None,
           listofhiscores : Union[None,list]  = None ):
        """ write the models in append mode in a single dictFile.
        :param dictFile: write to dictFile. If None, then self.pickleFile
        is used, but with ".dict" as extension.
        :param listofhiscores: either explicit list or None, in which case
        we use self.hiscores
        """
        self.log (f"Writing model to {dictFile}" )
        if dictFile==None:
            dictFile = self.pickleFile
        if dictFile.endswith(".cache"):
            dictFile = f"{dictFile[:-6]}.dict"
        if listofhiscores == None:
            listofhiscores = self.hiscores
        self.writeDictionariesToFile ( dictFile, listofhiscores )
        # for protomodel in listofhiscores:
# ma = Manipulator ( protomodel, initTestStats = True )
#            ma.writeDictFile ( outfile = dictFile, cleanOut=False,appendMode=True )
        #    Manipulator.writeDictionaryToFile ( dictFile, protomodel,
        #                                        appendMode=True )

    def writeListToPickle ( self, pickleFile : Union[None,str]=None,
            check : bool = True ):
        """ pickle the hiscore list.
        :param pickleFile: write to pickleFile. If None, then self.pickleFile
            is used.
        :param check: perform a check whether the file has changed?
        """
        if len ( self.hiscores ) == 0:
            self.log ( "hiscore list is empty will not write out" )
            return
        onlyNones = True
        for i in self.hiscores:
            if i != None:
                onlyNones = False
                break
        if onlyNones:
            self.log ( "hiscore list contains only nones" )
            return
        if pickleFile==None:
            pickleFile = self.pickleFile
        if check and os.path.exists ( self.pickleFile ):
            mtime = os.stat ( self.pickleFile ).st_mtime
            if mtime > self.mtime:
                self.pprint ( "while writing to pickle file I see that it has changed" )
                self.updateListFromPickle()
                return False
        shortname = helpers.simplifyUnixPath(pickleFile)
        self.log ( f"saving new hiscore list to {shortname}" )
        try:
            if self.backup and os.path.exists ( pickleFile ):
                subprocess.getoutput ( f"mv -f {pickleFile} old_{pickleFile}" )
            # self.clean()
            with open( pickleFile, "wb" ) as f:
                fcntl.flock ( f, fcntl.LOCK_EX )
                pickle.dump ( self.hiscores, f )
                pickle.dump ( time.asctime(), f )
                fcntl.flock ( f, fcntl.LOCK_UN )
                f.close()
            self.mtime = os.stat ( pickleFile ).st_mtime
            self.fileAttempts=0
            return True
        except OSError or BlockingIOError as e:
            self.fileAttempts+=1
            if self.fileAttempts>2:
                self.pprint ( f"error when writing ({self.fileAttempts}th attempt) pickle file {pickleFile} ({shortname}): {e}" )
            if self.fileAttempts<5: # try again
                time.sleep ( .2 )
                self.writeListToPickle( pickleFile, check )
            return False
        return False

    def newResult ( self, ma ):
        """ see if new result makes it into hiscore list. If yes, then add.
        :param ma: the manipulator object
        :returns: true, if it entered the hiscore list
        """
        def pprint ( value ):
            if type(value) in [ float, np.float64 ]:
                return f"{value:.2f}"
            return str(value)
        K = pprint ( ma.M.K )
        TL = pprint ( ma.M.TL )
        minK = pprint ( self.currentMinK() )
        saving = "yes" if self.save_hiscores else "no"
            
        self.log ( f"New result with K={K}, TL={TL}, needs to pass K>{minK}, saving: {saving}" )
        if not self.save_hiscores:
            return False
        if ma.M.K == None:
            return False # clearly out
        #if ma.M.K <= self.currentMinK():
        #    return False ## clearly out
        self.addResult ( ma )
        return True

if __name__ == "__main__":
    L=[ {"K": 7.6, "x": "d"}, {"K": 7.2, "x": "e"}, {"K": 7.1, "x": "f"} ]
    for x in [ 6.9, 6.7, 6.5, 6.4, 6.3, 6.2, 6.1, 6.0, 5.9, 5.8 ]:
        L.append ( { "K": x, "x": "blah" } )
    hi={"K": 7.3, "x": "new"}
    hilist = Hiscores ( 0, False )
    print ( hilist.insertHiscore( L, hi ) )
    # hilist.updateHiscoreFile( )
