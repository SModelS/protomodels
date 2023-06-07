#!/usr/bin/env python3

""" A class that centralizes access to the hiscore list over multiple threads.
"""

__all__ = [ "Hiscore" ]

import random, copy, pickle, os, fcntl, time, subprocess, colorama
from scipy import stats
from builder.manipulator import Manipulator
from tester.combiner import  Combiner
from ptools import helpers
from ptools.csetup import setup
from ptools import sparticleNames
from typing import Union
from os import PathLike

class Hiscore:
    """ encapsulates the hiscore list. """
    def __init__ ( self, walkerid: int = 0, save_hiscores: bool = False,
                   picklefile: PathLike="hiscore.hi", backup : bool = True, 
                   hiscores = None, predictor = None ):
        """ the constructor
        :param save_hiscores: if true, then assume you want to save, not just read.
        :param picklefile: path of pickle file name to connect hiscore list with
        :param backup: if True, make a backup pickle file old_<name>.pcl
        :param hiscores: if None, try to get them from file, if a list,
                         then these are the hiscore protomodels.
        """
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
        if hiscores == None:
            self.updateListFromPickle ( )
        else:
            self.hiscores = hiscores
            self.mtime = time.time()

    def currentMinZ ( self ):
        """ the current minimum Z to make it into the list. """
        if self.hiscores[-1] == None:
            return 0.
        return self.hiscores[-1].Z

    def currentMinK ( self, zeroIsMin=False ):
        """ the current minimum K to make it into the list.
        :param zeroIsMin:  if false, min k can become negative
        """
        if self.hiscores[-1] == None:
            if zeroIsMin:
                return 0.
            return -30.
        mk = -10.
        if hasattr ( self.hiscores[-1], "K" ):
            mk = self.hiscores[-1].K
        if zeroIsMin:
            return max ( mk, 0. )
        return mk

    def globalMaxZ ( self ):
        """ globally (across all walkers), the highest Z """
        ret = 0.
        if self.hiscores[0] != None:
            if self.hiscores[0].Z > ret:
                ret = self.hiscores[0].Z
        Zoldfile = "Zold.conf"
        if os.path.exists ( Zoldfile ):
            with open ( Zoldfile, "rt" ) as f:
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

    def similarDicts ( self, a : dict, b : dict ):
        """ are models a and b similar? """
        dK = abs ( a["K"] - b["K"] )
        if dK > 1e-5:
            return False
        dZ = abs ( a["Z"] - b["Z"] )
        if dZ > 1e-5:
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
                return L
            if oldhi["K"]>= K:
                ret.append ( oldhi )
            else:
                break
        ## now the new hiscore
        ret.append ( hi )
        ## now fill up
        for oldhi in L[len(ret)-1:]:
            ret.append ( oldhi )
        ret = ret[:10] ## cut off, max ten
        return ret

    @classmethod
    def fromDictionaryFile ( cls, path : PathLike ):
        assert False, "implement"

    def writeToHiscoreFile ( self, m : Manipulator ):
        """ we have a new hiscore, write to hiscores.dict
        :param m: manipulator
        """
        oldhiscores=[]
        fname = "hiscores.dict"
        if os.path.exists ( fname ):
            tryRead=0
            success=False
            ## stop loop at success or when tryRead is at least 5
            while (not success) and tryRead<5:
                tryRead+=1
                try:
                    with open ( fname, "rt" ) as h:
                        txt = h.read()
                        oldhiscores = eval( txt )
                        h.close()
                        success=True
                except SyntaxError as e:
                    time.sleep( .1+3*tryRead )
        D=m.M.dict() ## FIXME use writeDictFile instead!
        D["K"]=m.M.K
        D["Z"]=m.M.Z
        if hasattr ( m, "seed" ) and m.seed != None:
            D["seed"]=m.seed
        D["step"]=m.M.step
        D["timestamp"]=time.asctime()
        D["walkerid"]=m.M.walkerid
        # D=m.writeDictFile(outfile = None, ndecimals=6 )
        newlist = self.insertHiscore ( oldhiscores, D )
        self.pprint ( f"write model to {fname}" )
        with open ( fname, "wt" ) as f:
            f.write ( "[" )
            for ctr,l in enumerate(newlist):
                f.write ( "%s" % l )
                if ctr < len(newlist)-1:
                    f.write ( ",\n" % ( l ) )
            f.write ( "]\n" )
            f.close()
        with open ( "Kold.conf", "wt" ) as f:
            f.write ( "%f\n" % m.M.K  )
            f.close()
        with open ( "Kmin.conf", "wt" ) as f:
            f.write ( "%f\n" % newlist[-1]["K"] )
            f.close()
        return True

    def addResult ( self, ma ):
        """ add a result to the list
        :param ma: the manipulator object
        :returns: true, if result was added
        """
        if ma.M.K <= self.currentMinK( zeroIsMin = True ):
            return False ## doesnt pass minimum requirement
        if ma.M.K == 0.:
            return False ## just to be sure, should be taken care of above, though

        # Kold = self.globalMaxK()
        Kmin = self.globalMinK()
        # self.pprint ( f"adding results Kold is {Kold} Knew is {ma.M.K}" )
        ## FIXME we should only write into this file in the first maxstep/3 steps
        if ma.M.K > Kmin:
            self.pprint ( "WARNING we shouldnt write into hiscore file afte maxstep/3 steps!!" )
            self.writeToHiscoreFile( ma )
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
        from smodels.tools import runtime
        runtime._experimental = True

        #Make sure the model is backep up
        manipulator.backupModel()

        unfrozen = manipulator.M.unFrozenParticles( withLSP=False )
        oldZ = manipulator.M.Z
        oldK = manipulator.M.K
        particleContributions = {} ## save the scores for the non-discarded particles.
        #particleContributionsZ = {} ## save the scores for the non-discarded particles, Zs

        #Make sure predictor is accesible
        if not self.predictor:
            self.pprint( "asked to compute particle contributions to score, but predictor has not been set")
            return

        pidsnmasses = [ (x,manipulator.M.masses[x]) for x in unfrozen ]
        pidsnmasses.sort ( key=lambda x: x[1], reverse=True )
        for cpid,(pid,mass) in enumerate(pidsnmasses):
            self.pprint ( "computing contribution of %s (%.1f): [%d/%d]" % \
                   ( self.namer.asciiName(pid),
                     manipulator.M.masses[pid],(cpid+1),len(unfrozen) ) )

            #Remove particle and recompute SLHA file:
            manipulator.freezeParticle(pid, force=True )
            #Recompute cross-secions:
            manipulator.M.getXsecs()
            manipulator.M.K = 0.0
            manipulator.M.Z = 0.0
            self.predictor.predict( manipulator.M )
            if manipulator.M.K is None:
                self.pprint ( "when removing %s, K could not longer be computed. Setting to zero"% ( self.namer.asciiName(pid)))
                manipulator.M.K = 0.0
                manipulator.M.Z = 0.0
            if oldK <= 0:
                percK = 0.
            else:
                percK = ( manipulator.M.K - oldK ) / oldK
                self.pprint ( "when removing %s, K changed: %.3f -> %.3f (%.1f%s), Z: %.3f -> %.3f (%d evts)" % \
                    ( self.namer.asciiName(pid), oldK, manipulator.M.K, 100.*percK, "%", oldZ,manipulator.M.Z, manipulator.M.nevents ) )

            #Store the new Z and K values in the original model:
            particleContributions[pid]=manipulator.M.K
            #particleContributionsZ[pid]=manipulator.M.Z
            #Make sure to restore the model to its initial (full particle content) state
            manipulator.restoreModel()
            #Store contributions in the protomodel:
            manipulator.M.particleContributions = particleContributions
            #manipulator.M.particleContributionsZ = particleContributionsZ

        self.pprint ( "stored %d particle contributions" % len(manipulator.M.particleContributions) )

    def computeAnalysisContributions( self, manipulator ):
        """ compute the contributions to Z of the individual analyses
        :returns: the model with the analysic constributions attached as
                  .analysisContributions
        """

        try:
            self.pprint ( "Now computing analysis contributions" )
            self.pprint ( "Recompute the score. Old one at K=%.2f, Z=%.2f" % \
                          ( manipulator.M.K, manipulator.M.Z ) )
            contributionsZ = {}
            contributionsK = {}
            combiner = Combiner()
            dZtot, dKtot = 0., 0.
            bestCombo = copy.deepcopy ( manipulator.M.bestCombo )
            #self.pprint ( "we have %d entries in best combo" % len(bestCombo) )
            prior = combiner.computePrior ( manipulator.M )
            #self.pprint ( "the prior is %s" % prior )
            for ctr,pred in enumerate(bestCombo):
                #self.pprint ( "Now starting to compute for %d" % ctr )
                combo = bestCombo[:ctr]+bestCombo[ctr+1:]
                # combo = copy.deepcopy ( bestCombo )[:ctr]+copy.deepcopy ( bestCombo)[ctr+1:]
                #self.pprint ( "deep copy still worked: %d" % (len(combo)) )
                Z, muhat_ = combiner.getSignificance ( combo )
                #self.pprint ( "Z for %d is %s" % ( ctr, Z ) )
                K = combiner.computeK ( Z, prior )
                #self.pprint ( "K for %d is %s" % ( ctr, K ) )
                contributionsK [ ctr ] = K
            self.pprint ( "finished computing contributions" )

            contrsWithNames = {}
            for k,v in contributionsK.items():
                # self.pprint ( "contributionsK of %s reads %s" % ( k, v ) )
                contrsWithNames [ manipulator.M.bestCombo[k].analysisId() ] = v
            manipulator.M.analysisContributions = contrsWithNames
            self.pprint ( "stored %d analyses contributions" % len(manipulator.M.analysisContributions) )
        except Exception as e:
            self.pprint ( "in computeAnalysisContributions caught %s" % str(e) )

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
            self.pprint ( "loaded %d hiscores from %s." % \
                          ( nhs, self.pickleFile ) )
            # assert ( len(self.hiscores) == self.nkeep )
            self.fileAttempts=0
        except Exception as e:
        # except OSError or BlockingIOError or EOFError or pickle.UnpicklingError or TypeError as e:
            self.fileAttempts+=1
            if self.fileAttempts<20: # try again
                self.pprint ( "Exception[X] %s: type(%s), Waiting for %s file, %d" % (str(e),type(e),self.pickleFile,self.fileAttempts) )
                time.sleep ( (.2 + random.uniform(0.,1.))*self.fileAttempts )
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

    def writeListToDictFile ( self, dictFile=None ):
        """ write the models in append mode in a single dictFile.
        :param dictFile: write to dictFile. If None, then self.pickleFile
                         is used, but with ".dict" as extension.
        """
        if dictFile==None:
            dictFile = self.pickleFile
        if dictFile.endswith(".pcl"):
            dictFile = dictFile[:-4]+".py"
        f=open(dictFile,"wt")
        f.write("[")
        f.close()
        for protomodel in self.hiscores:
            ma = Manipulator ( protomodel )
            ma.writeDictFile ( outfile = dictFile, appendMode=True )
        f=open(dictFile,"at")
        f.write("]\n")
        f.close()

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
        self.pprint ( "saving new hiscore list to %s" % pickleFile )
        try:
            if self.backup:
                subprocess.getoutput ( "mv -f %s old_%s" % ( pickleFile, pickleFile ) )
            # self.clean()
            with open( pickleFile, "wb" ) as f:
                fcntl.flock ( f, fcntl.LOCK_EX )
                pickle.dump ( self.hiscores, f )
                pickle.dump ( time.asctime(), f )
                fcntl.flock ( f, fcntl.LOCK_UN )
                f.close()
            self.mtime = os.stat ( self.pickleFile ).st_mtime
            self.fileAttempts=0
            return True
        except OSError or BlockingIOError as e:
            self.fileAttempts+=1
            if self.fileAttempts>2:
                self.pprint ( f"error when writing ({self.fileAttempts}) pickle file: {e}" )
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
        self.pprint ( "New result with K=%.2f, Z=%.2f, needs to pass K>%.2f, saving: %s" % \
                ( ma.M.K, ma.M.Z, self.currentMinK(),
                  "yes" if self.save_hiscores else "no" ) )
        if not self.save_hiscores:
            return False
        if ma.M.K <= self.currentMinK():
            return False ## clearly out
        self.addResult ( ma )
        # self.writeListToPickle() ## and write it to pickle
        return True

    def pprint ( self, *args ):
        """ logging """
        print ( "[hiscore:%d] %s" % ( self.walkerid, " ".join(map(str,args))) )
        self.log ( *args )

    def log ( self, *args ):
        """ logging to file """
        # logfile = "walker%d.log" % self.walkerid
        logfile = "walker%d.log" % self.walkerid
        with open( logfile, "at" ) as f:
            tm = time.strftime("%b %d %H:%M:%S")
            f.write ( "[hiscore-%s] %s\n" % ( tm, " ".join(map(str,args)) ) )

if __name__ == "__main__":
    L=[ {"K": 7.6, "x": "d"}, {"K": 7.2, "x": "e"}, {"K": 7.1, "x": "f"} ]
    for x in [ 6.9, 6.7, 6.5, 6.4, 6.3, 6.2, 6.1, 6.0, 5.9, 5.8 ]:
        L.append ( { "K": x, "x": "blah" } )
    hi={"K": 7.3, "x": "new"}
    hilist = Hiscore ( 0, False )
    print ( hilist.insertHiscore( L, hi ) )
    hilist.writeToHiscoreFile( )
