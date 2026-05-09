#!/usr/bin/env python3

""" code that encapsulates the run environment,
all parameters that are specific to a run
"""

import os
import time
from ptools.helpers import py_dump, py_dumps
from base.loggerbase import LoggerBase
from base.pbase import openWithRetry

__all__ = [ "RunEnviron" ]

def dict_diff(d1, d2):
    """
    Return keys and values where two dictionaries differ.
    """
    diff = {}
    all_keys = set(d1) | set(d2)  # union of keys
    for k in all_keys:
        v1 = d1.get(k, "<missing>")
        v2 = d2.get(k, "<missing>")
        if v1 != v2:
            diff[k] = (v1, v2)
    return diff

class RunEnviron ( LoggerBase ):
    """ captures all parameters that pertain to a specific 'run'
    A 'run' is characterized a by a slew of similar walkers
    """
    def __init__ ( self, runDictFile : str = "run.dict" ):
        """
        :ivar dbpath(str): the database path
        :ivar dbver(str): the database version
        :ivar select(str): what txnames to select
        :ivar do_srcombine(bool): do sr-combinations
        :ivar forbiddenparticles(list): list of forbidden particles
        :ivar templateSLHA(str): template SLHA file
        :ivar allowN1N1Prod(bool): allow production of N1 N1
        :ivar susy_mode(bool): susy mode (penalize for unusual xsecs)
        :ivar use_initialiser(Union[None,str]): use given dict file
        for initialisation, or dont use initialiser (None)
        :ivar rundir(str): the run directory
        :ivar extrapolation_acceptance(float): the extrapolation acceptance
        in SModelS txnameObj
        """
        super ( RunEnviron, self ).__init__ ( "env" )
        self._db = None # we instantiate only when needed
        # at least for now, lets make this the default
        self.runDictFile = os.path.expanduser ( runDictFile )
        self.didReadRunDict = False # did we get the info from run.dict?
        self.run_dict = self.defaults()
        self.readRunDict()
        self.setAcceptedExtrapolationErrors ( 0. )
        self._setAttrs()


    def setAcceptedExtrapolationErrors ( self, err : float ):
        """ set the extrapolation errors of all results in the
        database to err
        :param err: default is 0.05
        """
        assert 0. <= err, "error {err} makes no sense"
        for er in self.database.expResultList:
            for ds in er.datasets:
                for txn in ds.txnameList:
                    txn.txnameData._accept_errors_upto = err
                    if txn.txnameDataExp is not None:
                        txn.txnameDataExp._accept_errors_upto = err
        self.extrapolation_acceptance = err

    @classmethod
    def create ( cls, **args ):
        """ create the run.dict, return the object """
        defaults = cls.defaults()
        import copy
        newdict = copy.deepcopy( defaults )
        runDictFile = "run.dict"
        check = True
        if "check" in args:
            check = args["check"]
            args.pop("check")
        if "runDictFile" in args:
            runDictFile = os.path.expanduser ( args["runDictFile"] )
            args.pop("runDictFile")
        oldret = None
        if os.path.exists ( runDictFile ):
            oldret = RunEnviron ( runDictFile )
        newdict.update ( **args )
        if "rundir" in newdict:
            path = os.path.abspath ( os.path.expanduser ( newdict["rundir"] ) )
            newdict["rundir"] = path
        old_dbver = "???"
        if oldret != None:
            ## dbver we check later
            old_dbver = oldret.run_dict["dbversion"]
            oldret.run_dict.pop ( "dbversion" )
            newdict.pop ( "dbversion" )
            if check and newdict != oldret.run_dict:
                print ( f"[RunEnviron] {runDictFile} differs from previous version:" )
                dd = dict_diff(oldret.run_dict,newdict)
                for k,v in dd.items():
                    line = f"[RunEnviron] {k:>16}: {v[0]} != {v[1]}"
                    print ( line )
                    if k not in [ "use_initialiser" ]:
                        print ( f"[RunEnviron] correct this!" )
                        raise Exception ( line )
                else:
                    print ( f"[RunEnviron] thats allowed!" )
        #    else: # newdict is compatible with old dict
        #        return oldret
        if not "dbversion" in newdict:
            newdict["dbversion"]=old_dbver
        py_dump ( newdict, runDictFile )
        ret = RunEnviron ( runDictFile )
        new_dbver = ret.databaseVersion
        if check and old_dbver != new_dbver:
            line= f"[RunEnviron] dbver changed from {old_dbver} to {new_dbver}"
            print ( line )
            raise Exception ( line )
        return ret

    @classmethod
    def new ( cls, **args ):
        """ create a new run.dict, override any old """
        runDictFile = "run.dict"
        if os.path.exists ( runDictFile ):
            import shutil
            shutil.move ( runDictFile, "run_old.dict" )
        args["check"]=False
        return cls.create ( **args )

    def moveRunDict ( self, dest : os.PathLike ):
        """ move this run.dict file """
        import shutil
        shutil.move ( self.runDictFile, dest )
        self.runDictFile = dest

    @classmethod
    def defaults ( obj ) -> dict:
        """ sets and returns the default values """
        defaults = { "dbpath": "official", "select": "all",
            "do_srcombine": True, "forbiddenparticles": [],
            "templateSLHA": "template_default.slha",
            "allowN1N1Prod": False, "susy_mode": False,
            "rundir": os.getcwd(), "strategy": "aggressive",
            "use_initialiser": None, "dbversion": "???",
            "extrapolation_acceptance": 0.
        }
        return defaults

    def __str__ ( self ):
        return py_dumps ( self.run_dict )

    def readRunDict ( self ):
        """ read the run.dict file, set runDict """
        if self.runDictFile is None or not os.path.exists ( self.runDictFile ):
            self.error ( f"no {self.runDictFile} exists" )
            self.error ( f"you can create one via RunEnviron.create()" )
            import sys; sys.exit()
        txt = ""
        try:
            while txt == "":
                with openWithRetry ( self.runDictFile, "rt" ) as f:
                    txt = f.read()
            d = eval ( txt )
            self.run_dict.update ( d )
            self.didReadRunDict = True
        except ( SyntaxError, ValueError ) as e:
            self.error (f"when parsing {os.path.abspath(self.runDictFile)} {type(e)}: {e}" )
            self.error (f"content: >>{txt}<<" )
            import sys; sys.exit(-1)
        self._setAttrs()

    @property
    def database(self):
        if self._db != None:
            return self._db
        from tester.combinationsmatrix import getYamlMatrix
        combinationsmatrix, status = getYamlMatrix()
        if not combinationsmatrix or status != 0:
            sys.exit("Combination matrix not loaded correctly in RunEnviron.")
        force_load = None
        if self.dbpath.endswith ( ".pcl" ):
            force_load = "pcl"
        if "/" in self.dbpath:
            ntries = 0
            while not os.path.exists ( self.dbpath ):
                ## give it a few tries
                ntries += 1
                time.sleep ( ntries * 5 )
                if ntries > 5:
                    break
        from smodels.experiment.databaseObj import Database
        self._db = Database ( self.dbpath, force_load = force_load,
               combinationsmatrix = combinationsmatrix )
        if "official" not in self.dbpath:
            from smodels_utils.helper.databaseManipulations import removeNonAggregatedFromDB
            self._db = removeNonAggregatedFromDB( self._db )
        return self._db

    @property
    def databaseVersion(self):
        return self.database.databaseVersion

    @property
    def templateSLHA(self):
        fname = os.path.dirname ( __file__ )
        ret = os.path.join ( fname, "..", "builder" , "templates", self.templateName )
        return os.path.abspath ( ret )

    def __eq__ ( self, other ):
        if type(other) != RunEnviron:
            return False
        return self.run_dict == other.run_dict

    def _setAttrs ( self ):
        for key, value in self.run_dict.items():
            if key == "templateSLHA":
                key = "templateName"
            if key == "rundir":
                value = os.path.expanduser ( value )
            setattr ( self, key, value )

if __name__ == "__main__":
    if os.path.exists ( "run.dict" ):
        environ = RunEnviron.create()
    else:
        environ = RunEnviron()
    #print ( "environment at environ. Try e.g. environ.templateSLHA" )
    import sys, IPython; IPython.embed( colors = "neutral" ); sys.exit()
