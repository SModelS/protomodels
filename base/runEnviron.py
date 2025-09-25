#!/usr/bin/env python3

""" code that encapsulates the run environment,
all parameters that are specific to a run
"""

import os
from ptools.helpers import py_dump, py_dumps

__all__ = [ "RunEnviron" ]

class RunEnviron:
    """ captures all parameters that pertain to a specific 'run' 
    A 'run' is characterized a by a slew of similar walkers
    """
    def __init__ ( self, dictfile : str = "run.dict" ):
        """
        :ivar dbpath(str): the database path
        :ivar select(str): what txnames to select
        :ivar do_srcombine(bool): do sr-combinations
        :ivar forbidden(list): list of forbidden particles
        :ivar templateSLHA(str): template SLHA file
        :ivar allowN1N1Prod(bool): allow production of N1 N1
        :ivar susy_model(bool): susy mode (penalize for unusual xsecs)
        :ivar use_initialiser(Union[None,str]): use given dict file
        :ivar rundir(str): the run directory
        for initialisation, or dont use initialiser (None)
        """
        self.dictfile = dictfile
        self.didReadRunDict = False # did we get the info from run.dict?
        self._runDict = {} 
        self._runDict = self.defaults()
        self.readRunDict()
        self._setAttrs()

    @classmethod
    def create ( cls, **args ):
        """ create the run.dict, return the object """
        defaults = cls.defaults()
        import copy
        newdict = copy.deepcopy( defaults )
        dictfile = "run.dict"
        if "dictfile" in args:
            dictfile = args["dictfile"]
            args.pop("dictfile")
        newdict.update ( **args )
        with open ( dictfile, "wt" ) as f:
            py_dump ( newdict, f )
        ret = RunEnviron ( dictfile )
        return ret

    @classmethod
    def defaults ( obj ) -> dict:
        """ sets and returns the default values """
        defaults = { "dbpath": "official", "select": "all",
            "do_srcombine": True, "forbidden": [],
            "templateSLHA": "template_default.slha",
            "allowN1NProd": False, "susy_mode": False,
            "rundir": os.getcwd(),
            "use_initialiser": None }
        return defaults

    def __str__ ( self ):
        return py_dumps ( self._runDict ) 

    def readRunDict ( self ):
        """ read the run.dict file, set runDict """
        if self.dictfile is None or not os.path.exists ( self.dictfile ):
            return
        self.didReadRunDict = True
        with open ( self.dictfile, "rt" ) as f:
            txt = f.read()
            d = eval ( txt )
            self._runDict.update ( d )
        self._setAttrs()

    @property
    def templateSLHA(self):
        return os.path.join ( os.path.dirname ( __file__ ), "templates", \
                              self.templateName )

    def _setAttrs ( self ):
        for key, value in self._runDict.items():
            if key == "templateSLHA":
                key = "templateName"
            setattr ( self, key, value )

if __name__ == "__main__":
    environ = RunEnviron()
    print ( "environment at environ. Try e.g. environ.templateSLHA" )
    import sys, IPython; IPython.embed( colors = "neutral" ); sys.exit()
