""" tools, convenience functions that come in handy when dealing with user
defined priors
"""

import glob
import hashlib
import importlib.util
import inspect
import os
import sys

## the global facility: name -> {"function", "description", "module", "path"}
PRIOR_REGISTRY = {}

def registerPrior ( function, name : str|None = None,
                    description : str|None = None,
                    module : str|None = None,
                    path : str|None = None ) -> dict:
    """ register a single prior function in the global PRIOR_REGISTRY

    :param function: the callable to register
    :param name: name used as the key in the registry; defaults to the
        function's __name__
    :param description: description of the prior; defaults to the function's
        docstring
    :param module: optional module name the function came from
    :param path: optional file path the function came from
    :returns: the registry entry
    """
    if not callable ( function ):
        raise TypeError ( f"can only register callables, not {type(function)}" )
    if name is None:
        name = getattr ( function, "__name__", None )
    if name is None:
        raise ValueError ( "cannot register a callable without a name" )
    if description is None:
        description = inspect.getdoc ( function ) or ""
    entry = { "function": function,
              "description": description,
              "module": module,
              "path": path }
    PRIOR_REGISTRY [ name ] = entry
    return entry

def _loadPriorModule ( path : str, reload : bool = False ) -> None:
    """ import one *_priors.py file and register the functions it defines """
    path = os.path.abspath ( path )
    digest = hashlib.sha1 ( path.encode ( "utf-8" ) ).hexdigest ()[:16]
    modname = f"_registered_priors_{digest}"
    module = sys.modules.get ( modname )
    if module is None or reload:
        spec = importlib.util.spec_from_file_location ( modname, path )
        if spec is None or spec.loader is None:
            print ( f"[prior_tools] warning: could not load prior file {path}" )
            return
        module = importlib.util.module_from_spec ( spec )
        sys.modules [ modname ] = module
        try:
            spec.loader.exec_module ( module )
        except Exception as e:
            sys.modules.pop ( modname, None )
            print ( f"[prior_tools] warning: could not import prior file {path}: {e}" )
            return
    for name, function in inspect.getmembers ( module, inspect.isfunction ):
        ## only take functions defined in this file, not imported helpers
        if function.__module__ != modname:
            continue
        registerPrior ( function, name = name, module = modname, path = path )

def discoverPriors ( directory : str|None = None,
                     pattern : str = "*_priors.py",
                     reload : bool = False ) -> dict:
    """ import all matching prior files and register their functions

    Every public function defined in a matching file is added to the global
    PRIOR_REGISTRY, using its docstring as the description. Files that cannot
    be imported are reported and skipped, so that one broken file does not
    prevent the remaining priors from being discovered.

    :param directory: directory to search; defaults to the directory of this
        file
    :param pattern: glob pattern for the files to import
    :param reload: if True, re-execute files that were already imported
    :returns: the updated PRIOR_REGISTRY
    """
    if directory is None:
        directory = os.path.dirname ( os.path.abspath ( __file__ ) )
    directory = os.path.abspath ( os.path.expanduser ( directory ) )
    if not os.path.isdir ( directory ):
        raise NotADirectoryError ( f"prior directory {directory} does not exist" )
    ## allow the prior files to import helper modules from their directory
    if directory not in sys.path:
        sys.path.insert ( 0, directory )
    for path in sorted ( glob.glob ( os.path.join ( directory, pattern ) ) ):
        _loadPriorModule ( path, reload = reload )
    return PRIOR_REGISTRY

def printRegisteredPriors ( priors : dict|None = None ) -> None:
    """ print the names and short descriptions of registered priors

    :param priors: registry to print; defaults to the global PRIOR_REGISTRY
    """
    if priors is None:
        priors = PRIOR_REGISTRY
    for name, entry in sorted ( priors.items () ):
        description = entry.get ( "description", "" ).splitlines ()
        first_line = description [ 0 ] if description else "(no description)"
        print ( f"  {name}: {first_line}" )

def makeUserPrior ( priorConfig : dict, verbose : bool = False,
                    extraPriors : list|None = None,
                    **commonKwargs ) -> callable:
    """ build a user prior from a dictionary of registered priors

    The keys of priorConfig are names in PRIOR_REGISTRY. The values configure
    the corresponding prior:

    - a dict: keyword arguments passed to that prior
    - True: use the prior without extra keyword arguments
    - None or False: do not use this prior

    Keyword arguments given to makeUserPrior itself are passed to every
    configured prior; per-prior arguments take precedence.

    The returned function has the usual user-prior signature
    (protomodel, nll=False). For nll=True it returns the summed penalty, for
    nll=False it returns exp(-penalty). Priors may return either a single
    number or a (penalty, reason) tuple.

    :param priorConfig: dictionary mapping prior names to their configuration
    :param verbose: print the individual reasons when the prior is evaluated
    :param extraPriors: optional callable or list of callables to register
        before building the prior
    :param commonKwargs: keyword arguments passed to every configured prior
    :returns: a function taking (protomodel, nll) and returning a float
    """
    if not isinstance ( priorConfig, dict ):
        raise TypeError ( f"priorConfig has to be a dict, not {type(priorConfig)}" )
    if extraPriors is not None:
        if callable ( extraPriors ):
            extraPriors = [ extraPriors ]
        for function in extraPriors:
            registerPrior ( function )
    if not PRIOR_REGISTRY or \
            any ( name not in PRIOR_REGISTRY for name in priorConfig ):
        discoverPriors ()
    configured = []
    for name, options in priorConfig.items():
        if name not in PRIOR_REGISTRY:
            available = ", ".join ( sorted ( PRIOR_REGISTRY.keys () ) )
            raise KeyError ( f"unknown prior '{name}'; available priors: {available}" )
        if options is None or options is False:
            continue
        if options is True:
            kwargs = {}
        elif isinstance ( options, dict ):
            kwargs = dict ( options )
        else:
            raise TypeError ( f"configuration for prior '{name}' has to be a "
                              f"dict, True, None or False, not {type(options)}" )
        configured.append ( ( name, kwargs ) )

    def userPrior ( protomodel, nll : bool = False ) -> float:
        """ the configured user prior """
        reasons = []
        penalty = 0.
        for name, kwargs in configured:
            callKwargs = dict ( commonKwargs )
            callKwargs.update ( kwargs )
            function = PRIOR_REGISTRY [ name ][ "function" ]
            result = function ( protomodel, **callKwargs )
            if isinstance ( result, ( tuple, list ) ):
                if len ( result ) != 2:
                    raise ValueError ( f"prior '{name}' returned {len(result)} values, "
                                      f"expected (penalty, reason)" )
                p, reason = result
            else:
                p, reason = result, ""
            penalty += p
            if reason not in [ None, "" ]:
                mergeReasons ( reasons, reason )
        if verbose:
            print ( f"[makeUserPrior] penalty={penalty:.3g} reason:" )
            if reasons:
                print ( " - " + "\n - ".join ( str ( r ) for r in reasons ) )
            else:
                print ( " - (no reasons)" )
        userPrior.lastReasons = reasons
        return return_penalty ( penalty, nll )

    names = ", ".join ( name for name, _ in configured )
    userPrior.__doc__ = f"user-configurable prior combining: {names}"
    userPrior.priorConfig = priorConfig
    userPrior.lastReasons = []
    return userPrior

def mergeReasons ( reasons : list, more_reasons : str|list ):
    """ very simple convenience method to merge reasons """
    if type(more_reasons) == str:
        reasons.append ( more_reasons )
        return
    for r in more_reasons:
        reasons.append ( r )
    
def return_penalty ( p : float, is_nll : bool ) -> float:
    """ very simple convenience to exponentiate penalty if needed """
    import numpy as np
    if is_nll:
        return p
    return np.exp ( -p )
