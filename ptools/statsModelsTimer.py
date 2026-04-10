#!/usr/bin/env python3

"""
.. module:: statsModelsTimer
   :synopsis: a facility to time the CPU costs of various stats models

.. moduleauthor:: Wolfgang Waltenberger <wolfgang.waltenberger@gmail.com>
   with the help of chatgpt

"""

import inspect
import functools
import re
import time, sys

statsForStats, curTimings = {}, {}

def printTimes():
    from operator import itemgetter
    sorted_pairs = sorted(statsForStats.items(), key=itemgetter(1), reverse=True)
    for k,v in sorted_pairs:
        print ( f"{k:20} {v:.2f}s" )

def reset():
    statsForStats, curTimings = {}, {}

def start(name, owner, args, kwargs):
    """ start the timing """
    # owner is instance for instance methods, class for classmethods, 
    # or None for staticmethods
    if False:
        print(f"[start] owner={owner} name={name} args={args} kwargs={kwargs}")
    dataObj = owner.dataObject
    if type ( dataObj) == list:
        return
    anaid = dataObj.globalInfo.id
    if not anaid in statsForStats:
        statsForStats[anaid]=0.
    curTimings[anaid]=time.time()

def stop(name, owner, result=None, exc=None):
    """ stop the timing """
    if exc is not None:
        print(f"[stop] onwer={owner} name={name} raised: exc={exc!r}")
        # Optionally transform or swallow the exception; here we re-raise
        raise exc
    if False:
        print(f"[stop] owner={owner} name={name} -> result={result!r}")
    dataObj = owner.dataObject
    if type ( dataObj) == list:
        return
    anaid = dataObj.globalInfo.id
    if not anaid in curTimings:
        logger.error ( f"timer wasnt started for {anaid}" )
        sys.exit()
    dt = time.time() - curTimings[anaid]
    statsForStats[anaid]+=dt
    curTimings.pop ( anaid )
    return result  # or return a transformed value

def _wrap_function(func, name, *, kind):  # kind in {'instance','class','static'}
    is_coro = inspect.iscoroutinefunction(func)

    def _call_with_hooks(*args, **kwargs):
        owner = args[0] if kind in ('instance', 'class') else None
        call_args = args[1:] if kind in ('instance', 'class') else args
        start(name, owner, call_args, kwargs)
        try:
            result = func(*args, **kwargs)
        except Exception as exc:
            return stop(name, owner, result=None, exc=exc)
        return stop(name, owner, result=result, exc=None)

    async def _call_with_hooks_async(*args, **kwargs):
        owner = args[0] if kind in ('instance', 'class') else None
        call_args = args[1:] if kind in ('instance', 'class') else args
        start(name, owner, call_args, kwargs)
        try:
            result = await func(*args, **kwargs)
        except Exception as exc:
            return stop(name, owner, result=None, exc=exc)
        return stop(name, owner, result=result, exc=None)

    wrapped = _call_with_hooks_async if is_coro else _call_with_hooks
    return functools.update_wrapper(wrapped, func)

def _matches(selector, name, attr, cls):
    """
    Returns True if (name, attr) matches selector.
    selector can be:
      - None: no match condition (treated as 'allow')
      - str: exact name match
      - iterable of strings: name must be in it
      - regex (compiled): selector.search(name) must match
      - callable: selector(name, attr, cls) -> bool
    """
    if selector is None:
        return True
    if callable(selector):
        return bool(selector(name, attr, cls))
    if isinstance(selector, str):
        return name == selector
    if hasattr(selector, "search"):  # regex
        return bool(selector.search(name))
    try:
        return name in selector  # list/tuple/set of names
    except TypeError:
        raise TypeError("selector must be None, str, iterable of names, regex, or callable")

def instrument_class(cls, *, only=None, exclude=None, include_private=False, include_dunder=False):
    """
    Patch cls in-place so its callables are routed via start/stop.

    Parameters:
      - only: filter for which names to instrument
              (None | str | iterable[str] | regex | callable(name, attr, cls) -> bool)
      - exclude: names to skip (same forms as 'only'); applied after 'only'
      - include_private: include names starting with a single underscore
      - include_dunder: include __dunder__ names
    """
    for name, attr in list(vars(cls).items()):
        # Visibility filters
        if not include_dunder and name.startswith('__') and name.endswith('__'):
            continue
        if not include_private and name.startswith('_') and not (name.startswith('__') and name.endswith('__')):
            continue

        # Name filters
        if only is not None and not _matches(only, name, attr, cls):
            continue
        if exclude is not None and _matches(exclude, name, attr, cls):
            continue

        # Wrap based on descriptor type
        if isinstance(attr, classmethod):
            func = attr.__func__
            setattr(cls, name, classmethod(_wrap_function(func, name, kind='class')))
        elif isinstance(attr, staticmethod):
            func = attr.__func__
            setattr(cls, name, staticmethod(_wrap_function(func, name, kind='static')))
        elif inspect.isfunction(attr):
            setattr(cls, name, _wrap_function(attr, name, kind='instance'))
        elif isinstance(attr, property):
            # Only instrument if the property name passes filters
            fget = attr.fget and _wrap_function(attr.fget, f"{name}.fget", kind='instance')
            fset = attr.fset and _wrap_function(attr.fset, f"{name}.fset", kind='instance')
            fdel = attr.fdel and _wrap_function(attr.fdel, f"{name}.fdel", kind='instance')
            setattr(cls, name, property(fget, fset, fdel, attr.__doc__))
    return cls

from smodels.statistics import statsTools
instrument_class ( statsTools.StatsComputer, only = { "likelihood" } )

