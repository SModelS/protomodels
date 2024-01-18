#!/usr/bin/env python3

"""
.. module:: installation
   :synopsis: a module for returning installation paths and version numbers.

.. moduleauthor:: Wolfgang Waltenberger <wolfgang.waltenberger@gmail.com>

"""

from typing import Union, Tuple

def version( return_tuple : bool = False) -> Union[str,Tuple]:
    """
    Return version number of the protomodels code

    :param return_tuple: Return a tuple of (major,minor,patch,...) instead of string
    :returns: version, either as string or as tuple
    """
    f = open("version" )
    l = f.readline()
    f.close()
    l = l.replace("\n", "")
    l.strip()
    if not return_tuple:
        return l
    import re
    ret = re.split("\.|-",l)
    for i,r in enumerate(ret):
        try:
            ret[i]=int(r)
        except ValueError as e:
            pass
    return tuple(ret)

def smodels_version( **args ):
    """ Return SModelS version """
    from smodels.installation import version
    return version ( **args )

