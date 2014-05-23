#!/usr/bin/python

"""
.. module:: setPath
   :synopsis: Sets the path such that e.g. from smodels.tools import toolBox works.
              correctly. Called as a script, the path is printed.

.. moduleauthor:: Wolfgang Waltenberger <wolfgang.waltenberger@gmail.com>

"""

from __future__ import print_function
import sys
import inspect
import os

def configure():
    """
    Get the path name of this file, remove set_path.py, remove the last
    subdir.

    The remaining string should be the base path name.

    """
    base = os.path.dirname(os.path.realpath(inspect.getabsfile(configure)))
    pos = base.rfind("/smodels")
    base = base[:pos + 1]
    sys.path.append(base)
    sys.path.append(base[:-9])
    return base


configure()


if __name__ == "__main__":
    """
    Called as a script, print out the path.

    """
    print("The following string is appended to the path variable:",
          configure())
