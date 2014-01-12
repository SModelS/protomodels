"""
.. module:: SModelS
    :synopsis: Intended as a potential main entry point, currently just for
               returning the SModelS version number.

.. moduleauthor:: Wolfgang Waltenberger <wolfgang.waltenberger@gmail.com>

"""

def installdir():
  """ return the software installation directory, by looking at location of this method """
  import os, inspect
  return os.path.dirname ( inspect.getabsfile(installdir) )

def version( astuple=False ):
  """ prints out version number of SModelS framework """
  f=open("%s/version" % installdir() )
  l=f.readline()
  f.close()
  if not astuple: return l
  T,C=l.split("/")
  A,B=T.split(".")
  return (int(A),int(B),C.strip())
  
