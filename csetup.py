#!/usr/bin/env python3

""" setup rundir, pythonpaths. Mostly for CLIP. """

import os, sys

def getDBPath ( dbpath, rundir ):
    """ obtain the database path, resolve <rundir> """
    if "<rundir>" in dbpath:
        dbpath = dbpath.replace("<rundir>",rundir+"/")
    dbpath = dbpath.replace("//","/")
    return dbpath

def setup( rundir = None, codedir = None ):
    """
    :param rundir: if not None, override the rundir defined per default
    :param codedir: if not None, override the codedir defined per default
    """
    # codedir = "/mnt/hephy/pheno/ww/git/"
    if codedir == None:
        codedir = "/scratch-cbe/users/wolfgan.waltenberger/git/"
    sys.path.insert(0, f"{codedir}/smodels/" )
    sys.path.insert(0, f"{codedir}/protomodels/" )
    if rundir != None:
        if not "/" in rundir[:-1]:
            rundir = f"/scratch-cbe/users/wolfgan.waltenberger/{rundir}"
        if not rundir.endswith("/"):
            rundir += "/"
        os.chdir ( rundir )
        return rundir
    home = os.environ["HOME"]
    if os.path.exists ( "./rundir.conf" ):
        with open ( "./rundir.conf" ) as f:
            rundir = f.read().strip()
            rundir = rundir.replace ( "~", home )
            os.chdir ( rundir )
        return rundir
    if os.path.exists ( f"{home}/rundir.conf" ):
        with open ( f"{home}/rundir.conf" ) as f:
            rundir = f.read().strip()
            rundir = rundir.replace ( "~", home )
            os.chdir ( rundir )
        return rundir
    cwd = __file__
    #cwd = os.getcwd()
    p1 = cwd.find("protomodels")
    if p1 > 0:
        cwd = cwd[:p1+11]
    ## print ( "cwd", cwd ) 
    return cwd
