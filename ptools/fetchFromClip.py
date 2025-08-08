#!/usr/bin/env python3

""" simple tool to fetch all sorts of files from clip """

import subprocess, sys, copy, argparse, os
from smodels_utils.helper.terminalcolors import *

def fetch ( files, rundir ):
    """ fetch the files in list """
    print ( f"{GREEN}fetching:{', '.join ( files )}{RESET}" )
    files = set ( files ) ## remove dupes
    basedir = os.environ['CODEDIR']
    for i in files:
        cmd=f"scp clip-login-1:{basedir}{rundir}/{i} ."
        print ( cmd )
        ret = subprocess.run(cmd.split(" "), stderr=sys.stderr, stdout=sys.stdout)

def main():
    import argparse
    argparser = argparse.ArgumentParser(
            description='file fetching utility, fetches from clip' )
    argparser.add_argument ( '-a', '--all', help='all files', action="store_true" )
    argparser.add_argument ( '--scan', help='files from scan', action="store_true" )
    argparser.add_argument ( '--hiscore', help='the hiscore', action="store_true" )
    argparser.add_argument ( '--states', help='the states', action="store_true" )
    argparser.add_argument ( '--pmodels', help='the pmodels', action="store_true" )
    argparser.add_argument ( '--png', help='the png files', action="store_true" )
    argparser.add_argument ( '--history', help='the history files', action="store_true" )
    argparser.add_argument ( '--database', help='the default.pcl database file', action="store_true" )
    argparser.add_argument ( '--file', help='the file <file>', type=str, default=None )
    argparser.add_argument ( '--dbdict', help='the database.dict file', action="store_true" )
    argparser.add_argument ( '--ssms', help='the ssm files', action="store_true" )
    argparser.add_argument ( '--llhds', help='the llhd files', action="store_true" )
    argparser.add_argument ( '--fake', help='the fake* databases', action="store_true" )
    argparser.add_argument ( '--slha', help='the hiscore.slha file', action="store_true" )
    argparser.add_argument ( '--copy', help='the copy of the hiscore file', 
                             action="store_true" )
    argparser.add_argument ( '-2', '--two', help='the second hiscore file', 
                             action="store_true" )
    argparser.add_argument ( '-R', '--rundir', help='name of remote rundir folder [rundir]', 
                             type=str, default=None )
    args = argparser.parse_args()
    if args.rundir == None:
        args.rundir = "rundir"
        if args.history:
            print ( "[fetchFromClip] did not supply a rundir but fetching history files, so choosing rundir.history" )
            args.rundir = "rundir.history"
    # files= [ "hiscore.cache" ]
    files = set()
    if args.file not in [ None, "", "None", "none" ]:
        files.add ( args.file )
    store = { "scan": [ "scanM*.pcl", "ssm*.pcl" ], 
              "states": [ "states.dict" ],
              "database": [ "default.pcl" ],
              "llhds": [ "llhd*pcl" ],
              "two": [ "hiscore2.hi" ],
              "fake": [ "fake*.pcl", "signal*.pcl" ],
              "copy": [ "hiscoreCopy.hi" ],
              "slha": [ "hiscore.slha" ],
              "dbdict": [ "d*.dict", "my.signal" ],
              "pmodels": [ "pmodel?.dict" ],
              "history": [ "history*.list" ],
              "png": [ "*.png" ],
              "hiscore": [ "hiscore.cache" ],
              "ssms": [ "ssm*.pcl" ]
    }
    for k,v in args.__dict__.items():
        if v == True or args.all:
            if k in store:
                for f in store[k]:
                    files.add ( f )
    if len(files) == 0 or args.all: ## the default
        files.add ( "hiscores.cache" )
    fetch ( files, args.rundir )

if __name__ == "__main__":
    main()
