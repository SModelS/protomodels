#!/usr/bin/env python3

""" create a rundir environment with all conveniences """

import os

def create ( path : str ):
    """ create a rundir at <path> """
    import subprocess
    print ( "create", path )
    if not os.path.exists ( path ):
        cmd = f"mkdir {path}"
        o = subprocess.getoutput ( cmd )
    if not os.path.exists ( f"{path}/protomodels" ):
        cmd = f"ln -s ~/git/protomodels {path}/"
        o = subprocess.getoutput ( cmd )
    if not os.path.exists ( f"{path}/printSimpleHiscoreList.py" ):
        cmd = f"ln -s ~/git/protomodels/snippets/printSimpleHiscoreList.py {path}/"
        o = subprocess.getoutput ( cmd )
    if not os.path.exists ( f"{path}/slurm_walk.py" ):
        cmd = f"ln -s ~/git/smodels-utils/clip/slurm_walk.py {path}/"
        o = subprocess.getoutput ( cmd )
    basedirpath = os.path.expanduser ( f"~/{os.path.basename(path)}" )
    if not os.path.exists ( basedirpath ):
        cmd = f"ln -s {path} {basedirpath}"
        o = subprocess.getoutput ( cmd )
    if not os.path.exists ( f"{path}/hiscoreCLI.py" ):
        cmd = f"ln -s ~/git/protomodels/ptools/hiscoreCLI.py {path}/"
        o = subprocess.getoutput ( cmd )
    if not os.path.exists ( f"{path}/clean.sh" ):
        cmd = f"ln -s ~/git/protomodels/snippets/clean.sh {path}/"
        o = subprocess.getoutput ( cmd )

if __name__ == "__main__":
    import sys
    if len ( sys.argv ) < 2:
        print ( "[createRundir] supply path" )
        sys.exit()
    create ( sys.argv[1] )
