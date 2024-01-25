#!/usr/bin/env python3

def create ( path : str ):
    """ create a rundir at <path> """
    import subprocess
    print ( "create", path )
    cmd = f"mkdir {path}"
    o = subprocess.getoutput ( cmd )
    cmd = f"ln -s ~/git/protomodels {path}/"
    o = subprocess.getoutput ( cmd )
    cmd = f"ln -s ~/git/protomodels/snippets/printSimpleHiscoreList.py {path}/"
    o = subprocess.getoutput ( cmd )
    cmd = f"ln -s ~/git/smodels-utils/clip/slurm_walk.py {path}/"
    o = subprocess.getoutput ( cmd )
    cmd = f"ln -s {path} ~/{os.path.basename{path}"
    o = subprocess.getoutput ( cmd )
    cmd = f"ln -s ~/git/protomodels/ptools/hiscoreCLI.py {path}/"
    o = subprocess.getoutput ( cmd )

if __name__ == "__main__":
    import sys
    if len ( sys.argv ) < 2:
        print ( "[createRundir] supply path" )
        sys.exit()
    create ( sys.argv[1] )
