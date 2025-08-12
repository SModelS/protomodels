#!/usr/bin/env python3

""" code that implements a file locking mechanism
"""

import os, signal, subprocess

ignore_locks = False
__locks__ = set()

def signal_handler(sig, frame):
    print('You pressed Ctrl+C, remove all locks!')
    for l in __locks__:
        cmd = "rm -f %s" % l
        subprocess.getoutput ( cmd )
        print ( cmd )
    sys.exit(0)

signal.signal(signal.SIGINT, signal_handler)

def lockfile ( basefile : os.PathLike ) -> os.PathLike:
    lock_file = "."+os.path.basename(basefile)+".lock"
    return lock_file

def lock ( filename : os.PathLike ) -> bool:
    """ lock the file filename, to make sure processes dont
    overwrite each other
    :returns: True if there is already a lock on it
    """
    import time, socket, random
    if ignore_locks:
        return False
    lock_file = lockfile ( filename )

    __locks__.add ( lock_file )
    if os.path.exists ( lock_file ):
        return True
    for i in range(5):
        try:
            with open ( lock_file, "wt" ) as f:
                f.write ( time.asctime()+","+socket.gethostname()+"\n" )
                f.close()
            return False
        except FileNotFoundError as e:
            t0 = random.uniform(2.,4.*i)
            print ( f"[locker] FileNotFoundError #{i} {e}. Sleep for {t0:.1f}s" )
            time.sleep( t0 )
    return True ## pretend there is a lock

def unlock ( filename : os.PathLike ) -> bool:
    """ unlock for topo and masses, to make sure processes dont
        overwrite each other """
    if ignore_locks:
        return
    lock_file = lockfile( filename )
    if lock_file in __locks__:
        __locks__.remove ( lock_file )
    if os.path.exists ( lock_file ):
        cmd = f"rm -f {lock_file}"
        subprocess.getoutput ( cmd )
