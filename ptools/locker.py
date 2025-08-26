#!/usr/bin/env python3

""" code that implements a file locking mechanism
"""

import os, signal, subprocess

ignore_locks = False
__locks__ = set()
    
old_handler = signal.getsignal(signal.SIGINT)

def signal_handler(sig, frame):
    if sig == signal.SIGINT:
        print( f'You pressed Ctrl+C, remove all locks! {sig}')
    for l in __locks__: ## remove always
        cmd = f"rm -f {l}"
        subprocess.getoutput ( cmd )
        print ( cmd )
    # os.kill ( os.getpid(), sig )
    if old_handler is signal.SIG_DFL:
        # Default behavior for SIGINT is to raise KeyboardInterrupt,
        # which usually exits with code 130.
        print("[locker] Exiting gracefully...")
        sys.exit(130)
    elif old_handler is signal.SIG_IGN:
        print("[locker] Old handler ignored SIGINT, continuing.")
    else:
        # Call previous custom handler
        old_handler(signum, frame)

signal.signal(signal.SIGINT, signal_handler)

def lockfile ( basefile : os.PathLike ) -> os.PathLike:
    lock_file = "."+os.path.basename(basefile)+".lock"
    return lock_file

def lock ( filename : os.PathLike ) -> bool:
    """ lock the file filename, to make sure processes dont
    overwrite each other

    :returns: True if it was able to lock
    """
    import time, socket, random
    if ignore_locks:
        return False
    if not os.path.exists ( filename ):
        # dont lock non-existing file
        return False
    lock_file = lockfile ( filename )

    ## a lock file exists already? wait!
    ctr = 0
    if os.path.exists ( lock_file ):
        while ( os.path.exists ( lock_file ) ):
            time.sleep ( .5*ctr + .2 )
            ctr += 1
            if ctr > 6: # we force an unlock after some time
                unlock ( filename )
    for i in range(5):
        try:
            with open ( lock_file, "wt" ) as f:
                f.write ( f"{{ 'time': '{time.asctime()}', 'host': '{socket.gethostname()}', 't': {time.time()} }}\n" )
                f.close()
            __locks__.add ( lock_file )
            return True
        except FileNotFoundError as e:
            t0 = random.uniform(2.,4.*i)
            print ( f"[locker] FileNotFoundError #{i} {e}. Sleep for {t0:.1f}s" )
            time.sleep( t0 )
    __locks__.add ( lock_file )
    return True ## pretend there is a lock

def unlock ( filename : os.PathLike ) -> bool:
    """ unlock filename, to make sure processes dont
        overwrite each other 

    :returns: true if there really was a lock 
    """
    #if ignore_locks:
    #    return
    lock_file = lockfile( filename )
    if lock_file in __locks__:
        __locks__.remove ( lock_file )
    if os.path.exists ( lock_file ):
        try:
            os.unlink ( lock_file )
            return True
        except FileNotFoundError as e:
            pass
        #cmd = f"rm -f {lock_file}"
        #subprocess.getoutput ( cmd )
    return False
