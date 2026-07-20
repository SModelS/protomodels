#!/usr/bin/env python3

"""File locking mechanism for multiprocess safety.

Provides file-based locking using lockfiles to prevent concurrent
processes from corrupting shared state (hiscores, databases, etc.).
"""

import os
import signal
import socket
import subprocess
import sys
import time
import random
from typing import Set

ignore_locks: bool = False
__locks__: Set[str] = set()

old_handler = signal.getsignal(signal.SIGINT)


def signal_handler(sig: int, frame) -> None:
    """Clean up all acquired locks on SIGINT and exit gracefully."""
    if sig == signal.SIGINT:
        print(
            f"You pressed Ctrl+C (sig {sig}), remove all {len(__locks__)} locks!"
        )
    for lock_path in __locks__:
        cmd = f"rm -f {lock_path}"
        subprocess.getoutput(cmd)
        print(cmd)
    if old_handler is signal.SIG_DFL:
        print("[locker] Exiting gracefully...")
        sys.exit(130)
    elif old_handler is signal.SIG_IGN:
        print("[locker] Old handler ignored SIGINT, continuing.")
    else:
        old_handler(sig, frame)


signal.signal(signal.SIGINT, signal_handler)


def lockfile(basefile: os.PathLike) -> str:
    """Return the lockfile path corresponding to *basefile*.

    :param basefile: The original file to protect.
    :returns: Path of the lock file (hidden, in the same directory).
    """
    return "." + os.path.basename(basefile) + ".lock"


def lock(filename: os.PathLike) -> bool:
    """Acquire a file lock for *filename*.

    Blocks until the lock is acquired or forced after repeated failures.
    Uses an exponential back-off when waiting for another process.

    :param filename: Path of the file to lock.
    :returns: ``True`` if the lock was acquired (or faked after timeout).
    """
    if ignore_locks:
        return False
    if not os.path.exists(filename):
        return False

    lf = lockfile(filename)

    ctr = 0
    if os.path.exists(lf):
        while os.path.exists(lf):
            time.sleep(0.5 * ctr + 0.2)
            ctr += 1
            if ctr > 6:
                unlock(filename)

    for i in range(5):
        try:
            with open(lf, "wt") as f:
                f.write(
                    f"{{ 'time': '{time.asctime()}', "
                    f"'host': '{socket.gethostname()}', "
                    f"'t': {time.time()} }}\n"
                )
            __locks__.add(lf)
            return True
        except FileNotFoundError:
            t0 = random.uniform(2.0, 4.0 * i)
            print(f"[locker] FileNotFoundError #{i}. Sleep for {t0:.1f}s")
            time.sleep(t0)

    __locks__.add(lf)
    return True  # pretend there is a lock


def unlock(filename: os.PathLike) -> bool:
    """Release the file lock for *filename*.

    :param filename: Path of the file whose lock should be released.
    :returns: ``True`` if there was a lock that got removed.
    """
    lf = lockfile(filename)
    __locks__.discard(lf)
    if os.path.exists(lf):
        try:
            os.unlink(lf)
            return True
        except FileNotFoundError:
            pass
    return False
