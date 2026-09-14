#!/usr/bin/env python3

""" protomodels base, it's all about the base, no treble
"""

from typing import IO, Optional
import time
import os

def prettyFileName ( path : os.PathLike ) -> os.PathLike:
    """ /home/walten/git/ -> ~/git/ """
    ret = path
    if os.getcwd() in ret:
        ret = ret.replace( os.getcwd(), "." )
    if os.environ["HOME"] in ret:
        ret = ret.replace( os.environ["HOME"], "~" )
    return ret

def openWithRetry( path: str, mode: str = "r", retries: int = 5,
                     delay: float = 1.0,) -> IO:
    """
    Open a file, retrying on transient I/O errors.

    This is useful for files on network filesystems (e.g., NFS, SMB)
    where `open()` may fail intermittently.

    Args:
        path: Path to the file to open.
        mode: File mode (same as built-in `open`).
        retries: Number of attempts before giving up.
        delay: Seconds to sleep between retries.

    Returns:
        An open file object.

    Raises:
        Exception: If all retries fail.
    """
    last_exc: Optional[Exception] = None

    for attempt in range(retries):
        try:
            return open(path, mode)
        except Exception as exc:
            last_exc = exc
            if attempt < retries - 1:
                time.sleep(delay*(attempt**2+1))
            else:
                raise last_exc

