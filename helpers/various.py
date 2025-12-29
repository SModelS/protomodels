#!/usr/bin/env python3

""" helpers and convenience functions of all sorts. """

def viewImage ( imgfile : str ):
    """ view an image on the terminal, if the timg tool exists. """
    from smodels_utils.plotting.mpkitty import timg
    timg ( imgfile )
    """
    import shutil
    if shutil.which ("timg") is not None:
        import subprocess
        o = subprocess.getoutput ( f"timg {imgfile}" )
        print ( o )
    """
