#!/usr/bin/env python3

""" take all models in Pmodels, agglomerate them into a simple list """


def agglomerate():
    import glob
    files = glob.glob ( "Pmodels/pmodel*.dict" )
    models = []
    for fl in files:
        with open ( fl, "rt" ) as f:
            m = eval ( f.read() )
            models.append ( m )
    from ptools.helpers import py_dump
    py_dump ( models, "pmodels.list" )


if __name__ == "__main__":
    agglomerate()
