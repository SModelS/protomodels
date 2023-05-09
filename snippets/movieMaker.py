#!/usr/bin/env python3

"""
.. module:: movieMaker
   :synopsis: a simple movie maker, makes a video clip that shows the evolution of
              the meta statistics over time

.. moduleauthor:: Wolfgang Waltenberger <wolfgang.waltenberger@gmail.com>

"""

import subprocess, os, sys

class MovieMaker:
    def __init__ ( self, dictfile : str ):
        self.dirname = "pics"
        self.dictfile = dictfile
        self.outfile = "out.mp4"
        self.checkDictFile()

    def mkdir ( self ):
        if os.path.exists ( self.dirname ):
            cmd = f"rm -r {self.dirname}"
            subprocess.getoutput ( cmd )
        cmd = f"mkdir {self.dirname}"
        subprocess.getoutput ( cmd )

    def mkffmpeg ( self ):
        framerate = 10
        cmd = f"ffmpeg -framerate {framerate} -y -pattern_type glob -i '{self.dirname}/p*.png'   -c:v libx264 -pix_fmt yuv420p {self.outfile}"
        subprocess.getoutput ( cmd )
        print ( f"[moveMaker] finished making {self.outfile}" )

    def checkDictFile ( self ):
        """ see if there are timestamps in the dict file """
        f = open ( self.dictfile )
        lines = f.readlines()
        f.close()
        data = eval("\n".join(lines[1:]))
        firstName, firstValues = list ( data.items() )[0]
        if not "timestamp" in firstValues:
            print ( f"[movieMaker] dictionary file {self.dictfile} seems to have no timestamps! Aborting." ) 
            sys.exit()
        return True

    def mkpics ( self ):
        for year in range(2016,2024):
            months = range(1,13)
            if year == 2016:
                months = range(9,13)
            if year == 2023:
                months = range(1,5)
            for month in months:
                date = f"{year}/{month:02d}/01"
                filename = f"{self.dirname}/p{date.replace('/','_')}.png"
                import sys
                sys.path.insert(0,"../../")
                from protomodels.plotting import plotDBDict
                topos = None
                # topos = "electroweakinos"
                poptions = { "topologies": topos }
                poptions["roughviz"]=False
                poptions["dictfile"] = self.dictfile
                options = {}
                options['ylabel']='# signal regions'
                options['plot_averages']= False
                poptions['verbose']=3
                options['plotStats']= False
                options['alwayslegend']= True
                options['yrange']=(0,300)
                poptions["options"] = options
                poptions["outfile"] = filename
                poptions["title"] = f"SModelS database: {date}"
                poptions["before"] = date
                plotter = plotDBDict.Plotter ( poptions )

    def create( self ):
        self.mkdir()
        self.mkpics()
        self.mkffmpeg()
        self.cleanUp()

    def cleanUp ( self ):
        """ remove all the individual pics """
        if os.path.exists ( self.dirname ):
            os.unlink ( self.dirname )

if __name__ == "__main__":
    import argparse
    argparser = argparse.ArgumentParser(description='Makes a video clip that shows the evolution of the meta statistics over time')
    argparser.add_argument ( '-d', '--dictfile',
                         help='the input dictfile [../timestamps.dict]',
                         type=str, default='../timestamps.dict' )
    args = argparser.parse_args()
    maker = MovieMaker( args.dictfile )
    maker.mkdir()
    maker.mkpics()
    maker.mkffmpeg()
