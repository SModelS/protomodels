#!/usr/bin/env python3

import subprocess, os

class MovieMaker:
    def __init__ ( self ):
        self.dirname = "pics"

    def mkdir ( self ):
        if os.path.exists ( self.dirname ):
            cmd = f"rm -r {self.dirname}"
            subprocess.getoutput ( cmd )
        cmd = f"mkdir {self.dirname}"
        subprocess.getoutput ( cmd )

    def mkffmpeg ( self ):
        framerate = 20
        outfile = "out.mp4"
        cmd = f"ffmpeg -framerate {framerate} -y -pattern_type glob -i '{self.dirname}/p*.png'   -c:v libx264 -pix_fmt yuv420p {outfile}"
        subprocess.getoutput ( cmd )

    def mkpics ( self ):
        for year in range(2016,2024):
            months = range(1,13)
            if year == 2016:
                months = range(9,13)
            for month in months:
                date = f"{year}/{month:02d}/01"
                filename = f"{self.dirname}/p{date.replace('/','_')}.png"
                import sys
                sys.path.insert(0,"../../")
                from protomodels.plotting import plotDBDict
                topos = None
                # topos = "electroweakinos"
                poptions = { "topologies": topos, "roughviz": False }
                poptions["dictfile"] = "../db222pre1timestamp.dict"
                poptions["options"] = {'ylabel':'# signal regions', 'plot_averages': False, 'plotStats': False }
                poptions["outfile"] = filename
                poptions["before"] = date
                plotter = plotDBDict.Plotter ( poptions )

    def create( self ):
        self.mkdir()
        self.mkpics()
        self.mkffmpeg()

if __name__ == "__main__":
    maker = MovieMaker()
    # maker.create()
    maker.mkffmpeg()
