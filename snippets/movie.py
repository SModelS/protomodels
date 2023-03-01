#!/usr/bin/env python3

import subprocess, os

class MovieMaker:
    def __init__ ( self ):
        self.dirname = "pics"

    def mkdir ( self ):
        if os.path.exists ( self.dirname ):
            cmd = f"rm -r {self.dirname}"
        cmd = f"mkdir {self.dirname}"
        subprocess.getoutput ( cmd )

    def mkffmpeg ( self ):
        cmd = f"ffmpeg -framerate 5 -y -pattern_type glob -i '{self.dirname}/p*.png'   -c:v libx264 -pix_fmt yuv420p out.mp4"
        subprocess.getoutput ( cmd )

    def mkpics ( self ):
        for year in range(2016,2024):
            months = range(1,13)
            if year == 2016:
                months = range(9,13)
            for month in months:
                date = f"{year}/{month:02d}/01"
                filename = f"{self.dirname}/p{date.replace('/','_')}.png"
                cmd = f'../plotting/plotDBDict.py -d ../db222pre1timestamp.dict --before "{date}" -o {filename} -t electroweakinos'
                subprocess.getoutput ( cmd )

    def create( self ):
        self.mkdir()
        self.mkpics()
        self.mkffmpeg()

if __name__ == "__main__":
    maker = MovieMaker()
    maker.create()
