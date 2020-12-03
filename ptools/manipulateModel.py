#!/usr/bin/env python3

import sys, os

class ModelManipulator:
    def __init__ ( self, inputfile, outputfile ):
        self.inputfile = inputfile
        self.outputfile = outputfile
        self.read()

    def pprint ( self, *args ):
        """ logging """
        print ( "[manipulateModel] %s" % ( " ".join(map(str,args))) )

    def read ( self ):
        if not os.path.exists ( self.inputfile ):
            self.pprint ( f"error, {self.inputfile} not found." )
            sys.exit()
        with open ( self.inputfile, "rt" ) as h:
            txt = h.read()
            self.model = eval ( txt )

    def write ( self ):
        with open ( self.outputfile, "wt" ) as h:
            h.write ( "%s\n" % self.model )
            h.close()

    def multiply ( self, factor ):
        """ multiply ssms with <factor> """
        for pids, v in self.model["ssmultipliers"].items():
            self.model["ssmultipliers"][pids]=v*factor
        comment = "multiplied with %.2f" % factor
        if "comment" in self.model:
            comment = self.model["comment"] + ", " + comment
        self.model["comment"] = comment

if __name__ == "__main__":
    import argparse
    argparser = argparse.ArgumentParser(
                        description='simple tool to manipulate a protomodel in a file' )
    argparser.add_argument ( '-i', '--input',
            help='inputfile [pmodel.py]',
            type=str, default="pmodel.py" )
    argparser.add_argument ( '-o', '--output',
            help='outputfile [out.py]',
            type=str, default="out.py" )
    argparser.add_argument ( '-m', '--multiply',
            help='multiply all ssms with this value [1.0]',
            type=float, default=1.0 )
    args = argparser.parse_args()

    manipulator = ModelManipulator ( args.input, args.output )
    if abs ( args.multiply - 1. ) > 1e-5:
        manipulator.multiply ( args.multiply )
    manipulator.write()
