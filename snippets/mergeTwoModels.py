#!/usr/bin/env python

from typing import Tuple

def merge( models : Tuple ):
    import time
    from ptools.hiscoreTools import mergeTwoModels
    model = mergeTwoModels ( *models )
    print ( "model", model )
    f=open("merged.dict","wt")
    f.write ( "# this model is due to a merge of the ewkino and the hadrons models\n" )
    f.write ( f"# created {time.asctime()}\n" )
    f.write ( str(model)+"\n" )
    f.close()

if __name__ == "__main__":
    import argparse 
    argparser = argparse.ArgumentParser(
        description='merge two protomodels' )
    argparser.add_argument ( '-1', '--model1', type=str,
        help="first protomodel [./pmodel1.dict]", 
        default="./pmodel1.dict" )
    argparser.add_argument ( '-2', '--model2', type=str,
        help="second protomodel [./pmodel2.dict]", 
        default="./pmodel2.dict" )
    args = argparser.parse_args()
    merge( ( args.model1, args.model2 ) )
