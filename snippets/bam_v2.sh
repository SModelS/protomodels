#!/bin/sh

./plotting/plotBAM.py -s 13 -e ATLAS --exclude "ATLAS-EXOT-*,ATLAS-SUSY-2018-22-multibin,ATLAS-SUSY-2018-16-hino,ATLAS-SUSY-2018-05-ewk" --rename "{'ATLAS-SUSY-2018-05-strong':'ATLAS-SUSY-2018-05'}" --show

./plotting/plotBAM.py -s 13 -e CMS --exclude "CMS-EXO-*" --show
