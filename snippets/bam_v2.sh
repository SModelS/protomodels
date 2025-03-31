#!/bin/sh

../plotting/plotBAM.py -s 13 -e ATLAS --exclude "ATLAS-EXOT-*,ATLAS-SUSY-2018-22-multibin,ATLAS-SUSY-2018-16-hino,ATLAS-SUSY-2018-05-ewk" --rename "{'ATLAS-SUSY-2018-05-strong':'ATLAS-SUSY-2018-05'}" --show --effmaps_only -t --title "ATLAS-"

# exclude cms-exo-19-* so that cms-exo-20-004 is in
../plotting/plotBAM.py -s 13 -e CMS --exclude "CMS-EXO-19-*" --show --effmaps_only -t --title "CMS-"
 
../plotting/plotBAM.py -s 8 -e ATLAS --show --effmaps_only -t --title "ATLAS-"

../plotting/plotBAM.py -s 8 -e CMS --show --effmaps_only -t --title "CMS-"

../plotting/plotDBDict.py -d ../share/300.dict --show -a '^ATLAS-EXOT-*,^CMS-EXO-19-*,^ATLAS-SUSY-2018-22-multibin,^ATLAS-SUSY-2018-16-hino,^ATLAS-SUSY-2018-05-ewk' -T "protomodels_v2: all" -o "metastats_all.png" --pvalues

../plotting/plotDBDict.py -d ../share/300.dict --show -a '^ATLAS-EXOT-*,^CMS-EXO-19-*,^ATLAS-SUSY-2018-22-multibin,^ATLAS-SUSY-2018-16-hino,^ATLAS-SUSY-2018-05-ewk' -t 'electroweakinos,darkmatter' -T "protomodels_v2: ewkinos + dark matter" -o "metastats_ewkinos.png" --pvalues
