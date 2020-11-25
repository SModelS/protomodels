#!/bin/sh

for i in `seq 1 19`; do ./expResModifier.py -R $RUNDIR.signal${i} -d $RUNDIR/original.pcl -s signal${i} -P $RUNDIR/pmodel.py --remove_orig --nofastlim --onlyvalidated --nosuperseded --symlink --seed 100$i -o $RUNDIR.signal${i}/filtered${i}.pcl ; done
