#!/bin/sh

for i in `seq 1 49`; do ./expResModifier.py -R $RUNDIR.fake${i} -d $RUNDIR/original.pcl -s fake${i} --remove_orig --nofastlim --onlyvalidated --nosuperseded --symlink -o $RUNDIR.fake${i}/filtered${i}.pcl ; done
