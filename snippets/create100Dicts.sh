#!/bin/sh
#
for i in `seq -f "%03.0f" 100`; do 
	[ -e $i.dict ] || {
		echo $i;
		./ptools/expResModifier.py -C -d official -o none -s $i; 
	};
done
