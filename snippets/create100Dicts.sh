#!/bin/sh
#
for i in `seq -f "%03.0f" 100`; do 
	[ -e dicts/$i.dict ] || {
		echo $i;
		./ptools/expResModifier.py -C -d official -o none -s $i; 
		mv $i.dict dicts/;
	};
done
