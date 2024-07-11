#!/bin/sh
#
for i in `seq -f "%03.0f" 100`; do 
	j="${i}f07";
	[ -e dicts/$j.dict ] || {
		echo $j;
#./ptools/expResModifier.py -C -d official -o none -s ${i}; 
		./ptools/expResModifier.py -f 0.7 -C -d official -o none -s ${j};
		mv $j.dict dicts/;
	};
done
