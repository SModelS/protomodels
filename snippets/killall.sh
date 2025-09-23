#!/bin/sh

for i in `cat jobs`; do
	echo $i;
	scancel $i;
done

mv jobs jobs.previous
