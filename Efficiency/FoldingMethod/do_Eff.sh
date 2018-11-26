#!/bin/bash

for((i=6; i<13; i++))
do
	root -l -b -q 'ScaledMCeff.C+('$i', 3, 3)'
	root -l -b -q 'ScaledMCeff.C+('$i', 2, 4)'
done
