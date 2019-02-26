#!/bin/bash

#for((i=0; i<13; i++))
#do
#	root -l -b -q 'GetScaleFactor.C+('$i')'
#done

root -l -b -q 'GetScaleFactor.C+(0)'
root -l -b -q 'GetScaleFactor.C+(5)'
