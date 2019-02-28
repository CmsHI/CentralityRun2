#!/bin/bash
#ScaledMCeff.C(Int_t ivar, Int_t Coin, Int_t Th, Int_t MCN, Int_t HFRange, Int_t RangeN)

root -l -b -q 'ScaledMCeff.C+(0, 2, 4, 0, 3, 3)'
root -l -b -q 'ScaledMCeff.C+(0, 2, 4, 0, 4, 4)'
root -l -b -q 'ScaledMCeff.C+(0, 2, 4, 1, 2, 2)'
root -l -b -q 'ScaledMCeff.C+(0, 2, 4, 1, 4, 4)'
root -l -b -q 'ScaledMCeff.C+(0, 2, 4, 2, 4, 4)'
root -l -b -q 'ScaledMCeff.C+(0, 2, 4, 2, 5, 5)'
root -l -b -q 'ScaledMCeff.C+(0, 2, 4, 3, 2, 2)'
root -l -b -q 'ScaledMCeff.C+(5, 2, 4, 0, 3, 3)'
root -l -b -q 'ScaledMCeff.C+(5, 2, 4, 0, 4, 4)'
root -l -b -q 'ScaledMCeff.C+(5, 2, 4, 1, 2, 2)'
root -l -b -q 'ScaledMCeff.C+(5, 2, 4, 1, 4, 4)'
root -l -b -q 'ScaledMCeff.C+(5, 2, 4, 2, 4, 4)'
root -l -b -q 'ScaledMCeff.C+(5, 2, 4, 2, 5, 5)'
root -l -b -q 'ScaledMCeff.C+(5, 2, 4, 3, 2, 2)'
